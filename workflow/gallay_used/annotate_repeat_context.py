
import argparse
import csv
import subprocess
import tempfile
from pathlib import Path


def parse_info(info):
    values = {}
    flags = set()

    for item in info.split(";"):
        if not item or item == ".":
            continue

        if "=" in item:
            key, value = item.split("=", 1)
            values[key] = value
        else:
            flags.add(item)

    return values, flags


def main():

    parser = argparse.ArgumentParser(
        description="Annotate SV repeat context using GIAB VCF and RepeatMasker BED."
    )

    parser.add_argument(
        "-v", "--vcf",
        required=True,
        help="Input VCF file (INS or DEL)"
    )

    parser.add_argument(
        "-r", "--repeat-bed",
        required=True,
        help="CHM13 RepeatMasker BED file"
    )

    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output TSV file (default: input filename + .repeat_context.tsv)"
    )

    parser.add_argument(
        "--bedtools",
        default="bedtools",
        help="Path to bedtools executable"
    )

    parser.add_argument(
        "--flank",
        type=int,
        default=1000,
        help="Flanking window size (default: 1000 bp)"
    )

    args = parser.parse_args()

    vcf_file = Path(args.vcf)
    repeat_bed = Path(args.repeat_bed)

    if not vcf_file.exists():
        raise FileNotFoundError(vcf_file)

    if not repeat_bed.exists():
        raise FileNotFoundError(repeat_bed)

    if args.flank <= 0:
        raise ValueError("Flank size must be positive.")

    if args.output:
        output_file = Path(args.output)
    else:
        name = vcf_file.name

        if name.endswith(".vcf.gz"):
            stem = name[:-7]
        elif name.endswith(".vcf"):
            stem = name[:-4]
        else:
            stem = vcf_file.stem

        output_file = vcf_file.with_name(
            stem + ".repeat_context.tsv"
        )

    # --------------------------------------
    # Read VCF
    # --------------------------------------

    records = []
    windows = []
    window_lookup = {}

    with open(vcf_file) as f:

        for line in f:

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 8:
                raise ValueError("Invalid VCF record.")

            chrom = fields[0]
            pos = int(fields[1])
            variant_id = fields[2]
            info = fields[7]

            info_values, info_flags = parse_info(info)

            svtype = info_values.get("SVTYPE")

            if svtype not in ("INS", "DEL"):
                continue

            if "SVLEN" not in info_values:
                raise ValueError(
                    f"Missing SVLEN at {chrom}:{pos}"
                )

            if "," in info_values["SVLEN"]:
                raise ValueError(
                    f"Multiple SVLEN values at {chrom}:{pos}. "
                    "Please use normalized VCF records."
                )

            svlen = abs(int(info_values["SVLEN"]))

            # GIAB VCF RepeatMasker annotation
            rm_annotation = "No"

            for key in ("RM_score", "RM_repeat", "RM_clsfam"):

                value = info_values.get(key)

                if value not in (None, "", "."):
                    rm_annotation = "Yes"
                    break

            # GIAB VCF TRF annotation
            trf_annotation = (
                "Yes"
                if "TRF" in info_flags or "TRF" in info_values
                else "No"
            )

            # Convert VCF coordinates to BED coordinates
            if svtype == "INS":

                sv_start = pos
                sv_end = pos

            else:

                # DEL: exclude the VCF anchor base
                sv_start = pos
                sv_end = pos + svlen

            sample_key = f"SV_{len(records) + 1:03d}"

            record = {
                "sample_key": sample_key,
                "CHROM": chrom,
                "POS": pos,
                "ID": variant_id,
                "SVTYPE": svtype,
                "SVLEN": svlen,
                "sv_start_0based": sv_start,
                "sv_end_0based": sv_end,
                "RM_annotation": rm_annotation,
                "TRF_annotation": trf_annotation
            }

            records.append(record)

            # Generate upstream and downstream windows
            for side, start, end in [
                (
                    "UP",
                    max(0, sv_start - args.flank),
                    sv_start
                ),
                (
                    "DOWN",
                    sv_end,
                    sv_end + args.flank
                )
            ]:

                window_id = f"{sample_key}_{side}"

                windows.append([
                    chrom,
                    start,
                    end,
                    window_id
                ])

                window_lookup[window_id] = (
                    sample_key,
                    side
                )

    if not records:
        raise ValueError("No INS or DEL records found in input VCF.")

    print(f"Total SV records: {len(records)}")
    print(f"INS: {sum(r['SVTYPE'] == 'INS' for r in records)}")
    print(f"DEL: {sum(r['SVTYPE'] == 'DEL' for r in records)}")

    # --------------------------------------
    # Run bedtools coverage
    # --------------------------------------

    coverage_results = {}

    with tempfile.TemporaryDirectory() as tmp:

        tmp = Path(tmp)

        window_file = tmp / "windows.bed"
        coverage_file = tmp / "coverage.tsv"

        with open(window_file, "w") as f:

            for chrom, start, end, window_id in windows:

                f.write(
                    f"{chrom}\t{start}\t{end}\t{window_id}\n"
                )

        with open(coverage_file, "w") as out:

            subprocess.run(
                [
                    args.bedtools,
                    "coverage",
                    "-a", str(window_file),
                    "-b", str(repeat_bed)
                ],
                stdout=out,
                check=True
            )

        with open(coverage_file) as f:

            for line in f:

                fields = line.rstrip("\n").split("\t")

                window_id = fields[3]

                covered_bp = int(fields[5])
                window_length = int(fields[6])
                fraction = float(fields[7])

                sample_key, side = window_lookup[window_id]

                coverage_results[(sample_key, side)] = {
                    "annotation": (
                        "Yes" if covered_bp > 0 else "No"
                    ),
                    "covered_bp": covered_bp,
                    "fraction": fraction
                }

    # --------------------------------------
    # Write output TSV
    # --------------------------------------

    columns = [
        "sample_key",
        "CHROM",
        "POS",
        "ID",
        "SVTYPE",
        "SVLEN",
        "sv_start_0based",
        "sv_end_0based",
        "RM_annotation",
        "TRF_annotation",
        "Repeat_upstream_1kb",
        "Repeat_upstream_bp",
        "Repeat_upstream_fraction",
        "Repeat_downstream_1kb",
        "Repeat_downstream_bp",
        "Repeat_downstream_fraction"
    ]

    with open(output_file, "w", newline="") as f:

        writer = csv.DictWriter(
            f,
            fieldnames=columns,
            delimiter="\t"
        )

        writer.writeheader()

        for record in records:

            key = record["sample_key"]

            up = coverage_results[(key, "UP")]
            down = coverage_results[(key, "DOWN")]

            row = dict(record)

            row.update({
                "Repeat_upstream_1kb": up["annotation"],
                "Repeat_upstream_bp": up["covered_bp"],
                "Repeat_upstream_fraction": round(
                    up["fraction"], 4
                ),

                "Repeat_downstream_1kb": down["annotation"],
                "Repeat_downstream_bp": down["covered_bp"],
                "Repeat_downstream_fraction": round(
                    down["fraction"], 4
                )
            })

            writer.writerow(row)

    # --------------------------------------
    # Summary
    # --------------------------------------

    print("\n========== Summary ==========")

    for field in ("RM_annotation", "TRF_annotation"):

        count = sum(
            r[field] == "Yes"
            for r in records
        )

        print(f"{field}: {count}/{len(records)}")

    for side in ("UP", "DOWN"):

        count = sum(
            coverage_results[(r["sample_key"], side)]["annotation"] == "Yes"
            for r in records
        )

        print(f"Repeat {side} 1kb: {count}/{len(records)}")

    print("\nOutput:", output_file)
    print("Completed!")


if __name__ == "__main__":
    main()
