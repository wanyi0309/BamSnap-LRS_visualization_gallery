#!/usr/bin/env python3
"""
Calculate BamSnap-LRS statistics directly from the same inputs as its CLI.

Example:
    # DNA
    python bamsnap_benchmark_stats.py dna \
        --bam sample.bam --regions variants.vcf \
        --fa reference.fa --padding 150 \
        --out-prefix results/dna.png \
        --bamsnap-src /path/to/BamSnap-LRS

    # RNA
    python bamsnap_benchmark_stats.py rna \
        --bam sample.rna.bam --regions genes.bed \
        --fa reference.fa --padding 1000 \
        --out-prefix results/rna.png \
        --bamsnap-src /path/to/BamSnap-LRS

    # Highlight (DNA by default; use --mode rna for RNA alignments)
    python bamsnap_benchmark_stats.py highlight \
        --bam sample.bam --regions variants.bed \
        --highlight-vcf phased_snps.vcf \
        --fa reference.fa --padding 1000 \
        --out-prefix results/highlight.png \
        --bamsnap-src /path/to/BamSnap-LRS

The script uses BamSnap-LRS's own region parser, read-fetching functions,
Highlight VCF parser, and pileup implementation. In Highlight mode it applies
the same default VCF-overlap read filter as the renderer before calculating
coverage and displayed-record counts.

Output: one row per region x BAM, with mean depth, maximum depth,
and the number of alignment records retained for display.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path
from typing import Dict, List, Optional

OUTPUT_COLUMNS = [
    "region_index",
    "region_id",
    "chrom",
    "plot_start",
    "plot_end",
    "target_start",
    "target_end",
    "bam",
    "mean_depth",
    "maximum_depth",
    "displayed_records",
]


def add_bamsnap_source(source_dir: Optional[str]) -> None:
    """Add a BamSnap-LRS source checkout to sys.path."""
    if not source_dir:
        return

    source_path = Path(source_dir).expanduser().resolve()
    candidates = [
        source_path,
        source_path / "src",
        source_path.parent if source_path.name == "bamsnap_lrs" else None,
    ]

    for candidate in candidates:
        if candidate is not None and (candidate / "bamsnap_lrs").is_dir():
            sys.path.insert(0, str(candidate))
            return

    raise FileNotFoundError(
        f"Cannot find bamsnap_lrs under {source_path}. "
        "Pass the BamSnap-LRS repository root (containing src/bamsnap_lrs), "
        "the src directory, or the bamsnap_lrs package directory."
    )


def import_bamsnap_functions(source_dir: Optional[str]):
    add_bamsnap_source(source_dir)
    try:
        from bamsnap_lrs.regions import parse_regions_file
        from bamsnap_lrs.reader import fetch_reads, fetch_rna_reads
        from bamsnap_lrs.highlight import parse_highlight_vcf
        from bamsnap_lrs.pileup import base_pileup
    except Exception as exc:
        raise ImportError(
            "Could not import BamSnap-LRS modules. Check --bamsnap-src and "
            "run this script in the Python environment containing BamSnap-LRS "
            "dependencies (including pysam)."
        ) from exc

    return parse_regions_file, fetch_reads, fetch_rna_reads, parse_highlight_vcf, base_pileup


def _open_text(path: Path):
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_region_ids(regions_path: str) -> List[str]:
    """
    Read optional BED/VCF labels in file order.

    Coordinates are still parsed exclusively by BamSnap-LRS's
    parse_regions_file(); these labels are only for making the output easier
    to trace back to the input records.
    """
    path = Path(regions_path)
    is_vcf = str(path).lower().endswith((".vcf", ".vcf.gz"))
    labels: List[str] = []

    try:
        with _open_text(path) as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")

                if is_vcf:
                    if len(fields) < 5:
                        continue
                    try:
                        int(fields[1])
                    except ValueError:
                        continue
                    label = fields[2].strip()
                    if not label or label == ".":
                        label = f"{fields[0]}:{fields[1]}"
                    labels.append(label)
                else:
                    if len(fields) < 3:
                        continue
                    try:
                        start = int(fields[1])
                        end = int(fields[2])
                    except ValueError:
                        continue
                    if start < 0 or end <= start:
                        continue
                    label = fields[3].strip() if len(fields) > 3 else ""
                    labels.append(label or f"{fields[0]}:{start}-{end}")
    except OSError:
        # The authoritative parser below will produce the detailed file error.
        return []

    return labels


def choose_output_path(args) -> Path:
    if args.stats_out:
        return Path(args.stats_out).expanduser().resolve()

    if args.out_prefix:
        prefix = Path(args.out_prefix).expanduser()
        if prefix.suffix.lower() in {".png", ".svg", ".pdf"}:
            # Reuse the plotting command's prefix, but never write statistics
            # into an image file.
            output = prefix.with_name(prefix.stem + ".stats.tsv")
        elif prefix.suffix.lower() in {".tsv", ".txt"}:
            output = prefix
        else:
            output = Path(str(prefix) + ".stats.tsv")
        return output.resolve()

    if args.analysis_mode == "dna":
        return Path("bamsnap_benchmark_stats.tsv").resolve()
    return Path(f"bamsnap_benchmark_stats_{args.analysis_mode}.tsv").resolve()


def filter_highlight_reads(reads, highlight_sites, no_hap_filter: bool = False):
    """Apply the same read-overlap filter used by the Highlight renderer.

    The renderer retains an alignment record when its reference span overlaps
    at least one position from the parsed Highlight VCF. This intentionally
    uses the alignment span (read.start/read.end), not only callable A/C/G/T
    bases; a VCF site inside a deletion or reference skip still satisfies the
    display filter, matching the renderer.
    """
    if no_hap_filter or not highlight_sites:
        return reads

    filter_positions = sorted(site.pos for site in highlight_sites)
    if not filter_positions:
        return reads

    first_pos = filter_positions[0]
    last_pos = filter_positions[-1]

    def covers_any(read):
        if read.end <= first_pos or read.start > last_pos:
            return False
        return any(read.start <= pos < read.end for pos in filter_positions)

    return [read for read in reads if covers_any(read)]


def calculate_depth_stats(reads, start: int, end: int, base_pileup):
    pile = base_pileup(reads, start, end)
    expected_length = end - start

    if len(pile) != expected_length:
        raise RuntimeError(
            f"base_pileup returned {len(pile)} positions for an interval "
            f"of {expected_length} bp ({start}-{end})."
        )

    depths = [int(position.get("depth", 0)) for position in pile]
    if not depths:
        return 0.0, 0

    mean_depth = sum(depths) / len(depths)
    maximum_depth = max(depths)
    return mean_depth, maximum_depth


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate BamSnap-LRS mean depth, maximum depth, and displayed "
            "alignment-record counts for DNA, RNA, or VCF Highlight mode "
            "using the same inputs and read filters as BamSnap-LRS."
        )
    )
    parser.add_argument(
        "analysis_mode",
        nargs="?",
        choices=["dna", "rna", "highlight"],
        default="dna",
        help="Analysis mode: dna, rna, or highlight (default: dna).",
    )
    parser.add_argument(
        "--bam",
        required=True,
        nargs="+",
        action="extend",
        help="BAM/CRAM path(s); may be supplied more than once, as in BamSnap-LRS.",
    )
    parser.add_argument(
        "--regions",
        required=True,
        help="BED/VCF/VCF.GZ file; parsed by BamSnap-LRS regions.parse_regions_file().",
    )
    parser.add_argument("--fa", help="Reference FASTA; used as in the BamSnap-LRS CLI.")
    parser.add_argument(
        "--padding",
        type=int,
        default=None,
        help="Padding passed directly to BamSnap-LRS parse_regions_file().",
    )
    parser.add_argument("--mapq", type=int, default=0, help="Minimum MAPQ (default: 0).")
    parser.add_argument(
        "--show-supp",
        action="store_true",
        help="Include supplementary alignments (default: hidden, matching BamSnap-LRS).",
    )
    parser.add_argument(
        "--show-secondary",
        action="store_true",
        help="Include secondary alignments (default: hidden, matching BamSnap-LRS).",
    )
    parser.add_argument("--use-md", action="store_true")
    parser.add_argument("--use-cs", action="store_true")
    parser.add_argument(
        "--bed",
        default=None,
        help="RNA feature-annotation BED accepted for command compatibility; not used in read statistics.",
    )
    parser.add_argument(
        "--gff",
        default=None,
        help="Gene-annotation GFF/GTF accepted for command compatibility; not used in read statistics.",
    )
    parser.add_argument(
        "--highlight-vcf",
        default=None,
        help="VCF/VCF.GZ used by BamSnap-LRS Highlight mode; required when mode=highlight.",
    )
    parser.add_argument(
        "--highlight-samples",
        nargs="+",
        default=None,
        help="Optional VCF sample subset accepted for Highlight CLI compatibility.",
    )
    parser.add_argument(
        "--mode",
        dest="highlight_data_mode",
        choices=["dna", "rna"],
        default="dna",
        help="Underlying sequencing mode for Highlight: dna (default) or rna.",
    )
    parser.add_argument(
        "--no-hap-filter",
        action="store_true",
        help="Keep reads without overlap to any Highlight VCF site (default: filter them).",
    )
    parser.add_argument(
        "--no-hap-sort",
        action="store_true",
        help="Accepted for Highlight CLI compatibility; sorting does not affect these statistics.",
    )
    parser.add_argument(
        "--hap-layout",
        choices=["packed", "row"],
        default="packed",
        help="Accepted for Highlight CLI compatibility; layout does not affect these statistics.",
    )
    parser.add_argument(
        "--bamsnap-src",
        default=None,
        help="BamSnap-LRS repository root, src/ directory, or bamsnap_lrs package directory.",
    )

    # Output controls. --out-prefix is accepted so the original plotting
    # command can be reused with minimal edits.
    parser.add_argument(
        "--out-prefix",
        default=None,
        help=(
            "Optional plotting-style output prefix. If it ends in .png/.svg/.pdf, "
            "statistics are written beside it as <stem>.stats.tsv."
        ),
    )
    parser.add_argument(
        "--stats-out",
        "--output",
        dest="stats_out",
        default=None,
        help="Explicit statistics TSV output path; takes precedence over --out-prefix.",
    )

    # Rendering-only flags accepted for compatibility with the example command.
    # They do not change these three statistics and are intentionally ignored.
    parser.add_argument("--show-axis", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--show-coverage", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--detail", choices=["low", "mid", "high"], default="mid",
                        help=argparse.SUPPRESS)
    parser.add_argument("--overview-detail", choices=["hide", "show"], default="hide",
                        help=argparse.SUPPRESS)
    parser.add_argument("--width", type=int, default=1200, help=argparse.SUPPRESS)
    parser.add_argument("--read-height", type=int, default=6, help=argparse.SUPPRESS)

    args = parser.parse_args()

    if args.mapq < 0:
        parser.error("--mapq must be >= 0")
    if args.padding is not None and args.padding < 0:
        parser.error("--padding must be >= 0")
    if not args.bam:
        parser.error("At least one --bam file is required")

    for bam_path in args.bam:
        if not Path(bam_path).is_file():
            parser.error(f"BAM/CRAM file not found: {bam_path}")

    if args.fa and not Path(args.fa).is_file():
        parser.error(f"Reference FASTA not found: {args.fa}")

    if any(path.lower().endswith(".cram") for path in args.bam) and not args.fa:
        parser.error("CRAM input requires --fa, matching BamSnap-LRS.")

    regions_path = Path(args.regions)
    if not regions_path.is_file():
        parser.error(f"Regions file not found: {args.regions}")

    if args.analysis_mode == "highlight":
        if not args.highlight_vcf:
            parser.error("Highlight mode requires --highlight-vcf")
        if not Path(args.highlight_vcf).is_file():
            parser.error(f"Highlight VCF not found: {args.highlight_vcf}")

    try:
        (parse_regions_file, fetch_reads, fetch_rna_reads,
         parse_highlight_vcf, base_pileup) = import_bamsnap_functions(args.bamsnap_src)

        # Use the exact BamSnap-LRS region parser and padding implementation.
        regions = parse_regions_file(
            str(regions_path),
            padding=args.padding,
        )
        if not regions:
            raise ValueError(f"No valid regions found in {args.regions}")

        region_ids = read_region_ids(str(regions_path))
        # If malformed input lines make the label count differ from the
        # authoritative BamSnap-LRS parser, fall back to generated labels so
        # labels can never become shifted relative to the parsed coordinates.
        if len(region_ids) != len(regions):
            region_ids = [f"region_{i}" for i in range(1, len(regions) + 1)]
        output_path = choose_output_path(args)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        total = len(regions) * len(args.bam)
        completed = 0
        results: List[Dict[str, str]] = []

        for region_index, region in enumerate(regions, start=1):
            chrom, plot_start, plot_end, target_start, target_end = region

            for bam_path in args.bam:
                completed += 1
                region_id = (
                    region_ids[region_index - 1]
                    if region_index - 1 < len(region_ids)
                    else f"region_{region_index}"
                )

                print(
                    f"[{completed}/{total}] {region_id} | {Path(bam_path).name} | "
                    f"{chrom}:{plot_start}-{plot_end}",
                    file=sys.stderr,
                )

                # Highlight is a wrapper around the DNA/RNA read-fetching
                # pipeline. Its underlying sequencing mode is selected with
                # --mode (dna by default), matching the BamSnap-LRS CLI.
                data_mode = (
                    args.highlight_data_mode
                    if args.analysis_mode == "highlight"
                    else args.analysis_mode
                )
                fetch_function = fetch_rna_reads if data_mode == "rna" else fetch_reads
                reads = fetch_function(
                    bam_path,
                    chrom,
                    plot_start,
                    plot_end,
                    mapq_min=args.mapq,
                    show_supp=args.show_supp,
                    show_secondary=args.show_secondary,
                    use_md=args.use_md,
                    use_cs=args.use_cs,
                    use_ref=bool(args.fa),
                    fa_path=args.fa,
                )

                if args.analysis_mode == "highlight":
                    highlight_sites, highlight_samples = parse_highlight_vcf(
                        args.highlight_vcf,
                        chrom,
                        plot_start,
                        plot_end,
                    )
                    if args.highlight_samples:
                        missing_samples = [
                            sample for sample in args.highlight_samples
                            if sample not in highlight_samples
                        ]
                        if missing_samples:
                            print(
                                f"WARNING: Highlight samples not found in VCF: {missing_samples}",
                                file=sys.stderr,
                            )
                        selected_samples = [
                            sample for sample in args.highlight_samples
                            if sample in highlight_samples
                        ]
                        for site in highlight_sites:
                            site.sample_calls = {
                                name: call for name, call in site.sample_calls.items()
                                if name in selected_samples
                            }

                    before_filter = len(reads)
                    reads = filter_highlight_reads(
                        reads,
                        highlight_sites,
                        no_hap_filter=args.no_hap_filter,
                    )
                    print(
                        f"  Highlight VCF sites: {len(highlight_sites)}; "
                        f"records retained after Highlight filter: {len(reads)}/{before_filter}",
                        file=sys.stderr,
                    )

                # Coverage is calculated after Highlight filtering, matching
                # the renderer, which filters each track before drawing its
                # per-track coverage track and read pileup.
                mean_depth, maximum_depth = calculate_depth_stats(
                    reads,
                    plot_start,
                    plot_end,
                    base_pileup,
                )

                results.append(
                    {
                        "region_index": str(region_index),
                        "region_id": region_id,
                        "chrom": chrom,
                        "plot_start": str(plot_start),
                        "plot_end": str(plot_end),
                        "target_start": str(target_start),
                        "target_end": str(target_end),
                        "bam": bam_path,
                        "mean_depth": f"{mean_depth:.3f}",
                        "maximum_depth": str(maximum_depth),
                        "displayed_records": str(len(reads)),
                    }
                )

        with output_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=OUTPUT_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(results)

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(
        f"\nFinished: {len(regions)} regions x {len(args.bam)} BAM(s) "
        f"= {len(results)} rows.",
        file=sys.stderr,
    )
    print(f"Statistics written to: {output_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
