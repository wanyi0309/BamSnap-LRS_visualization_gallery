mkdir -p regions plots logs

# Split BED into individual regions
rm -f regions/Highlight_*.bed

awk -v dir="regions" '
BEGIN {FS=OFS="\t"}
NF && $1 !~ /^#/ {
    n++
    file=sprintf("%s/Highlight_%02d.bed", dir, n)
    print > file
    close(file)
}
' ../random_50_regions.bed

# Run each region separately
for region_bed in regions/Highlight_*.bed
do
    id=$(basename "$region_bed" .bed)

    echo "Running ${id} ..."

    if time -v \
        -o "logs/${id}.time.txt" \
        python bamsnap-lrs highlight \
        --highlight-vcf ../HG002.phased_het_snps.vcf.gz \
        --bam HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam \
        --regions "$region_bed" \
        --out-prefix "plots/${id}.png" \
        --fa hg38.fa \
        --no-hap-filter \
        --padding 0 \
        --mapq 0 \
        > "logs/${id}.log" 2>&1
    then
        echo "${id}: SUCCESS"
    else
        echo "${id}: FAILED"
    fi

done

echo "All Highlight regions completed."


{
    printf "sv_ID\tRuntime_s\tMaximum_RSS_GiB\n"

    for f in logs/Highlight_*.time.txt
    do
        region=$(basename "$f" .time.txt)
        region=${region//_/}

        awk -v region="$region" '
        /Elapsed \(wall clock\) time/ {
            t=$NF
            n=split(t, a, ":")

            if (n == 2)
                runtime = a[1]*60 + a[2]
            else if (n == 3)
                runtime = a[1]*3600 + a[2]*60 + a[3]

            has_time=1
        }

        /Maximum resident set size \(kbytes\)/ {
            memory=$NF/(1024*1024)
            has_memory=1
        }

        END {
            if (has_time && has_memory)
                printf "%s\t%.3f\t%.3f\n", region, runtime, memory
            else
                printf "%s\tNA\tNA\n", region
        }
        ' "$f"

    done
} > Highlight_runtime_memory.tsv
