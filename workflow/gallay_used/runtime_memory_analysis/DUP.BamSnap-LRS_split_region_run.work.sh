#!/bin/bash
mkdir -p regions plots logs

# Split BED into individual regions
# rm -f regions/DUP_*.bed

awk -v dir="regions" '
BEGIN {FS=OFS="\t"}
NF && $1 !~ /^#/ {
    n++
    file=sprintf("%s/DUP_%02d.bed", dir, n)
    print > file
    close(file)
}
' ../DUP_sample_50.bed

# Run each region separately
for region_bed in regions/DUP_*.bed
do
    id=$(basename "$region_bed" .bed)

    echo "Running ${id} ..."

    time -v \
        -o "logs/${id}.time.txt" \
        python bamsnap-lrs dna \
        --bam SRR16005301.ont.hg38.sorted.bam \
              SRR8955953.pacbio.hg38.sorted.bam \
        --regions "$region_bed" \
        --out-prefix "plots/${id}.png" \
        --fa hg38.fa \
        --show-axis \
        --show-coverage \
        --padding 2000 \
        --show-supp \
        > "logs/${id}.log" 2>&1

    if [ $? -eq 0 ]; then
        echo "${id}: SUCCESS"
    else
        echo "${id}: FAILED"
    fi

done

echo "All DUP regions completed."


{
    printf "sv_ID\tRuntime_s\tMaximum_RSS_GiB\n"

    for f in logs/DUP_*.time.txt
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
} > DUP_runtime_memory.tsv

cat DUP_runtime_memory.tsv
