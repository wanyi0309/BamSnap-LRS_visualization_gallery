#!/usr/bin/env bash

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BAMSNAP_LRS="${BAMSNAP_LRS:-bamsnap-lrs}"
DESIGN_DIR="${HERE}/design"
OUT_DIR="${HERE}/highlight_synthetic_out"
PLOT_DIR="${HERE}/plots"

mkdir -p "${DESIGN_DIR}" "${OUT_DIR}" "${PLOT_DIR}"
unzip -oq "${HERE}/highlight_synthetic_fixture_design.zip" -d "${DESIGN_DIR}"

python "${DESIGN_DIR}/generate_highlight_synthetic.py" \
  --design-dir "${DESIGN_DIR}" \
  --out-dir "${OUT_DIR}"

# Generate the 11 single-sample (S1) fixture plots.
"${BAMSNAP_LRS}" highlight \
  --highlight-vcf "${OUT_DIR}/highlight_synthetic.vcf.gz" \
  --bam "${OUT_DIR}/highlight_synthetic.bam" \
  --regions "${OUT_DIR}/highlight_synthetic_regions.bed" \
  --out-prefix "${PLOT_DIR}/SYN.png" \
  --fa "${OUT_DIR}/highlight_synthetic.fa" \
  --padding 0 \
  --mapq 0 \
  --no-hap-filter \
  --show-supp \
  --track-title S1

# Re-render F05 with the same BAM interpreted against both VCF samples.
# This replaces the single-sample F05 plot with the two-track S1/S2 example.
"${BAMSNAP_LRS}" highlight \
  --bam "${OUT_DIR}/highlight_synthetic.bam" \
  --bam "${OUT_DIR}/highlight_synthetic.bam" \
  --highlight-vcf "${OUT_DIR}/highlight_synthetic.vcf.gz" \
  --highlight-samples S1 S2 \
  --track-title S1 S2 \
  --pos chrSYN:800-1000 \
  --fa "${OUT_DIR}/highlight_synthetic.fa" \
  --out "${PLOT_DIR}/SYN_chrSYN_800_1000.png"

echo "Synthetic data and plots generated in: ${OUT_DIR} and ${PLOT_DIR}"
