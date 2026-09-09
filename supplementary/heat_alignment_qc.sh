#!/usr/bin/env bash
set -euo pipefail
run=${1:?Supply one run accession}
: "${ANNOTATION_BED:?}"
bam="results/heat/aligned/${run}_Aligned.sortedByCoord.out.bam"
test -s "$bam.csi"
if [[ ! -e "$bam.bai" ]]; then
    ln -s "$(basename "$bam").csi" "$bam.bai"
fi
[[ "$(readlink -f "$bam.bai")" == "$(readlink -f "$bam.csi")" ]]
mkdir -p results/heat/rseqc
infer_experiment.py -i "$bam" -r "$ANNOTATION_BED" -s 200000 -q 30 > "results/heat/rseqc/$run.infer_experiment.txt" 2>&1
python supplementary/heat_alignment_qc_distribution.py -i "$bam" -r "$ANNOTATION_BED" > "results/heat/rseqc/$run.read_distribution.txt" 2>&1
geneBody_coverage.py -i "$bam" -r "$ANNOTATION_BED" -o "results/heat/rseqc/$run" > "results/heat/rseqc/$run.coverage.log" 2>&1
