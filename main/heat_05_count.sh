#!/usr/bin/env bash
set -euo pipefail
model_filter=${1:?Supply one model ID}
: "${GTF:?}"
mkdir -p results/heat/counts
bams=()
strands=()
while IFS=$'\t' read -r model sample run genotype condition replicate tissue time read1 read2; do
    [[ "$model" == "$model_filter" ]] || continue
    strand=$(awk -F '\t' -v run="$run" '$2==run && $NF=="passed" {print $3}' config/heat_strands.tsv)
    [[ "$strand" =~ ^[012]$ ]] || { echo "Missing reviewed strand: $run" >&2; exit 1; }
    bam="results/heat/aligned/${run}_Aligned.sortedByCoord.out.bam"
    test -s "$bam.csi"
    samtools quickcheck -v "$bam"
    bams+=("$bam")
    strands+=("$strand")
done < <(tail -n +2 config/heat_samples.tsv)
[[ ${#bams[@]} -ge 4 ]]
strand_list=$(IFS=,; echo "${strands[*]}")
test ! -e "results/heat/counts/$model_filter.txt"
featureCounts -T 20 -p --countReadPairs -s "$strand_list" -Q 10 -t exon -g gene_id -a "$GTF" \
    -o "results/heat/counts/$model_filter.txt" "${bams[@]}"
