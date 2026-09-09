#!/usr/bin/env bash
set -euo pipefail
run_filter=${1:?Supply one run accession}
matched=false
mkdir -p results/heat/trimmed results/heat/fastp results/heat/fastqc/trimmed
while IFS=$'\t' read -r model sample run genotype condition replicate tissue time read1 read2; do
    [[ "$run" == "$run_filter" ]] || continue
    matched=true
    test -s "results/heat/fastqc/raw/${run}_1_fastqc.zip"
    test -s "results/heat/fastqc/raw/${run}_2_fastqc.zip"
    test ! -e "results/heat/trimmed/${run}_1.fastq.gz"
    poly_g=--disable_trim_poly_g
    [[ "$model" == D_* || "$model" == F_* ]] && poly_g=--trim_poly_g
    [[ -f config/heat_polyg_models.txt ]] && grep -qx "$model" config/heat_polyg_models.txt && poly_g=--trim_poly_g
    fastp -i "$read1" -I "$read2" -o "results/heat/trimmed/${run}_1.fastq.gz" -O "results/heat/trimmed/${run}_2.fastq.gz" \
        --detect_adapter_for_pe --qualified_quality_phred 15 --unqualified_percent_limit 40 --length_required 36 \
        --thread 20 "$poly_g" --json "results/heat/fastp/$run.json" --html "results/heat/fastp/$run.html"
    fastqc --threads 2 --outdir results/heat/fastqc/trimmed "results/heat/trimmed/${run}_1.fastq.gz" "results/heat/trimmed/${run}_2.fastq.gz"
done < <(tail -n +2 config/heat_samples.tsv)
$matched
