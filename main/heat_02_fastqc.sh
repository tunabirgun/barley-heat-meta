#!/usr/bin/env bash
set -euo pipefail
run_filter=${1:?Supply one run accession}
matched=false
mkdir -p results/heat/fastqc/raw
while IFS=$'\t' read -r model sample run genotype condition replicate tissue time read1 read2; do
    [[ "$run" == "$run_filter" ]] || continue
    matched=true
    test -s "results/heat/integrity/$run.tsv"
    fastqc --threads 2 --outdir results/heat/fastqc/raw "$read1" "$read2"
done < <(tail -n +2 config/heat_samples.tsv)
$matched
