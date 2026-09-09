#!/usr/bin/env bash
set -euo pipefail
run=${1:?Supply one run accession}
: "${STAR_INDEX:?}" "${STAR_TMP:?}"
test -s "results/heat/fastqc/trimmed/${run}_1_fastqc.zip"
test -s "results/heat/fastqc/trimmed/${run}_2_fastqc.zip"
test ! -e "results/heat/aligned/${run}_Aligned.sortedByCoord.out.bam"
mkdir -p results/heat/aligned results/heat/resources "$STAR_TMP"
/usr/bin/time -v -o "results/heat/resources/${run}_STAR.txt" \
    STAR --genomeDir "$STAR_INDEX" --readFilesIn "results/heat/trimmed/${run}_1.fastq.gz" "results/heat/trimmed/${run}_2.fastq.gz" \
    --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --runThreadN 20 --runRNGseed 777 \
    --limitBAMsortRAM 32000000000 --outFilterMultimapNmax 10 --outFilterMismatchNoverReadLmax 1.0 --twopassMode None \
    --outTmpDir "$STAR_TMP/$run" --outFileNamePrefix "results/heat/aligned/${run}_"
samtools quickcheck -v "results/heat/aligned/${run}_Aligned.sortedByCoord.out.bam"
samtools index -@ 4 -c "results/heat/aligned/${run}_Aligned.sortedByCoord.out.bam"
samtools flagstat -@ 4 "results/heat/aligned/${run}_Aligned.sortedByCoord.out.bam" > "results/heat/aligned/$run.flagstat.txt"
awk '$4=="primary" && $5=="mapped" {print $1+$3}' "results/heat/aligned/$run.flagstat.txt" > "results/heat/aligned/$run.primary_mapped.count"
awk '$4=="singletons" {print $1+$3}' "results/heat/aligned/$run.flagstat.txt" > "results/heat/aligned/$run.primary_singleton.count"
test -s "results/heat/aligned/$run.primary_mapped.count"
test -s "results/heat/aligned/$run.primary_singleton.count"
