#!/usr/bin/env bash
set -euo pipefail
: "${GENOME:?Set the decompressed reference FASTA}" "${GTF:?Set the decompressed reference GTF}"
: "${STAR_INDEX:?Set the index directory}" "${ANNOTATION_BED:?Set the output BED12 path}"
test -s "$GENOME"
test -s "$GTF"
test ! -e "$STAR_INDEX/Genome"
test ! -e "$ANNOTATION_BED"
mkdir -p "$STAR_INDEX" "$(dirname "$ANNOTATION_BED")"
STAR --runMode genomeGenerate --genomeDir "$STAR_INDEX" --genomeFastaFiles "$GENOME" \
    --sjdbGTFfile "$GTF" --sjdbOverhang 150 --genomeSAsparseD 2 \
    --runThreadN 4 --limitGenomeGenerateRAM 40000000000
gffread "$GTF" --bed -o "$ANNOTATION_BED"
test -s "$STAR_INDEX/Genome"
test -s "$ANNOTATION_BED"
