# Barley heat-response analysis

**Work in progress for a research article.** Core workflow for read processing, differential expression, cross-study synthesis and candidate heat-response gene tables.

## Workflow

Use Linux or WSL2. Run commands from the repository root and inspect outputs between stages. Input locations and study definitions are in `config/`. Set `GENOME`, `GTF`, `STAR_INDEX`, `STAR_TMP` and `ANNOTATION_BED` before reference preparation and alignment. Run one alignment at a time.

```bash
micromamba create -f environment.yml
micromamba activate barley-heat
```

| Stage | Scripts, in execution order |
| --- | --- |
| Reference inputs | `main/heat_00_resources.py FILE`, then `bash main/heat_00_reference.sh` after decompressing the FASTA/GTF |
| Read inputs | `main/heat_00_download.py RUN`, then `main/heat_01_integrity.py RUN` |
| Read QC and alignment | `bash main/heat_02_fastqc.sh RUN`, `bash main/heat_03_trim.sh RUN`, `bash main/heat_04_align.sh RUN` |
| Strand review and counts | `bash supplementary/heat_alignment_qc.sh RUN`, then `bash main/heat_05_count.sh MODEL` |
| Sample QC and differential expression | `supplementary/heat_sample_qc.R MODEL`, `main/heat_06_deseq2.R MODEL`, `main/heat_07_shrink.R MODEL` |
| Collect results and pathways | `tables/heat_01_results.R`, then `supplementary/heat_pathways.R` |
| Cross-family synthesis | `supplementary/heat_synthesis_01_families.R`, then `supplementary/heat_synthesis_02_models.R` |
| Atlas and candidate tables | `supplementary/01_atlas_expression.py config/atlas_panel_samples.tsv results/atlas_panel`, `supplementary/02_gene_mapping.py`, `tables/heat_07_synthesis_revised.py` |

Run `.py` files with `python` and `.R` files with `Rscript --vanilla`. `RUN` and `MODEL` refer to entries in the sample and model tables. Results are written under `results/`. Download the atlas resources listed in `config/resource_downloads.tsv` before the atlas stage. R package versions are recorded in `renv.lock`.

## Permission

The analysis scripts and workflows may be copied or used only by editors and reviewers of journals considering the associated article, solely for confidential editorial assessment and peer review. All other copying, reuse, modification or distribution requires prior written permission from the authors. See [LICENSE](LICENSE). Third-party data and software retain their own terms.
