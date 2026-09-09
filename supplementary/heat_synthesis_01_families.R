# R 4.5.2; base R. Run from publication_scripts/.
args <- commandArgs(TRUE)
out <- if (length(args)) args[1] else 'results/heat/synthesis_revised'
stopifnot(!dir.exists(out))
contrasts <- read.delim('config/heat_contrasts.tsv', quote = '', check.names = FALSE)
models <- read.delim('config/heat_models.tsv', quote = '', check.names = FALSE)
genes <- read.delim('config/gene_descriptions.tsv', quote = '', check.names = FALSE)$gene_id
stopifnot(!anyDuplicated(genes), !anyDuplicated(contrasts$contrast_id), !anyDuplicated(models$model_id))
heat <- contrasts[contrasts$comparison_family == 'heat' & contrasts$contrast_id != 'B_average_heat', ]
heat$family <- models$publication_family_id[match(heat$model_id, models$model_id)]
stopifnot(nrow(heat) > 0, !anyNA(heat$family))
path <- 'results/heat/inference/ordinary_results.tsv'
header <- names(read.delim(path, nrows = 0, check.names = FALSE))
keep <- c('gene_id', 'contrast_id', 'model_id', 'publication_family_id', 'analysis_tier', 'multiplicity_status', 'log2FoldChange', 'lfcSE', 'pvalue', 'ordinary_valid', 'beta_converged')
stopifnot(all(keep %in% header))
classes <- ifelse(header %in% c('log2FoldChange', 'lfcSE', 'pvalue'), 'numeric', ifelse(header %in% keep, 'character', 'NULL'))
rows <- read.delim(path, quote = '', check.names = FALSE, colClasses = classes)
rows <- rows[rows$contrast_id %in% heat$contrast_id, ]
stopifnot(setequal(rows$gene_id, genes), !anyDuplicated(rows[c('gene_id', 'contrast_id')]))
stopifnot(nrow(rows) == length(genes) * nrow(heat))
stopifnot(all(rows$ordinary_valid %in% c('True', 'TRUE', 'False', 'FALSE')))
stopifnot(all(rows$analysis_tier == 'conditional_exploratory'), all(rows$multiplicity_status == 'final_declared_conditional_family'))
stopifnot(all(rows$model_id == heat$model_id[match(rows$contrast_id, heat$contrast_id)]))
rows$family <- heat$family[match(rows$contrast_id, heat$contrast_id)]
stopifnot(all(rows$publication_family_id == rows$family))
rows <- rows[rows$ordinary_valid %in% c('True', 'TRUE'), ]
stopifnot(all(is.finite(rows$log2FoldChange)), all(is.finite(rows$lfcSE) & rows$lfcSE > 0))
stopifnot(all(rows$beta_converged %in% c('True', 'TRUE')), all(is.finite(rows$pvalue) & rows$pvalue >= 0 & rows$pvalue <= 1))
rows$w <- 1 / rows$lfcSE^2
fam <- aggregate(cbind(wsum = w, wx = w * log2FoldChange, wse = w * lfcSE, k = 1) ~ gene_id + family, rows, sum)
families <- unique(heat$family)
grid <- expand.grid(gene_id = genes, family = families, stringsAsFactors = FALSE)
fam <- merge(grid, fam, by = c('gene_id', 'family'), all.x = TRUE, sort = TRUE)
expected <- table(heat$family)
fam$expected_contrasts <- as.integer(expected[fam$family])
fam$k[is.na(fam$k)] <- 0
fam$complete <- fam$k == fam$expected_contrasts
fam$estimate <- ifelse(fam$complete, fam$wx / fam$wsum, NA_real_)
fam$se_independence <- ifelse(fam$complete, sqrt(1 / fam$wsum), NA_real_)
# For positive normalized weights a, Var(sum(a*y)) <= (sum(a*SE))^2 for any correlations.
fam$se <- ifelse(fam$complete, fam$wse / fam$wsum, NA_real_)
fam$z <- fam$estimate / fam$se
fam$p_up <- ifelse(fam$complete, pnorm(fam$z, lower.tail = FALSE), 1)
fam$p_down <- ifelse(fam$complete, pnorm(fam$z), 1)
fam$pvalue <- pmin(1, 2 * pmin(fam$p_up, fam$p_down))
fam$q_within_family <- ave(fam$pvalue, fam$family, FUN = function(p) p.adjust(p, 'BH'))
fam$status <- ifelse(fam$complete, 'complete_covariance_upper_bound', 'incomplete_no_fixed_context_estimate')
fam$analysis_tier <- 'conditional_exploratory'
stopifnot(all(fam$se[fam$complete] >= fam$se_independence[fam$complete] - 1e-12))
dir.create(out, recursive = TRUE)
write.table(fam, file.path(out, 'family_effects.tsv'), sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
write.table(heat, file.path(out, 'included_contrasts.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
design <- read.delim('config/heat_study_design.tsv', quote = '', check.names = FALSE)
metadata <- unique(design[design$model_id %in% heat$model_id, c('publication_family_id', 'model_id', 'tissue', 'heat_regime')])
write.table(metadata, file.path(out, 'context_metadata.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, 'family_sessionInfo.txt'))
