# R 4.5.2. Complete declared conditional gene-by-contrast families.
heat_collect <- function() {
  models <- read.delim('config/heat_models.tsv',check.names=FALSE)
  contrasts <- read.delim('config/heat_contrasts.tsv',check.names=FALSE,quote='')
  genes <- read.delim('config/gene_descriptions.tsv')$gene_id
  stopifnot(!anyDuplicated(models$model_id),!anyDuplicated(contrasts$contrast_id),!anyDuplicated(genes),
    setequal(contrasts$model_id,models$model_id))
  results <- stabilized <- vectors <- audits <- list()
  for (model_id in models$model_id) {
    path <- file.path('results/heat/inference/models',model_id)
    z <- read.delim(file.path(path,'ordinary_results.tsv'),check.names=FALSE)
    cr <- contrasts$contrast_id[contrasts$model_id==model_id]
    stopifnot(all(z$model_id==model_id),setequal(z$contrast_id,cr),
      all(z$analysis_tier=='conditional_exploratory'),!anyDuplicated(z[c('gene_id','contrast_id')]),
      nrow(z)==length(genes)*length(cr),all(is.na(z$project_family_q)),
      all(z$multiplicity_status=='pending_all_declared_models'))
    for (id in cr) stopifnot(setequal(z$gene_id[z$contrast_id==id],genes))
    expected <- contrasts[match(z$contrast_id,contrasts$contrast_id),,drop=FALSE]
    model <- models[models$model_id==model_id,,drop=FALSE]
    stopifnot(identical(z$comparison_family,expected$comparison_family),
      identical(z$contrast_label,expected$contrast_label),all(z$publication_family_id==model$publication_family_id),
      all(z$contrast_role==model$contrast_role),is.logical(z$ordinary_valid),is.logical(z$beta_converged),
      is.logical(z$test_eligible),!anyNA(z$ordinary_valid),!anyNA(z$test_eligible),
      all(!z$ordinary_valid | (z$test_eligible & z$beta_converged %in% TRUE &
        is.finite(z$log2FoldChange) & is.finite(z$lfcSE) & z$lfcSE>0 & is.finite(z$wald_stat) &
        is.finite(z$pvalue) & z$pvalue>=0 & z$pvalue<=1)),
      all(is.na(z$pvalue[!z$ordinary_valid])),all(is.na(z$wald_stat[!z$ordinary_valid])))
    s <- read.delim(file.path('results/heat/shrink',model_id,'stabilized_display.tsv'),check.names=FALSE)
    stopifnot(!anyDuplicated(s[c('gene_id','contrast_id')]),nrow(s)==nrow(z),
      all(s$analysis_tier=='conditional_exploratory'),all(s$model_id==model_id),
      identical(paste(z$gene_id,z$contrast_id),paste(s$gene_id,s$contrast_id)),
      all(is.na(s$stabilized_log2FoldChange[!z$ordinary_valid])),all(is.na(s$stabilized_lfcSE[!z$ordinary_valid])))
    results[[model_id]] <- z
    stabilized[[model_id]] <- s
    vectors[[model_id]] <- read.delim(file.path(path,'coefficient_vectors.tsv'))
    audits[[model_id]] <- read.delim(file.path(path,'model_audit.tsv'))
  }
  ordinary <- do.call(rbind,results)
  for (family in unique(ordinary$comparison_family)) {
    rows <- ordinary$comparison_family==family
    valid <- rows & ordinary$ordinary_valid %in% TRUE & ordinary$beta_converged %in% TRUE &
      is.finite(ordinary$pvalue) & ordinary$pvalue>=0 & ordinary$pvalue<=1
    denominator <- sum(rows & ordinary$test_eligible)
    ordinary$family_n_tests[rows] <- denominator
    ordinary$project_family_q[valid] <- p.adjust(ordinary$pvalue[valid],'BH',n=denominator)
  }
  ordinary$multiplicity_status <- 'final_declared_conditional_family'
  output <- list('results/heat/inference/ordinary_results.tsv'=ordinary,
    'results/heat/shrink/stabilized_display.tsv'=do.call(rbind,stabilized),
    'results/heat/inference/coefficient_vectors.tsv'=do.call(rbind,vectors),
    'results/heat/inference/model_audit.tsv'=do.call(rbind,audits))
  if (any(file.exists(names(output)))) stop('Collected outputs already exist; preserve before rerunning')
  for (path in names(output)) write.table(output[[path]],path,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
  capture.output(sessionInfo(),file='results/heat/inference/collection_sessionInfo.txt')
  invisible(ordinary)
}

if (sys.nframe()==0L) heat_collect()
