# R 4.5.2; DESeq2 1.50.2. One serial worker.
source('supplementary/heat_sample_qc.R')

heat_vector <- function(text,results_names) {
  x <- unlist(jsonlite::fromJSON(text),use.names=TRUE)
  if (!is.numeric(x) || is.null(names(x)) || anyDuplicated(names(x)) ||
      any(!is.finite(x)) || !all(names(x) %in% results_names) || !any(x!=0) || abs(sum(x))>1e-10)
    stop('Invalid named, nonzero, sum-zero contrast vector')
  v <- setNames(numeric(length(results_names)),results_names)
  v[names(x)] <- x
  v
}

heat_ordinary <- function(res,dds,genes,keep,model,contrast,nmin) {
  x <- as.data.frame(res)
  conv <- S4Vectors::mcols(dds)$betaConv
  if (is.null(conv) || length(conv)!=nrow(x)) conv <- rep(NA,nrow(x))
  valid <- conv %in% TRUE & is.finite(x$log2FoldChange) & is.finite(x$lfcSE) & x$lfcSE>0 &
    is.finite(x$stat) & is.finite(x$pvalue) & x$pvalue>=0 & x$pvalue<=1
  reason <- ifelse(valid,NA_character_,ifelse(is.na(conv),'coefficient_convergence_unknown',
    ifelse(!conv,'coefficient_nonconvergence',ifelse(!is.finite(x$lfcSE) | x$lfcSE<=0,
      'invalid_standard_error',ifelse(is.na(x$pvalue),'outlier_filter_or_invalid_wald_test','invalid_ordinary_test')))))
  cooks <- S4Vectors::mcols(dds)$maxCooks
  if (is.null(cooks)) cooks <- rep(NA_real_,nrow(x))
  z <- data.frame(gene_id=rownames(x),baseMean=x$baseMean,log2FoldChange=x$log2FoldChange,
    lfcSE=x$lfcSE,wald_stat=ifelse(valid,x$stat,NA_real_),pvalue=ifelse(valid,x$pvalue,NA_real_),
    native_padj=ifelse(valid,x$padj,NA_real_),ordinary_valid=valid,ordinary_na_reason=reason,
    beta_converged=conv,max_cooks=cooks)
  z <- z[match(genes,z$gene_id),,drop=FALSE]
  z$gene_id <- genes
  z$ordinary_valid[!keep] <- FALSE
  z$ordinary_na_reason[!keep] <- 'prefiltered_low_counts'
  z$test_eligible <- keep
  z$conditional_ci_low <- ifelse(z$ordinary_valid,z$log2FoldChange-qnorm(.975)*z$lfcSE,NA_real_)
  z$conditional_ci_high <- ifelse(z$ordinary_valid,z$log2FoldChange+qnorm(.975)*z$lfcSE,NA_real_)
  data.frame(model[rep(1,nrow(z)),c('model_id','publication_family_id','analysis_tier','contrast_role','uncertainty')],
    contrast[rep(1,nrow(z)),c('contrast_id','comparison_family','contrast_label')],
    z,smallest_cell_n=nmin,row.names=NULL,check.names=FALSE)
}

heat_deseq2 <- function(model_id) {
  models <- read.delim('config/heat_models.tsv',check.names=FALSE)
  contrasts <- read.delim('config/heat_contrasts.tsv',check.names=FALSE,quote='')
  admission <- read.delim('config/heat_admission.tsv',check.names=FALSE)
  stopifnot(length(model_id)==1L,model_id %in% models$model_id,
    !anyDuplicated(contrasts$contrast_id),setequal(contrasts$model_id,models$model_id),
    all(contrasts$comparison_family %in% c('heat','baseline','interaction','heatgenotype')),
    !anyDuplicated(admission$model_id))
  admission <- admission[admission$model_id==model_id,,drop=FALSE]
  stopifnot(nrow(admission)==1L,
    all(admission$analysis_tier=='conditional_exploratory'),all(admission$qc_review=='pass'),
    !anyNA(admission$review_note),all(nzchar(trimws(admission$review_note))))
  input <- heat_input(model_id)
  qc <- file.path('results/heat/sample_qc',model_id)
  stopifnot(all(file.exists(file.path(qc,c('pca.tsv','sample_correlations.tsv','sample_metrics.tsv','prefilter.tsv')))))
  prior <- read.delim(file.path(qc,'prefilter.tsv'))
  metrics <- read.delim(file.path(qc,'sample_metrics.tsv'))
  stopifnot(identical(prior$gene_id,rownames(input$counts)),identical(prior$prefilter_pass,input$keep),
    identical(metrics$run,colnames(input$counts)),
    identical(as.numeric(metrics$library_assigned_fragments),as.numeric(colSums(input$counts))))
  cr <- contrasts[contrasts$model_id==model_id,,drop=FALSE]
  for (vector in cr$coefficient_vector)
    heat_vector(vector,colnames(model.matrix(~0+cell,input$samples)))
  out <- 'results/heat/inference'
  if (dir.exists(file.path(out,'models',model_id)))
    stop('Model output already exists; preserve and review before rerunning')
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  ordinary <- vectors <- list()
  model_out <- file.path(out,'models',model_id)
  dir.create(model_out,recursive=TRUE)
  dds <- DESeqDataSetFromMatrix(input$counts[input$keep,,drop=FALSE],input$samples,~0+cell)
  set.seed(812745)
  dds <- DESeq(dds,quiet=TRUE,parallel=FALSE,BPPARAM=BiocParallel::SerialParam(progressbar=FALSE))
  saveRDS(dds,file.path(model_out,'dds.rds'))
  heat_write(input$samples,file.path(model_out,'samples.tsv'))
  normalized <- sweep(input$counts,2,sizeFactors(dds),'/')
  heat_write(data.frame(gene_id=rownames(normalized),normalized,check.names=FALSE),
    file.path(model_out,'normalized_counts.tsv'))
  heat_write(data.frame(run=names(sizeFactors(dds)),size_factor=unname(sizeFactors(dds))),
    file.path(model_out,'size_factors.tsv'))
  for (i in seq_len(nrow(cr))) {
    v <- heat_vector(cr$coefficient_vector[i],resultsNames(dds))
    res <- results(dds,contrast=v,alpha=.05,parallel=FALSE)
    saveRDS(res,file.path(model_out,paste0(cr$contrast_id[i],'.rds')))
    ordinary[[cr$contrast_id[i]]] <- heat_ordinary(res,dds,rownames(input$counts),input$keep,
      input$model,cr[i,,drop=FALSE],input$smallest_cell_n)
    vectors[[cr$contrast_id[i]]] <- data.frame(model_id,contrast_id=cr$contrast_id[i],
      results_name=names(v),coefficient=unname(v))
  }
  audit <- data.frame(model_id,publication_family_id=input$model$publication_family_id,
    analysis_tier='conditional_exploratory',samples=nrow(input$samples),cells=nlevels(input$samples$cell),
    smallest_cell_n=input$smallest_cell_n,input_genes=nrow(input$counts),prefilter_pass=sum(input$keep),
    converged=sum(S4Vectors::mcols(dds)$betaConv %in% TRUE),model_rank=qr(model.matrix(~0+cell,input$samples))$rank,
    dispersion_fit=attr(dispersionFunction(dds),'fitType'))
  ordinary <- do.call(rbind,ordinary)
  ordinary$project_family_q <- NA_real_
  ordinary$family_n_tests <- NA_integer_
  ordinary$multiplicity_status <- 'pending_all_declared_models'
  heat_write(ordinary,file.path(model_out,'ordinary_results.tsv'))
  heat_write(do.call(rbind,vectors),file.path(model_out,'coefficient_vectors.tsv'))
  heat_write(audit,file.path(model_out,'model_audit.tsv'))
  capture.output(sessionInfo(),file=file.path(model_out,'sessionInfo.txt'))
  invisible(ordinary)
}

if (sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  if (length(args)!=1L) stop('Supply exactly one model_id')
  heat_deseq2(args[1])
}
