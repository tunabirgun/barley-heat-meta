# R 4.5.2; DESeq2 1.50.2; ashr 2.2-63. One serial worker.
suppressPackageStartupMessages(library(DESeq2))

heat_shrink <- function(model_id) {
  models <- read.delim('config/heat_models.tsv')
  stopifnot(length(model_id)==1L,model_id %in% models$model_id)
  ordinary <- read.delim(file.path('results/heat/inference/models',model_id,'ordinary_results.tsv'),check.names=FALSE)
  stopifnot(all(ordinary$analysis_tier=='conditional_exploratory'),
    !anyDuplicated(ordinary[c('model_id','contrast_id','gene_id')]))
  output <- file.path('results/heat/shrink',model_id)
  if (dir.exists(output)) stop('Shrinkage output already exists; preserve and review before rerunning')
  dir.create(output,recursive=TRUE)
  tables <- list()
  for (id in unique(ordinary$contrast_id)) {
    rows <- ordinary[ordinary$contrast_id==id,,drop=FALSE]
    stopifnot(length(unique(rows$model_id))==1L)
    model_dir <- file.path('results/heat/inference/models',rows$model_id[1])
    dds <- readRDS(file.path(model_dir,'dds.rds'))
    res <- readRDS(file.path(model_dir,paste0(id,'.rds')))
    keep <- rows$ordinary_valid %in% TRUE
    valid_genes <- rows$gene_id[keep]
    stopifnot(all(valid_genes %in% rownames(dds)),identical(rownames(dds),rownames(res)))
    z <- rows[,c('gene_id','model_id','publication_family_id','contrast_id','comparison_family',
      'contrast_role','analysis_tier','uncertainty'),drop=FALSE]
    z$stabilized_log2FoldChange <- z$stabilized_lfcSE <- NA_real_
    z$stabilization_status <- ifelse(keep,'pending',rows$ordinary_na_reason)
    if (length(valid_genes)) {
      index <- match(valid_genes,rownames(dds))
      set.seed(812745)
      fit <- tryCatch(lfcShrink(dds[index,],res=res[index,],type='ashr',parallel=FALSE,quiet=TRUE),error=identity)
      if (inherits(fit,'error')) {
        z$stabilization_status[keep] <- paste0('failed: ',gsub('[\t\r\n]+',' ',conditionMessage(fit)))
      } else {
        finite <- is.finite(fit$log2FoldChange) & is.finite(fit$lfcSE) & fit$lfcSE>=0
        z$stabilized_log2FoldChange[keep] <- ifelse(finite,fit$log2FoldChange,NA_real_)
        z$stabilized_lfcSE[keep] <- ifelse(finite,fit$lfcSE,NA_real_)
        z$stabilization_status[keep] <- ifelse(finite,'ashr_posterior_mean_and_sd','invalid_shrinkage_output')
        saveRDS(fit,file.path(output,paste0(id,'.rds')))
      }
    }
    tables[[id]] <- z
  }
  write.table(do.call(rbind,tables),file.path(output,'stabilized_display.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
  capture.output(sessionInfo(),file=file.path(output,'sessionInfo.txt'))
  invisible(do.call(rbind,tables))
}

if (sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  if (length(args)!=1L) stop('Supply exactly one model_id')
  heat_shrink(args[1])
}
