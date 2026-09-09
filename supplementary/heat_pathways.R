# R 4.5.2; fgsea 1.36.2. One serial worker.
suppressPackageStartupMessages(library(fgsea))

heat_pathways <- function() {
  ordinary <- read.delim('results/heat/inference/ordinary_results.tsv',check.names=FALSE)
  contrasts <- read.delim('config/heat_contrasts.tsv',check.names=FALSE,quote='')
  annotation <- read.delim('config/gene_go.tsv',check.names=FALSE)
  stopifnot(all(ordinary$analysis_tier=='conditional_exploratory'),
    all(ordinary$multiplicity_status=='final_declared_conditional_family'),
    setequal(ordinary$contrast_id,contrasts$contrast_id),
    identical(ordinary$comparison_family,contrasts$comparison_family[match(ordinary$contrast_id,contrasts$contrast_id)]),
    !anyDuplicated(ordinary[c('contrast_id','gene_id')]),
    !anyDuplicated(annotation[c('gene_id','pathway_id')]),
    !anyNA(annotation$gene_id),!anyNA(annotation$pathway_id))
  universe <- unique(annotation[c('pathway_id','pathway_name')])
  stopifnot(!anyDuplicated(universe$pathway_id))
  sets <- split(annotation$gene_id,annotation$pathway_id)
  out <- 'results/heat/pathways'
  if (dir.exists(out)) stop('Pathway output already exists; preserve and review before rerunning')
  dir.create(out,recursive=TRUE)
  tables <- audits <- warnings <- list()
  for (id in unique(ordinary$contrast_id)) {
    rows <- ordinary[ordinary$contrast_id==id,,drop=FALSE]
    valid <- rows$ordinary_valid %in% TRUE & rows$beta_converged %in% TRUE &
      is.finite(rows$wald_stat) & is.finite(rows$pvalue)
    ranks <- setNames(rows$wald_stat[valid],rows$gene_id[valid])
    ranks <- ranks[order(-ranks,names(ranks))]
    mapped <- lengths(lapply(sets,intersect,names(ranks)))
    eligible <- names(sets)[mapped>=15L & mapped<=500L]
    result <- data.frame(pathway_id=character(),pval=numeric(),padj=numeric(),log2err=numeric(),
      ES=numeric(),NES=numeric(),size=integer(),leadingEdge=character())
    messages <- character()
    if (length(eligible)) {
      set.seed(812745)
      result <- withCallingHandlers(as.data.frame(fgsea::fgsea(pathways=sets[eligible],stats=ranks,
        minSize=15L,maxSize=500L,eps=0,nproc=1L,BPPARAM=BiocParallel::SerialParam(progressbar=FALSE))),
        warning=function(w) { messages <<- c(messages,conditionMessage(w)); invokeRestart('muffleWarning') })
      saveRDS(result,file.path(out,paste0(id,'.rds')))
      result$leadingEdge <- vapply(result$leadingEdge,paste,'',collapse=';')
      names(result)[names(result)=='pathway'] <- 'pathway_id'
    }
    z <- merge(universe,result,by='pathway_id',all.x=TRUE,sort=FALSE)
    names(z)[names(z)=='pval'] <- 'pvalue'
    names(z)[names(z)=='padj'] <- 'native_padj'
    z$mapped_set_size <- unname(mapped[z$pathway_id])
    z$test_eligible <- z$pathway_id %in% eligible
    z$status <- ifelse(!z$test_eligible,'not_eligible_mapped_size',
      ifelse(is.finite(z$pvalue) & z$pvalue>=0 & z$pvalue<=1 & is.finite(z$NES),'tested','test_failed'))
    z$pvalue[z$status!='tested'] <- z$native_padj[z$status!='tested'] <- NA_real_
    for (field in c('model_id','publication_family_id','contrast_id','comparison_family',
      'contrast_role','analysis_tier','uncertainty')) {
      stopifnot(length(unique(rows[[field]]))==1L)
      z[[field]] <- rows[[field]][1]
    }
    tables[[id]] <- z
    audits[[id]] <- data.frame(contrast_id=id,analysis_tier='conditional_exploratory',
      valid_ranked_genes=length(ranks),annotated_ranked_genes=sum(names(ranks) %in% annotation$gene_id),
      tied_rank_values=sum(duplicated(unname(ranks))),positive_ranks=sum(ranks>0),negative_ranks=sum(ranks<0),
      eligible_gene_sets=length(eligible),tested_gene_sets=sum(z$status=='tested'),seed=812745L,workers=1L)
    if (length(messages)) warnings[[id]] <- data.frame(contrast_id=id,warning=gsub('[\t\r\n]+',' ',messages))
  }
  tables <- do.call(rbind,tables)
  tables$project_family_q <- NA_real_
  tables$family_n_tests <- NA_integer_
  for (family in unique(tables$comparison_family)) {
    selected <- tables$comparison_family==family
    valid <- selected & tables$status=='tested'
    denominator <- sum(selected & tables$test_eligible)
    tables$family_n_tests[selected] <- denominator
    tables$project_family_q[valid] <- p.adjust(tables$pvalue[valid],'BH',n=denominator)
  }
  write.table(tables,file.path(out,'go_enrichment.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
  write.table(do.call(rbind,audits),file.path(out,'rank_audit.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
  warning_table <- if (length(warnings)) do.call(rbind,warnings) else data.frame(contrast_id=character(),warning=character())
  write.table(warning_table,file.path(out,'warnings.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
  capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
  invisible(tables)
}

if (sys.nframe()==0L) heat_pathways()
