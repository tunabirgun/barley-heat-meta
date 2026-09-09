# R 4.5.2; DESeq2 1.50.2. One serial worker.
suppressPackageStartupMessages(library(DESeq2))

heat_input <- function(model_id) {
  samples <- read.delim('config/heat_samples.tsv', check.names=FALSE)
  models <- read.delim('config/heat_models.tsv', check.names=FALSE)
  genes <- read.delim('config/gene_descriptions.tsv', check.names=FALSE)$gene_id
  stopifnot(all(c('model_id','sample','run','genotype','condition','replicate','tissue','time') %in% names(samples)),
    !anyDuplicated(samples$sample), !anyDuplicated(samples$run), !anyNA(samples$sample), !anyNA(samples$run),
    !anyDuplicated(models$model_id), setequal(samples$model_id, models$model_id),
    !anyDuplicated(genes), !anyNA(genes), all(nzchar(genes)),
    all(models$analysis_tier=='conditional_exploratory'), all(models$model_formula=='~ 0 + cell'))
  model <- models[models$model_id==model_id,,drop=FALSE]
  stopifnot(nrow(model)==1L)
  samples <- samples[samples$model_id==model_id,,drop=FALSE]
  input <- read.delim(file.path('results/heat/counts',paste0(model_id,'.txt')), comment.char='#', check.names=FALSE)
  stopifnot(identical(names(input)[1:6],c('Geneid','Chr','Start','End','Strand','Length')),
    !anyDuplicated(input$Geneid), setequal(input$Geneid,genes), ncol(input)==nrow(samples)+6L)
  counts <- as.matrix(input[,7:ncol(input),drop=FALSE])
  rownames(counts) <- input$Geneid
  colnames(counts) <- sub('_Aligned.sortedByCoord.out.bam$','',basename(colnames(counts)))
  if (!is.numeric(counts) || any(!is.finite(counts) | counts<0 | counts!=round(counts)))
    stop('Counts must be finite non-negative integers')
  stopifnot(!anyDuplicated(colnames(counts)),setequal(colnames(counts),samples$run),
    all(samples$condition %in% c('control','heat')), !anyNA(samples$genotype),
    !anyNA(samples$replicate), all(nzchar(samples$genotype)))
  samples <- samples[match(colnames(counts),samples$run),,drop=FALSE]
  rownames(samples) <- samples$run
  samples$cell <- factor(paste(samples$genotype,samples$condition,sep='_'))
  stopifnot(!anyDuplicated(samples[c('cell','replicate')]),all(table(samples$cell)>=2L))
  counts <- counts[match(genes,rownames(counts)),,drop=FALSE]
  nmin <- min(table(samples$cell))
  keep <- unname(rowSums(counts>=10L)>=nmin)
  if (sum(keep)<2L) stop('Too few genes after the count prefilter')
  design <- model.matrix(~0+cell,samples)
  stopifnot(qr(design)$rank==ncol(design),nrow(design)>ncol(design))
  list(counts=counts,samples=samples,model=model,keep=keep,smallest_cell_n=unname(nmin))
}

heat_write <- function(x,path) write.table(x,path,sep='\t',quote=FALSE,row.names=FALSE,na='NA')

heat_sample_qc <- function(model_id) {
  input <- heat_input(model_id)
  out <- file.path('results/heat/sample_qc',model_id)
  if (dir.exists(out)) stop('QC output already exists; preserve and review before rerunning')
  dds <- DESeqDataSetFromMatrix(input$counts[input$keep,,drop=FALSE],input$samples,~0+cell)
  transformed <- assay(varianceStabilizingTransformation(dds,blind=TRUE))
  gene_var <- apply(transformed,1,var)
  selected <- head(order(-gene_var,rownames(transformed)),500L)
  pca <- prcomp(t(transformed[selected,,drop=FALSE]),center=TRUE,scale.=FALSE)
  if (!is.finite(sum(pca$sdev^2)) || sum(pca$sdev^2)<=0) stop('PCA has no finite variation')
  dir.create(out,recursive=TRUE)
  heat_write(data.frame(input$samples,pca$x,check.names=FALSE),file.path(out,'pca.tsv'))
  heat_write(data.frame(component=paste0('PC',seq_along(pca$sdev)),variance=pca$sdev^2/sum(pca$sdev^2)),
    file.path(out,'pca_variance.tsv'))
  heat_write(data.frame(run=colnames(transformed),cor(transformed),check.names=FALSE),
    file.path(out,'sample_correlations.tsv'))
  heat_write(data.frame(gene_id=rownames(input$counts),prefilter_pass=input$keep,
    samples_at_least_10=rowSums(input$counts>=10L),required_samples=input$smallest_cell_n),
    file.path(out,'prefilter.tsv'))
  heat_write(data.frame(gene_id=rownames(transformed)[selected],variance=gene_var[selected]),
    file.path(out,'pca_genes.tsv'))
  heat_write(data.frame(input$samples,library_assigned_fragments=colSums(input$counts),
    detected_genes=colSums(input$counts>0)),file.path(out,'sample_metrics.tsv'))
  heat_write(data.frame(model_id,analysis_tier='conditional_exploratory',qc_review='pending',
    review_note='Review technical diagnostics and biological sample identity; PCA separation alone is not an exclusion criterion.'),
    file.path(out,'review_input.tsv'))
  saveRDS(transformed,file.path(out,'vst.rds'))
  capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
}

if (sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  if (length(args)!=1L) stop('Supply exactly one model_id')
  heat_sample_qc(args[1])
}
