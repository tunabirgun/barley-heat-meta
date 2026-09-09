# R 4.5.2; metafor 5.0-1. Run after heat_synthesis_01_families.R.
args <- commandArgs(TRUE)
out <- if (length(args)) args[1] else 'results/heat/synthesis_revised'
stopifnot(!file.exists(file.path(out, 'partial_conjunction.tsv')))
fam <- read.delim(file.path(out, 'family_effects.tsv'), quote = '', check.names = FALSE)
stopifnot(!anyDuplicated(fam[c('gene_id', 'family')]), all(fam$p_up >= 0 & fam$p_up <= 1), all(fam$p_down >= 0 & fam$p_down <= 1))
genes <- sort(unique(fam$gene_id)); families <- sort(unique(fam$family)); n <- length(families)
stopifnot(n >= 2, nrow(fam) == length(genes) * n)
key <- paste(fam$gene_id, fam$family, sep = '\t')
grid <- expand.grid(family = families, gene_id = genes, stringsAsFactors = FALSE)
ordered <- fam[match(paste(grid$gene_id, grid$family, sep = '\t'), key), ]
up <- matrix(ordered$p_up, ncol = n, byrow = TRUE)
down <- matrix(ordered$p_down, ncol = n, byrow = TRUE)
pc <- data.frame(gene_id = genes, families_available = rowSums(matrix(ordered$complete, ncol = n, byrow = TRUE)))
previous_up <- previous_down <- rep(0, length(genes))
p_columns <- character()
for (u in 2:n) {
  p_up <- pmin(1, (n-u+1) * apply(up, 1, sort)[u, ])
  p_down <- pmin(1, (n-u+1) * apply(down, 1, sort)[u, ])
  previous_up <- pmax(previous_up, p_up); previous_down <- pmax(previous_down, p_down)
  pc[[paste0('p_up_at_least_', u)]] <- previous_up
  pc[[paste0('p_down_at_least_', u)]] <- previous_down
  column <- paste0('p_at_least_', u)
  pc[[column]] <- pmin(1, 2 * pmin(previous_up, previous_down))
  pc[[paste0('direction_at_least_',u)]] <- ifelse(previous_up < previous_down, 'up', ifelse(previous_down < previous_up, 'down', 'undetermined'))
  p_columns <- c(p_columns, column)
}
q <- matrix(p.adjust(as.matrix(pc[p_columns]), 'BH'), nrow = nrow(pc))
for (i in seq_along(p_columns)) pc[[sub('^p_', 'q_', p_columns[i])]] <- q[, i]
pc$direction <- pc[[paste0('direction_at_least_', n)]]
pc$adjustment <- 'two_directions_Bonferroni_then_BH_all_genes_and_u'
write.table(pc, file.path(out, 'partial_conjunction.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)

split_fam <- split(fam[fam$complete, ], fam$gene_id[fam$complete])
results <- vector('list', length(genes)); issues <- list()
for (i in seq_along(genes)) {
  g <- genes[i]; f <- split_fam[[g]]
  r <- data.frame(gene_id = g, k_families = if (is.null(f)) 0 else nrow(f), reml_estimate = NA_real_, reml_se = NA_real_,
    reml_ci_low = NA_real_, reml_ci_high = NA_real_, reml_pvalue = NA_real_, reml_tau2 = NA_real_, reml_I2 = NA_real_,
    df = NA_real_, dl_estimate = NA_real_, status = 'fewer_than_two_complete_families')
  if (!is.null(f) && nrow(f) >= 2) {
    warnings <- character()
    fit <- tryCatch(withCallingHandlers(metafor::rma.uni(yi = f$estimate, sei = f$se, method = 'REML', test = 'adhoc', control = list(stepadj = .5, maxiter = 1000)),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning') }), error = function(e) e)
    if (inherits(fit, 'error')) {
      r$status <- 'fit_failed'; issues[[length(issues)+1]] <- data.frame(gene_id=g, issue=conditionMessage(fit))
    } else {
      r$reml_estimate <- as.numeric(fit$b); r$reml_se <- fit$se
      r$reml_ci_low <- fit$ci.lb; r$reml_ci_high <- fit$ci.ub; r$reml_pvalue <- fit$pval
      r$reml_tau2 <- fit$tau2; r$reml_I2 <- fit$I2; r$df <- fit$ddf; r$status <- if (length(warnings)) 'estimated_with_warning' else 'estimated'
      w <- 1/f$se^2; fixed <- sum(w*f$estimate)/sum(w)
      tau <- max(0,(sum(w*(f$estimate-fixed)^2)-(nrow(f)-1))/(sum(w)-sum(w^2)/sum(w)))
      wr <- 1/(f$se^2+tau); r$dl_estimate <- sum(wr*f$estimate)/sum(wr)
      if (length(warnings)) issues[[length(issues)+1]] <- data.frame(gene_id=g, issue=paste(unique(warnings),collapse='; '))
    }
  }
  results[[i]] <- r
}
re <- do.call(rbind, results)
re$reml_q <- p.adjust(ifelse(is.na(re$reml_pvalue), 1, re$reml_pvalue), 'BH')
write.table(re, file.path(out, 'random_effects.tsv'), sep='\t', quote=FALSE, row.names=FALSE, na='NA')
issue_table <- if (length(issues)) do.call(rbind, issues) else data.frame(gene_id=character(), issue=character())
issue_table$issue <- gsub('[\r\n\t]+', ' ', issue_table$issue)
write.table(issue_table, file.path(out, 'fit_issues.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
audit <- data.frame(genes=length(genes), families=n, genes_with_all_families=sum(pc$families_available==n),
  genes_all_families_q05=sum(pc[[paste0('q_at_least_',n)]]<=.05), reml_estimated=sum(!is.na(re$reml_estimate)),
  reml_failed=sum(re$status=='fit_failed'), warning_genes=nrow(issue_table), inference='modified_Knapp_Hartung_covariance_bound',
  moderator_status='withdrawn_heterogeneous_contexts_and_six_confounded_families')
write.table(audit, file.path(out, 'audit.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
capture.output(sessionInfo(), file=file.path(out, 'sessionInfo.txt'))
