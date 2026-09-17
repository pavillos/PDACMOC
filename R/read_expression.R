#' @title Collapse duplicated genes
#'
#' @author Villoslada-Blanco, Pablo
#'
#' @description
#' Average the rows of an expression matrix that share the same gene ID (the same
#' rule applied when Ensembl ID versions are discarded). If there are no duplicated
#' IDs, the IDs are just set as row names.
#'
#' @param expmat Matrix or data frame of RNAseq raw counts (genes in rows)
#' @param gene_ids Character vector with the gene ID of each row
#'
#' @return expmat with one row per gene ID

collapse.duplicated.genes <- function(expmat, gene_ids) {

  gene_ids <- as.character(gene_ids)
  if (!anyDuplicated(gene_ids)) {
    rownames(expmat) <- gene_ids
    return(expmat)
  }

  n_duplicated <- sum(duplicated(gene_ids))
  sums <- rowsum(as.matrix(expmat), gene_ids, reorder = FALSE)
  counts <- rowsum(rep(1, length(gene_ids)), gene_ids, reorder = FALSE)
  collapsed <- as.data.frame(sums / as.vector(counts))
  message(n_duplicated, ' rows with duplicated gene IDs were averaged')

  return(collapsed)

}

#' @title Read expression file
#'
#' @author Villoslada-Blanco, Pablo
#'
#' @description
#' Read a tab-separated file with gene IDs in the first column and samples in the
#' remaining columns. Rows with duplicated gene IDs are averaged. The header may
#' include a name for the gene ID column or omit it.
#'
#' @param path Path to the tsv file
#'
#' @return data frame with genes in rows and samples in columns

read.expression.file <- function(path) {

  first_lines <- readLines(path, n = 2, warn = FALSE)
  if (length(first_lines) < 2) {
    stop('The file must have a header and at least one gene')
  }
  header <- strsplit(first_lines[1], '\t', fixed = TRUE)[[1]]
  n_fields <- length(strsplit(first_lines[2], '\t', fixed = TRUE)[[1]])

  if (n_fields == length(header) + 1) {
    # Header without a name for the gene ID column
    expmat <- read.delim(path, header = FALSE, skip = 1, check.names = FALSE)
    colnames(expmat) <- c('GeneID', gsub('^"|"$', '', header))
  } else {
    expmat <- read.delim(path, header = TRUE, check.names = FALSE)
  }
  if (ncol(expmat) < 2) {
    stop('The file must have gene IDs in the first column and at least one sample column')
  }

  gene_ids <- expmat[[1]]
  expmat <- collapse.duplicated.genes(expmat[, -1, drop = FALSE], gene_ids)

  return(expmat)

}
