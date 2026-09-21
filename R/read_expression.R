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
#' @usage collapse.duplicated.genes(expmat, gene_ids)

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

#' @title Read counts
#'
#' @author Villoslada-Blanco, Pablo
#'
#' @description
#' Read a file of raw counts with gene IDs in the first column and samples in the
#' remaining columns, as the Shiny app does. The separator (tab, comma or semicolon)
#' is detected from the header, so .tsv, .csv and .txt files are accepted. Rows with
#' duplicated gene IDs are averaged. The header may include a name for the gene ID
#' column or omit it. The result can be passed to \code{omni.classify}.
#'
#' @param path Path to the file
#'
#' @return data frame with genes in rows and samples in columns
#'
#' @examples
#' \dontrun{
#' counts <- read.counts('my_counts.csv')
#' classification <- omni.classify(counts, gene_id = 'EnsemblID')
#' }
#'
#' @aliases read.expression.file
#' @export

read.counts <- function(path) {

  first_lines <- readLines(path, n = 2, warn = FALSE)
  if (length(first_lines) < 2) {
    stop('The file must have a header and at least one gene')
  }

  # Separator: tab first (as in previous versions), then comma, then semicolon
  sep <- if (grepl('\t', first_lines[1], fixed = TRUE)) {
    '\t'
  } else if (grepl(',', first_lines[1], fixed = TRUE)) {
    ','
  } else if (grepl(';', first_lines[1], fixed = TRUE)) {
    ';'
  } else {
    stop('The separator could not be detected: use tabs, commas or semicolons between columns')
  }
  read_file <- function(...) {
    if (sep == '\t') read.delim(path, check.names = FALSE, ...)
    else read.table(path, sep = sep, quote = '"', comment.char = '', fill = TRUE,
                    check.names = FALSE, ...)
  }

  header <- strsplit(first_lines[1], sep, fixed = TRUE)[[1]]
  n_fields <- length(strsplit(first_lines[2], sep, fixed = TRUE)[[1]])

  if (n_fields == length(header) + 1) {
    # Header without a name for the gene ID column
    expmat <- read_file(header = FALSE, skip = 1)
    colnames(expmat) <- c('GeneID', gsub('^"|"$', '', header))
  } else {
    expmat <- read_file(header = TRUE)
  }
  if (ncol(expmat) < 2) {
    stop('The file must have gene IDs in the first column and at least one sample column')
  }

  gene_ids <- expmat[[1]]
  expmat <- collapse.duplicated.genes(expmat[, -1, drop = FALSE], gene_ids)

  return(expmat)

}

# Former internal name, kept so existing scripts keep working
read.expression.file <- read.counts
