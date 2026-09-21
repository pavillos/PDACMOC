#' @title PDAConsensus label
#'
#' @author Villoslada-Blanco, Pablo
#'
#' @description Creates the consistently styled PDAConsensus classifier label.
#'
#' @param suffix Optional plain-text suffix appended to the classifier name.
#'
#' @return A \code{shiny::HTML} object containing the italicized label and suffix.

pdaconsensus.label <- function(suffix = NULL) {
  suffix <- if (is.null(suffix) || !nzchar(suffix)) "" else paste0(" ", suffix)
  shiny::HTML(paste0("<em>PDAConsensus</em>", suffix))
}
