#' @title Report progress
#'
#' @author Villoslada-Blanco, Pablo
#'
#' @description
#' Write a classification step to the progress file given by the option
#' \code{PDACMOC.progress_file}, if set. Used when the classification runs in a
#' background process so that the Shiny app can show which step is running. It only
#' reports progress: it does not change any result.
#'
#' @param message Description of the step
#' @param value Fraction of the work done when the step starts (0 to 1)
#'
#' @return Nothing

report.progress <- function(message, value) {
  progress_file <- getOption('PDACMOC.progress_file')
  if (!is.null(progress_file)) {
    cat(sprintf('%.2f\t%s\t%s\n', value, format(Sys.time(), '%H:%M:%S'), message),
        file = progress_file, append = TRUE)
  }
  invisible(NULL)
}
