#' @title Usage counter of the Shiny app
#'
#' @author Villoslada-Blanco, Pablo
#'
#' @description
#' Counts the classifications finished in the Shiny app and the samples classified,
#' for the usage badge of the README. Nothing about the users or their data is stored:
#' only the two totals, the date the count started and the date of the last update.
#' The counter is off unless \code{options(PDACMOC.stats_dir = <folder>)} is set, which
#' is done only on the public server. The folder then holds \code{usage.json}, which
#' the app serves read-only at \code{/stats/usage.json}.
#'
#' @name usage
NULL

usage.file <- function() {
  dir <- getOption('PDACMOC.stats_dir')
  if (is.null(dir)) NULL else file.path(dir, 'usage.json')
}

#' @describeIn usage Current totals (zero when the file does not exist yet)
#' @param file Path of usage.json
read.usage <- function(file) {
  usage <- list(classifications = 0, samples = 0, since = format(Sys.Date()), updated = '')
  if (file.exists(file)) {
    text <- paste(readLines(file, warn = FALSE), collapse = '')
    for (field in names(usage)) {
      pattern <- sprintf('"%s"\\s*:\\s*"?([^",}]*)"?', field)
      value <- regmatches(text, regexec(pattern, text))[[1]]
      if (length(value) == 2) {
        usage[[field]] <- if (field %in% c('classifications', 'samples')) as.numeric(value[2]) else value[2]
      }
    }
    # an unreadable file is left as it is rather than restarting the count from zero
    if (anyNA(c(usage$classifications, usage$samples)) || !grepl('"classifications"', text)) {
      stop('usage file could not be read')
    }
  }
  usage
}

#' @describeIn usage Add one finished classification of \code{n_samples} samples
#' @param n_samples Number of samples classified
record.usage <- function(n_samples, file = usage.file()) {
  if (is.null(file)) return(invisible(NULL))
  # the counter must never stop the app
  tryCatch({
    usage <- read.usage(file)
    usage$classifications <- usage$classifications + 1
    usage$samples <- usage$samples + n_samples
    write.usage(usage, file)
  }, error = function(e) NULL)
  invisible(NULL)
}

write.usage <- function(usage, file) {
  usage$updated <- format(Sys.time(), '%Y-%m-%d %H:%M')
  json <- sprintf('{"classifications": %.0f, "samples": %.0f, "since": "%s", "updated": "%s"}',
                  usage$classifications, usage$samples, usage$since, usage$updated)
  # write under another name first so the file is never read half-written
  writeLines(json, paste0(file, '.tmp'))
  file.rename(paste0(file, '.tmp'), file)
}

#' @describeIn usage Sentence for the footer of the app, or NULL when the counter is off
usage.sentence <- function(file = usage.file()) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  usage <- tryCatch(read.usage(file), error = function(e) NULL)
  if (is.null(usage)) return(NULL)
  since <- tryCatch(format(as.Date(usage$since), '%B %Y'), error = function(e) usage$since)
  plural <- function(n, word) sprintf('%s %s%s', format(n, big.mark = ','), word, if (n == 1) '' else 's')
  sprintf('Used for %s (%s) since %s', plural(usage$classifications, 'classification'),
          plural(usage$samples, 'sample'), since)
}

# Serve the stats folder at /stats (only usage.json is written there)
serve.usage <- function() {
  dir <- getOption('PDACMOC.stats_dir')
  if (!is.null(dir)) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    file <- usage.file()
    if (!file.exists(file)) write.usage(read.usage(file), file)
    shiny::addResourcePath('stats', normalizePath(dir))
  }
}
