#' @title Background classification
#'
#' @author Villoslada-Blanco, Pablo
#'
#' @description
#' Run \code{omni.classify} in a separate R process so that the Shiny app stays
#' responsive and a classification can be cancelled. Each classification ("job") lives
#' in its own temporary folder: the input, a progress file written by
#' \code{report.progress}, the package messages, and the result or the error. At most
#' \code{getOption('PDACMOC.max_jobs', 2)} jobs run at the same time in one app; later
#' ones wait for a free slot.
#'
#' @name background
NULL

# Jobs running in this R process (folder -> process id)
.background_jobs <- new.env(parent = emptyenv())

max.background.jobs <- function() getOption('PDACMOC.max_jobs', 2)

process.alive <- function(pid) {
  !is.na(pid) && file.exists(file.path('/proc', pid))
}

read.pid <- function(dir) {
  f <- file.path(dir, 'pid')
  if (file.exists(f)) suppressWarnings(as.integer(readLines(f, n = 1, warn = FALSE))) else NA_integer_
}

running.background.jobs <- function() {
  dirs <- ls(.background_jobs)
  alive <- vapply(dirs, function(d) {
    pid <- read.pid(d)
    # a job that has not written its pid yet is still starting
    is.na(pid) || process.alive(pid)
  }, logical(1))
  sum(alive)
}

#' @describeIn background Rough running time in minutes, from timings measured on the
#' server (fixed start-up plus a cost per sample, larger with stroma)
#' @param n_samples Number of samples
#' @param batch,stroma Options of the classification
estimate.minutes <- function(n_samples, batch, stroma) {
  0.4 + (if (batch) 4 else 0.5) + n_samples * (if (stroma) 0.12 else 0.03)
}

#' @describeIn background Prepare a job; it starts as soon as a slot is free
#' @param expmat Expression matrix as read by \code{read.counts}
#' @param args Named list of arguments for \code{omni.classify} (without expmat)
background.submit <- function(expmat, args) {
  dir <- tempfile('pdacmoc_job_')
  dir.create(dir)
  saveRDS(list(expmat = expmat, args = args), file.path(dir, 'input.rds'))
  list(dir = dir, submitted = Sys.time(), started = NULL, n_samples = ncol(expmat),
       estimate = estimate.minutes(ncol(expmat), isTRUE(args$batch), isTRUE(args$stroma)))
}

# Start the worker process of a job
background.start <- function(job) {
  python <- tryCatch(reticulate::py_config()$python, error = function(e) '')
  rscript <- file.path(R.home('bin'), 'Rscript')
  expr <- sprintf(
    "%s PDACMOC:::background.worker('%s')",
    if (nzchar(python)) sprintf("reticulate::use_python('%s', required = TRUE);", python) else '',
    job$dir)
  env <- c(sprintf('R_LIBS=%s', paste(.libPaths(), collapse = ':')),
           sprintf('LC_ALL=%s', Sys.getenv('LC_ALL', 'C.UTF-8')))
  # setsid puts the worker and its parallel children in one process group,
  # so that cancelling stops all of them
  system2('setsid', c(shQuote(rscript), '-e', shQuote(expr)), wait = FALSE,
          stdout = file.path(job$dir, 'worker.log'), stderr = file.path(job$dir, 'worker.log'),
          env = env)
  assign(job$dir, TRUE, envir = .background_jobs)
  job$started <- Sys.time()
  job
}

#' @describeIn background State of a job; starts it when a slot is free
#' @param job Job returned by \code{background.submit}
background.status <- function(job) {
  dir <- job$dir
  if (is.null(job$started)) {
    if (running.background.jobs() < max.background.jobs()) {
      job <- background.start(job)
      return(list(job = job, state = 'running', progress = NULL, messages = character(0)))
    }
    return(list(job = job, state = 'waiting', progress = NULL, messages = character(0)))
  }
  progress <- NULL
  f <- file.path(dir, 'progress.tsv')
  if (file.exists(f)) {
    lines <- readLines(f, warn = FALSE)
    lines <- lines[nzchar(lines)]
    if (length(lines) > 0) {
      parts <- strsplit(lines, '\t', fixed = TRUE)
      progress <- data.frame(value = as.numeric(vapply(parts, `[`, '', 1)),
                             time = vapply(parts, `[`, '', 2),
                             step = vapply(parts, `[`, '', 3), stringsAsFactors = FALSE)
    }
  }
  m <- file.path(dir, 'messages.txt')
  messages <- if (file.exists(m)) readLines(m, warn = FALSE) else character(0)
  state <- if (file.exists(file.path(dir, 'result.rds'))) {
    'done'
  } else if (file.exists(file.path(dir, 'error.txt'))) {
    'error'
  } else if (!is.na(read.pid(dir)) && !process.alive(read.pid(dir))) {
    'lost'
  } else if (is.na(read.pid(dir)) &&
             as.numeric(difftime(Sys.time(), job$started, units = 'secs')) > 120) {
    'lost'
  } else {
    'running'
  }
  list(job = job, state = state, progress = progress, messages = messages)
}

#' @describeIn background Result of a finished job
background.result <- function(job) readRDS(file.path(job$dir, 'result.rds'))

#' @describeIn background Error message of a failed job
background.error <- function(job) {
  f <- file.path(job$dir, 'error.txt')
  if (file.exists(f)) return(paste(readLines(f, warn = FALSE), collapse = ' '))
  log <- file.path(job$dir, 'worker.log')
  detail <- if (file.exists(log)) utils::tail(readLines(log, warn = FALSE), 3) else character(0)
  paste(c('the background process stopped unexpectedly', detail), collapse = ' ')
}

#' @describeIn background Stop a job and all its processes
background.cancel <- function(job) {
  pid <- read.pid(job$dir)
  if (process.alive(pid)) {
    # negative pid: the whole process group created by setsid
    system2('kill', c('-TERM', paste0('-', pid)), stdout = FALSE, stderr = FALSE)
  }
  background.cleanup(job)
}

#' @describeIn background Forget a job and delete its folder
background.cleanup <- function(job) {
  if (exists(job$dir, envir = .background_jobs, inherits = FALSE)) {
    rm(list = job$dir, envir = .background_jobs)
  }
  unlink(job$dir, recursive = TRUE)
}

#' @describeIn background Entry point of the worker process
#' @param job_dir Folder of the job
background.worker <- function(job_dir) {
  writeLines(as.character(Sys.getpid()), file.path(job_dir, 'pid'))
  options(PDACMOC.progress_file = file.path(job_dir, 'progress.tsv'))
  input <- readRDS(file.path(job_dir, 'input.rds'))
  message_file <- file.path(job_dir, 'messages.txt')
  result <- tryCatch(
    withCallingHandlers(
      do.call(omni.classify, c(list(expmat = input$expmat), input$args)),
      message = function(m) {
        cat(trimws(conditionMessage(m)), '\n', sep = '', file = message_file, append = TRUE)
        invokeRestart('muffleMessage')
      }),
    error = function(e) e)
  if (inherits(result, 'error')) {
    writeLines(conditionMessage(result), file.path(job_dir, 'error.txt'))
  } else {
    # write under another name first so the app never reads a half-written file
    saveRDS(result, file.path(job_dir, 'result.part'))
    file.rename(file.path(job_dir, 'result.part'), file.path(job_dir, 'result.rds'))
  }
  invisible(NULL)
}
