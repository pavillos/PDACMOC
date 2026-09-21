#' @title Run panel
#'
#' @author Villoslada-Blanco, Pablo
#'
#' @description
#' Content of the Run card of the Shiny app for each state of a classification:
#' idle, waiting for a free slot, running (steps, progress, elapsed and estimated
#' time, Cancel button), finished, cancelled or failed.
#'
#' @param status List with the state and, while running, the job and its progress
#' @param df Uploaded expression matrix, or NULL
#' @param estimate Estimated minutes for the uploaded file and the selected options
#' @param messages Data frame of messages (type and text) about the file and the classification
#'
#' @return An HTML tag

run.panel <- function(status, df, estimate, messages = NULL) {
  tagList(run.state(status, df, estimate), message.list(messages))
}

# Messages from reading the file and from the classification, below the state
message.list <- function(messages) {
  if (is.null(messages) || nrow(messages) == 0) return(NULL)
  div(class = 'run-messages',
      div(class = 'run-messages-title', 'Messages'),
      tags$ul(lapply(seq_len(nrow(messages)), function(i) {
        tags$li(class = messages$type[i],
                span(class = 'icon', if (messages$type[i] == 'error') '!' else 'i'),
                messages$text[i])
      })))
}

run.state <- function(status, df, estimate) {

  file_line <- if (!is.null(df)) {
    div(class = 'helper', style = 'margin: .75rem 0 0',
        sprintf('Your file: %d samples. Estimated time with the selected options: about %s.',
                ncol(df), minutes.label(estimate)))
  }
  idle_steps <- div(
    class = 'empty-state',
    div(class = 'step', span(class = 'num', '1'), 'Upload a .tsv or .csv file of raw counts.'),
    div(class = 'step', span(class = 'num', '2'), 'Choose the gene ID type and the classifiers.'),
    div(class = 'step', span(class = 'num', '3'), 'Press Run classification. You can keep using the app in another tab while it runs.'))

  switch(
    status$state,
    idle = tagList(idle_steps, file_line),

    done = tagList(
      div(class = 'run-banner done',
          sprintf('Classification of %d samples finished in %s.', status$n_samples,
                  minutes.label(status$minutes)),
          actionLink('goToResults', 'See the results')),
      step.list(status$progress, 'done')),

    cancelled = tagList(div(class = 'run-banner', 'Classification cancelled.'),
                        step.list(status$progress, 'stopped'), idle_steps, file_line),

    error = tagList(div(class = 'run-banner error', paste('The classification could not be completed:', status$message)),
                    step.list(status$progress, 'stopped'), idle_steps, file_line),

    waiting = tagList(
      div(class = 'run-head', span(class = 'run-title', 'Waiting for a free slot'),
          actionButton('cancelButton', 'Cancel', class = 'btn-outline-primary btn-sm')),
      div(class = 'helper', style = 'margin-top:.5rem',
          'Other classifications are running on the server. Yours starts automatically as soon as one finishes.')),

    running = {
      job <- status$job
      progress <- status$progress
      elapsed <- as.numeric(difftime(Sys.time(), job$started %||% job$submitted, units = 'mins'))
      value <- if (is.null(progress)) 0.02 else max(0.02, utils::tail(progress$value, 1))
      steps <- step.list(progress, 'now')
      tagList(
        div(class = 'run-head',
            span(class = 'run-title', sprintf('Classifying %d samples', job$n_samples)),
            span(class = 'run-time', sprintf('%s elapsed · about %s in total',
                                             elapsed.label(elapsed), minutes.label(job$estimate)))),
        div(class = 'run-bar', div(style = sprintf('width:%.0f%%', 100 * value))),
        steps %||% div(class = 'helper', 'Starting the classification...'),
        div(class = 'run-actions',
            span(class = 'helper', 'You can keep using the app in another tab.'),
            actionButton('cancelButton', 'Cancel', class = 'btn-outline-primary btn-sm')))
    }
  )
}

# Steps reported by the classification; the last one is marked as current, done or
# the step where it stopped
step.list <- function(progress, last) {
  if (is.null(progress)) return(NULL)
  steps <- progress$step[progress$step != 'Classification started']
  steps <- sub('\\. This may take a while\\.\\.\\.', '', steps)
  if (length(steps) == 0) return(NULL)
  tags$ul(class = 'run-steps', lapply(seq_along(steps), function(i) {
    tags$li(class = if (i == length(steps)) last else 'done', tags$i(class = 'dot'), steps[i])
  }))
}

minutes.label <- function(minutes) {
  if (minutes < 1.5) '1 minute' else sprintf('%.0f minutes', minutes)
}

# Elapsed time as m:ss
elapsed.label <- function(minutes) {
  seconds <- round(minutes * 60)
  sprintf('%d:%02d', seconds %/% 60, seconds %% 60)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
