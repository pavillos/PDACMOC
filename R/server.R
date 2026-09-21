#' @title Server
#' 
#' @author Villoslada-Blanco, Pablo
#' 
#' @description Server for the Shiny app
#' 
#' @import DT
#' @import shiny
#' @import shinydashboard
#' @import shinyjs
#' @import shinythemes
#' 
#' @param input input of the server
#' @param output output of the server
#' @param session session of the server
#' 
#' @return none

server <- function(input, output, session) {
  
  options(shiny.maxRequestSize = 1000*1024^2)
  
  # Published performance, drawn to match the colour mode chosen by the user
  dark_mode <- reactive({ identical(input$colour_mode, 'dark') })
  output$tumor_img <- renderPlot(accuracy.radar(published.accuracy()$tumor, dark_mode()),
                                 bg = 'transparent', res = 96)
  output$stroma_img <- renderPlot(accuracy.radar(published.accuracy()$stroma, dark_mode()),
                                  bg = 'transparent', res = 96)
  
  output$cnio_logo <- renderImage(list(src = system.file('logos', 'CNIO.jpg', package = 'PDACMOC'),
                                       height = '34px', alt = 'CNIO'), deleteFile = FALSE)
  output$gmeg_logo <- renderImage(list(src = system.file('logos', 'GMEG.png', package = 'PDACMOC'),
                                       height = '40px', alt = 'GMEG'), deleteFile = FALSE)
  
  output$pdacmoc_logo <- renderImage(list(src = system.file('logos', 'PDACMOC_logo.png', package = 'PDACMOC'),
                                          height = '36px', alt = 'PDACMOC'), deleteFile = FALSE)
  output$pdaconsensus_logo <- renderImage(list(src = system.file('logos', 'PDAConsensus_logo.png', package = 'PDACMOC'),
                                               height = '26px', alt = 'PDAConsensus'), deleteFile = FALSE)
  
  # Messages shown in the Run card of the Classify tab. Those about the file stay until
  # another file is uploaded; those of a classification until the next one starts.
  messages <- reactiveVal(data.frame(source = character(0), type = character(0), text = character(0)))
  add.message <- function(text, type = 'info', source = 'file') {
    messages(rbind(isolate(messages()), data.frame(source = source, type = type, text = text)))
  }
  clear.messages <- function(source) {
    current <- isolate(messages())
    messages(current[current$source != source, , drop = FALSE])
  }

  # Messages from the package while reading the file (e.g. averaged genes)
  with.messages <- function(expr) {
    withCallingHandlers(expr, message = function(m) {
      add.message(trimws(conditionMessage(m)))
      invokeRestart('muffleMessage')
    })
  }
  
  # Summary cards before a classification finishes
  summary_ids <- c('overview_collisson', 'overview_moffitt', 'overview_bailey', 'overview_puleo',
                   'overview_chan', 'overview_consensus', 'overview_moffitt_stroma',
                   'overview_maurer_stroma', 'overview_consensus_stroma')
  for (id in summary_ids) {
    local({
      output_id <- id
      output[[output_id]] <- renderUI(div(class = 'summary-empty', 'No classification yet.'))
    })
  }

  uploaded <- reactiveVal(NULL)
  
  observeEvent(input$tsvfile, {
    
    uploaded(NULL)
    shinyjs::disable("runButton")
    clear.messages('file')

    if (!grepl("\\.(tsv|csv|txt)$", input$tsvfile$name, ignore.case = TRUE)) {
      add.message("Invalid file format. Please upload a .tsv, .csv or .txt file.", type = "error")
      return(NULL)
    }

    df <- tryCatch(with.messages(read.counts(input$tsvfile$datapath)),
                   error = function(e) {
                     add.message(paste0("The file could not be read: ", conditionMessage(e)), type = "error")
                     NULL
                   })
    if (is.null(df)) {
      return(NULL)
    }
    
    uploaded(df)
    shinyjs::enable("runButton")
    
  })
  
  # Hide or show the options while a classification runs
  lock.controls <- function(locked) {
    for (id in c('batch', 'gene_id', 'tumor_classifiers', 'stroma', 'stroma_classifiers')) {
      runjs(sprintf('document.getElementById("%s").style.display = "%s";', id, if (locked) 'none' else ''))
    }
    runjs(sprintf('document.getElementById("tsvfile").disabled = %s;', if (locked) 'true' else 'false'))
    if (locked) shinyjs::disable('runButton') else shinyjs::enable('runButton')
  }

  # Fill the Results and Summary tabs and the downloads from a finished classification
  show.results <- function(result, settings) {
    # Results and summary of one classifier (display only; downloads use `result`)
    show.classifier <- function(table_id, summary_id, results, classifier, threshold) {
      output[[table_id]] <- renderDT(classification.table(results, classifier, threshold))
      output[[summary_id]] <- renderUI(subtype.summary(results, classifier, threshold))
    }

    if ('Collisson' %in% settings$tumor_classifiers) {
      show.classifier('collisson', 'overview_collisson', result$Tumor$Collisson, 'Collisson', 50)
    }
    if ('Moffitt' %in% settings$tumor_classifiers) {
      show.classifier('moffitt', 'overview_moffitt', result$Tumor$Moffitt, 'Moffitt', 70)
    }
    if ('Bailey' %in% settings$tumor_classifiers) {
      show.classifier('bailey', 'overview_bailey', result$Tumor$Bailey, 'Bailey', 33)
    }
    if ('Puleo' %in% settings$tumor_classifiers) {
      show.classifier('puleo', 'overview_puleo', result$Tumor$Puleo, 'Puleo', 30)
    }
    if ('Chan-Seng-Yue' %in% settings$tumor_classifiers) {
      show.classifier('chan', 'overview_chan', result$Tumor$`Chan-Seng-Yue`, 'Chan-Seng-Yue', 30)
    }
    if ('PDAConsensus' %in% settings$tumor_classifiers) {
      show.classifier('consensus', 'overview_consensus', result$Tumor$PDAConsensus, 'PDAConsensus', 70)
    }

    if (settings$stroma) {
      output$proportions <- renderDT(datatable(result$Proportions, class = 'compact hover',
                                               options = list(pageLength = 25, dom = 'ftip')))

      if ('Moffitt' %in% settings$stroma_classifiers) {
        show.classifier('moffitt_stroma', 'overview_moffitt_stroma', result$Stroma$Moffitt, 'Stroma Moffitt', 70)
      }
      if ('Maurer' %in% settings$stroma_classifiers) {
        show.classifier('maurer_stroma', 'overview_maurer_stroma', result$Stroma$Maurer, 'Stroma Maurer', 70)
      }
      if ('PDAConsensus' %in% settings$stroma_classifiers) {
        show.classifier('consensus_stroma', 'overview_consensus_stroma', result$Stroma$PDAConsensus,
                        'Stroma PDAConsensus', 70)
      }
    }

    if ('Collisson' %in% settings$tumor_classifiers) {
      shinyjs::enable('downloadCollisson')
      output$downloadCollisson <- downloadHandler(
        filename = function() {paste0(format(Sys.Date(), '%y%m%d'), '_collisson_classification.csv')},
        content = function(file) {
          write.csv(result$Tumor$Collisson, file, row.names = TRUE)
        }
      )
    }
    if ('Moffitt' %in% settings$tumor_classifiers) {
      shinyjs::enable('downloadMoffitt')
      output$downloadMoffitt <- downloadHandler(
        filename = function() {paste0(format(Sys.Date(), '%y%m%d'), '_moffitt_classification.csv')},
        content = function(file) {
          write.csv(result$Tumor$Moffitt, file, row.names = TRUE)
        }
      )
    }
    if ('Bailey' %in% settings$tumor_classifiers) {
      shinyjs::enable('downloadBailey')
      output$downloadBailey <- downloadHandler(
        filename = function() {paste0(format(Sys.Date(), '%y%m%d'), '_bailey_classification.csv')},
        content = function(file) {
          write.csv(result$Tumor$Bailey, file, row.names = TRUE)
        }
      )
    }
    if ('Puleo' %in% settings$tumor_classifiers) {
      shinyjs::enable('downloadPuleo')
      output$downloadPuleo <- downloadHandler(
        filename = function() {paste0(format(Sys.Date(), '%y%m%d'), '_puleo_classification.csv')},
        content = function(file) {
          write.csv(result$Tumor$Puleo, file, row.names = TRUE)
        }
      )
    }
    if ('Chan-Seng-Yue' %in% settings$tumor_classifiers) {
      shinyjs::enable('downloadChan-Seng-Yue')
      output$`downloadChan-Seng-Yue` <- downloadHandler(
        filename = function() {paste0(format(Sys.Date(), '%y%m%d'), '_chan-seng-yue_classification.csv')},
        content = function(file) {
          write.csv(result$Tumor$`Chan-Seng-Yue`, file, row.names = TRUE)
        }
      )
    }
    if ('PDAConsensus' %in% settings$tumor_classifiers) {
      shinyjs::enable('downloadConsensus')
      output$downloadConsensus <- downloadHandler(
        filename = function() {paste0(format(Sys.Date(), '%y%m%d'), '_PDAConsensus_classification.csv')},
        content = function(file) {
          write.csv(result$Tumor$PDAConsensus, file, row.names = TRUE)
        }
      )
    }
    
    if (settings$stroma) {
      shinyjs::enable('downloadProportions')
      output$downloadProportions <- downloadHandler(
        filename = function() {paste0(format(Sys.Date(), '%y%m%d'), '_proportions.csv')},
        content = function(file) {
          write.csv(result$Proportions, file, row.names = TRUE)
        }
      )
      
      if ('Moffitt' %in% settings$stroma_classifiers) {
        shinyjs::enable('downloadMoffittstroma')
        output$downloadMoffittstroma <- downloadHandler(
          filename = function() {paste0(format(Sys.Date(), '%y%m%d'), '_moffitt_stroma_classification.csv')},
          content = function(file) {
            write.csv(result$Stroma$Moffitt, file, row.names = TRUE)
          }
        )
      }
      if ('Maurer' %in% settings$stroma_classifiers) {
        shinyjs::enable('downloadMaurerstroma')
        output$downloadMaurerstroma <- downloadHandler(
          filename = function() {paste0(format(Sys.Date(), '%y%m%d'), '_maurer_stroma_classification.csv')},
          content = function(file) {
            write.csv(result$Stroma$Maurer, file, row.names = TRUE)
          }
        )
      }
      if ('PDAConsensus' %in% settings$stroma_classifiers) {
        shinyjs::enable('downloadConsensusstroma')
        output$downloadConsensusstroma <- downloadHandler(
          filename = function() {paste0(format(Sys.Date(), '%y%m%d'), '_PDAConsensus_stroma_classification.csv')},
          content = function(file) {
            write.csv(result$Stroma$PDAConsensus, file, row.names = TRUE)
          }
        )
      }
    }
  }

  # Background classification: one job per session, polled every second
  job <- reactiveVal(NULL)
  run_status <- reactiveVal(list(state = 'idle'))
  shown_messages <- reactiveVal(0)

  observeEvent(input$runButton, {
    df <- uploaded()
    req(df, is.null(job()))
    settings <- list(batch = input$batch == 'Yes', gene_id = input$gene_id,
                     tumor_classifiers = input$tumor_classifiers,
                     stroma = input$stroma == 'Yes', stroma_classifiers = input$stroma_classifiers)
    args <- list(batch = settings$batch, gene_id = settings$gene_id,
                 classifier = settings$tumor_classifiers)
    if (settings$stroma) {
      args <- c(args, list(stroma = TRUE, stroma_classifier = settings$stroma_classifiers))
    }
    lock.controls(TRUE)
    shown_messages(0)
    clear.messages('run')
    new_job <- background.submit(df, args)
    new_job$settings <- settings
    job(new_job)
    run_status(list(state = 'waiting', job = new_job))
  })

  observe({
    current <- job()
    req(current)
    invalidateLater(1000)
    status <- background.status(current)
    if (!identical(status$job$started, current$started)) {
      status$job$settings <- current$settings
      current <- status$job
      job(current)
    }

    # package messages (imputed genes, averaged duplicates...) in the Run card
    n_shown <- isolate(shown_messages())
    if (length(status$messages) > n_shown) {
      for (m in status$messages[(n_shown + 1):length(status$messages)]) {
        add.message(m, source = 'run')
      }
      shown_messages(length(status$messages))
    }

    if (status$state %in% c('waiting', 'running')) {
      run_status(list(state = status$state, job = current, progress = status$progress))
    } else if (status$state == 'done') {
      result <- background.result(current)
      background.cleanup(current)
      job(NULL)
      lock.controls(FALSE)
      minutes <- as.numeric(difftime(Sys.time(), current$submitted, units = 'mins'))
      show.results(result, current$settings)
      shinyjs::enable('resetButton')
      record.usage(current$n_samples)
      run_status(list(state = 'done', minutes = minutes, n_samples = current$n_samples,
                      progress = status$progress))
      bslib::nav_select('main_nav', 'Results')
    } else {
      message <- background.error(current)
      background.cleanup(current)
      job(NULL)
      lock.controls(FALSE)
      run_status(list(state = 'error', message = message, progress = status$progress))
    }
  })

  observeEvent(input$cancelButton, {
    current <- job()
    req(current)
    background.cancel(current)
    job(NULL)
    lock.controls(FALSE)
    run_status(list(state = 'cancelled', progress = isolate(run_status())$progress))
  })

  observeEvent(input$goToResults, bslib::nav_select('main_nav', 'Results'))

  # Do not leave processes running when the browser tab is closed
  session$onSessionEnded(function() {
    current <- isolate(job())
    if (!is.null(current)) background.cancel(current)
  })

  output$run_panel <- renderUI({
    status <- run_status()
    df <- uploaded()
    estimate <- if (!is.null(df)) {
      estimate.minutes(ncol(df), identical(input$batch, 'Yes'), identical(input$stroma, 'Yes'))
    }
    run.panel(status, df, estimate, messages())
  })

  observeEvent(input$resetButton, {
    session$reload()
  })
  
}
