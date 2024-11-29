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
  
  output$tumor_img <- renderImage(list(src = system.file('graphs_and_tables', 'balanced_accuracy_tumor.png', package = 'PDACMOC'),
                                       width = '500px'), deleteFile = FALSE)
  output$stroma_img <- renderImage(list(src = system.file('graphs_and_tables', 'balanced_accuracy_stroma.png', package = 'PDACMOC'),
                                        width = '500px'), deleteFile = FALSE)
  
  output$cnio_logo <- renderImage(list(src = system.file('logos', 'CNIO.jpg', package = 'PDACMOC'),
                                       width = '100%'), deleteFile = FALSE)
  output$gmeg_logo <- renderImage(list(src = system.file('logos', 'GMEG.png', package = 'PDACMOC'),
                                       width = '75%'), deleteFile = FALSE)
  
  output$pdacmoc_logo <- renderImage(list(src = system.file('logos', 'PDACMOC_logo.png', package = 'PDACMOC'),
                                          width = '75%'), deleteFile = FALSE)
  output$pdaconsensus_logo <- renderImage(list(src = system.file('logos', 'PDAConsensus_logo.png', package = 'PDACMOC'),
                                               width = '80%'), deleteFile = FALSE)
  
  observeEvent(input$tsvfile, {
    
    if (!grepl("\\.tsv$", input$tsvfile$name, ignore.case = TRUE)) {
      showNotification("Invalid file format. Please upload a tsv file.", type = "error")
      shinyjs::disable("runButton")
      return(NULL)
    }
    
    shinyjs::enable("runButton")
    df <- read.csv(input$tsvfile$datapath, sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)
    
    observeEvent(input$runButton, {
      
      runjs('document.getElementById("tsvfile").disabled = "true";')
      runjs('document.getElementById("batch").style.display = "none";')
      runjs('document.getElementById("gene_id").style.display = "none";')
      runjs('document.getElementById("tumor_classifiers").style.display = "none";')
      runjs('document.getElementById("stroma").style.display = "none";')
      runjs('document.getElementById("stroma_classifiers").style.display = "none";')
      runjs('document.getElementById("runButton").disabled = "true";')
            
      start_time <- Sys.time()
      if (input$stroma == 'Yes') {
        if (input$batch == 'Yes') {batch_opt <- TRUE} else if (input$batch == 'No') {batch_opt <- FALSE}
        result <- omni.classify(df, batch = batch_opt, gene_id = input$gene_id,
                                classifier = input$tumor_classifiers,
                                stroma = TRUE,
                                stroma_classifier = input$stroma_classifiers)
      } else if (input$stroma == 'No') {
        if (input$batch == 'Yes') {batch_opt <- TRUE} else if (input$batch == 'No') {batch_opt <- FALSE}
        result <- omni.classify(df, batch = batch_opt, gene_id = input$gene_id,
                                classifier = input$tumor_classifiers)
      }

      end_time <- Sys.time()
      time_in_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
      if (time_in_seconds >= 60) {
        time_in_minutes <- round(time_in_seconds / 60, 2)
        showNotification(paste0("Classification has finished in ", time_in_minutes, " minutes"), duration = 60)
      } else {
        showNotification(paste0("Classification has finished in ", round(time_in_seconds, 2), " seconds"), duration = 60)
      }
      
      if ('Collisson' %in% input$tumor_classifiers) {
        output$collisson <- renderDT(datatable(result$Tumor$Collisson,
                                               options = list(columnDefs = list(list(className = 'dt-center', targets = 1:2)))) %>%
                                       formatStyle('Probability',
                                                   backgroundColor = styleInterval(49, c('red', 'default'))))
        matrix_collisson <- as.matrix(table(result$Tumor$Collisson$Predicted.subtype))
        colnames(matrix_collisson) <- 'Collisson'
        output$overview_collisson <- renderTable(matrix_collisson, rownames = TRUE)
      }
      if ('Moffitt' %in% input$tumor_classifiers) {
        output$moffitt <- renderDT(datatable(result$Tumor$Moffitt,
                                             options = list(columnDefs = list(list(className = 'dt-center', targets = 1:2)))) %>%
                                     formatStyle('Probability',
                                                 backgroundColor = styleInterval(69, c('red', 'default'))))
        matrix_moffitt <- as.matrix(table(result$Tumor$Moffitt$Predicted.subtype))
        colnames(matrix_moffitt) <- 'Moffitt'
        output$overview_moffitt <- renderTable(matrix_moffitt, rownames = TRUE)
      }
      if ('Bailey' %in% input$tumor_classifiers) {
        output$bailey <- renderDT(datatable(result$Tumor$Bailey,
                                            options = list(columnDefs = list(list(className = 'dt-center', targets = 1:2)))) %>%
                                    formatStyle('Probability',
                                                backgroundColor = styleInterval(32, c('red', 'default'))))
        matrix_bailey <- as.matrix(table(result$Tumor$Bailey$Predicted.subtype))
        colnames(matrix_bailey) <- 'Bailey'
        output$overview_bailey <- renderTable(matrix_bailey, rownames = TRUE)
      }
      if ('Puleo' %in% input$tumor_classifiers) {
        output$puleo <- renderDT(datatable(result$Tumor$Puleo,
                                           options = list(columnDefs = list(list(className = 'dt-center', targets = 1:2)))) %>%
                                   formatStyle('Probability',
                                               backgroundColor = styleInterval(29, c('red', 'default'))))
        matrix_puleo <- as.matrix(table(result$Tumor$Puleo$Predicted.subtype))
        colnames(matrix_puleo) <- 'Puleo'
        output$overview_puleo <- renderTable(matrix_puleo, rownames = TRUE)
      }
      if ('Chan-Seng-Yue' %in% input$tumor_classifiers) {
        output$chan <- renderDT(datatable(result$Tumor$`Chan-Seng-Yue`,
                                          options = list(columnDefs = list(list(className = 'dt-center', targets = 1:2)))) %>%
                                  formatStyle('Probability',
                                              backgroundColor = styleInterval(29, c('red', 'default'))))
        matrix_chan <- as.matrix(table(result$Tumor$`Chan-Seng-Yue`$Predicted.subtype))
        colnames(matrix_chan) <- 'Chan-Seng-Yue'
        output$overview_chan <- renderTable(matrix_chan, rownames = TRUE)
      }
      if ('PDAConsensus' %in% input$tumor_classifiers) {
        output$consensus <- renderDT(datatable(result$Tumor$PDAConsensus,
                                               options = list(columnDefs = list(list(className = 'dt-center', targets = 1:3)))) %>%
                                     formatStyle('Probability',
                                                   backgroundColor = styleInterval(69, c('red', 'default'))) %>%
                                     formatRound('NonClassicalScore', 4) %>%
                                     formatStyle('NonClassicalScore',
                                                 background = styleColorBar(result$Tumor$PDAConsensus$NonClassicalScore, '#E02438'),
                                                 backgroundSize = '95% 50%',
                                                 backgroundRepeat = 'no-repeat',
                                                 backgroundPosition = 'right'))
        matrix_consensus <- as.matrix(table(result$Tumor$PDAConsensus$Predicted.subtype))
        colnames(matrix_consensus) <- 'PDAConsensus'
        output$overview_consensus <- renderTable(matrix_consensus, rownames = TRUE)
      }
      
      if (input$stroma == 'Yes') {
        output$proportions <- renderDT(datatable(result$Proportions,
                                                 options = list(columnDefs = list(list(className = 'dt-center', targets = 1:7)))))
        
        if ('Moffitt' %in% input$stroma_classifiers) {
          output$moffitt_stroma <- renderDT(datatable(result$Stroma$Moffitt,
                                                      options = list(columnDefs = list(list(className = 'dt-center', targets = 1:2)))) %>%
                                                      formatStyle('Probability',
                                                                  backgroundColor = styleInterval(69, c('red', 'default'))))
          matrix_moffitt_stroma <- as.matrix(table(result$Stroma$Moffitt$Predicted.subtype))
          colnames(matrix_moffitt_stroma) <- 'Moffitt'
          output$overview_moffitt_stroma <- renderTable(matrix_moffitt_stroma, rownames = TRUE)
        }
        if ('Maurer' %in% input$stroma_classifiers) {
          output$maurer_stroma <- renderDT(datatable(result$Stroma$Maurer,
                                                     options = list(columnDefs = list(list(className = 'dt-center', targets = 1:2)))) %>%
                                             formatStyle('Probability',
                                                         backgroundColor = styleInterval(69, c('red', 'default'))))
          matrix_maurer <- as.matrix(table(result$Stroma$Maurer$Predicted.subtype))
          colnames(matrix_maurer) <- 'Maurer'
          output$overview_maurer_stroma <- renderTable(matrix_maurer, rownames = TRUE)
        }
        if ('PDAConsensus' %in% input$stroma_classifiers) {
          output$consensus_stroma <- renderDT(datatable(result$Stroma$PDAConsensus,
                                                        options = list(columnDefs = list(list(className = 'dt-center', targets = 1:3)))) %>%
                                              formatStyle('Probability',
                                                          backgroundColor = styleInterval(69, c('red', 'default'))) %>%
                                              formatRound('ActivatedECMScore', 4) %>%
                                              formatStyle('ActivatedECMScore',
                                                          background = styleColorBar(result$Stroma$PDAConsensus$ActivatedECMScore, '#CE18A2'),
                                                          backgroundSize = '95% 50%',
                                                          backgroundRepeat = 'no-repeat',
                                                          backgroundPosition = 'right'))
          matrix_consensus_stroma <- as.matrix(table(result$Stroma$PDAConsensus$Predicted.subtype))
          colnames(matrix_consensus_stroma) <- 'PDAConsensus'
          output$overview_consensus_stroma <- renderTable(matrix_consensus_stroma, rownames = TRUE)
        }
      }
      
      if ('Collisson' %in% input$tumor_classifiers) {
        shinyjs::enable('downloadCollisson')
        output$downloadCollisson <- downloadHandler(
          filename = function() {paste0(format(Sys.Date(), '%y%m%d'), '_collisson_classification.csv')},
          content = function(file) {
            write.csv(result$Tumor$Collisson, file, row.names = TRUE)
          }
        )
      }
      if ('Moffitt' %in% input$tumor_classifiers) {
        shinyjs::enable('downloadMoffitt')
        output$downloadMoffitt <- downloadHandler(
          filename = function() {paste0(format(Sys.Date(), '%y%m%d'), '_moffitt_classification.csv')},
          content = function(file) {
            write.csv(result$Tumor$Moffitt, file, row.names = TRUE)
          }
        )
      }
      if ('Bailey' %in% input$tumor_classifiers) {
        shinyjs::enable('downloadBailey')
        output$downloadBailey <- downloadHandler(
          filename = function() {paste0(format(Sys.Date(), '%y%m%d'), '_bailey_classification.csv')},
          content = function(file) {
            write.csv(result$Tumor$Bailey, file, row.names = TRUE)
          }
        )
      }
      if ('Puleo' %in% input$tumor_classifiers) {
        shinyjs::enable('downloadPuleo')
        output$downloadPuleo <- downloadHandler(
          filename = function() {paste0(format(Sys.Date(), '%y%m%d'), '_puleo_classification.csv')},
          content = function(file) {
            write.csv(result$Tumor$Puleo, file, row.names = TRUE)
          }
        )
      }
      if ('Chan-Seng-Yue' %in% input$tumor_classifiers) {
        shinyjs::enable('downloadChan-Seng-Yue')
        output$`downloadChan-Seng-Yue` <- downloadHandler(
          filename = function() {paste0(format(Sys.Date(), '%y%m%d'), '_chan-seng-yue_classification.csv')},
          content = function(file) {
            write.csv(result$Tumor$`Chan-Seng-Yue`, file, row.names = TRUE)
          }
        )
      }
      if ('PDAConsensus' %in% input$tumor_classifiers) {
        shinyjs::enable('downloadConsensus')
        output$downloadConsensus <- downloadHandler(
          filename = function() {paste0(format(Sys.Date(), '%y%m%d'), '_PDAConsensus_classification.csv')},
          content = function(file) {
            write.csv(result$Tumor$PDAConsensus, file, row.names = TRUE)
          }
        )
      }
      
      if (input$stroma == 'Yes') {
        shinyjs::enable('downloadProportions')
        output$downloadProportions <- downloadHandler(
          filename = function() {paste0(format(Sys.Date(), '%y%m%d'), '_proportions.csv')},
          content = function(file) {
            write.csv(result$Proportions, file, row.names = TRUE)
          }
        )
        
        if ('Moffitt' %in% input$stroma_classifiers) {
          shinyjs::enable('downloadMoffittstroma')
          output$downloadMoffittstroma <- downloadHandler(
            filename = function() {paste0(format(Sys.Date(), '%y%m%d'), '_moffitt_stroma_classification.csv')},
            content = function(file) {
              write.csv(result$Stroma$Moffitt, file, row.names = TRUE)
            }
          )
        }
        if ('Maurer' %in% input$stroma_classifiers) {
          shinyjs::enable('downloadMaurerstroma')
          output$downloadMaurerstroma <- downloadHandler(
            filename = function() {paste0(format(Sys.Date(), '%y%m%d'), '_maurer_stroma_classification.csv')},
            content = function(file) {
              write.csv(result$Stroma$Maurer, file, row.names = TRUE)
            }
          )
        }
        if ('PDAConsensus' %in% input$stroma_classifiers) {
          shinyjs::enable('downloadConsensusstroma')
          output$downloadConsensusstroma <- downloadHandler(
            filename = function() {paste0(format(Sys.Date(), '%y%m%d'), '_PDAConsensus_stroma_classification.csv')},
            content = function(file) {
              write.csv(result$Stroma$Consensus, file, row.names = TRUE)
            }
          )
        }
      }
      
    })

    observeEvent(input$resetButton, {
      session$reload()
    })
    
  })
}
