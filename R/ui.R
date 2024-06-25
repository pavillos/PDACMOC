#' @title UI
#' 
#' @author Villoslada-Blanco, Pablo
#' 
#' @description UI for the Shiny app
#' 
#' @import DT
#' @import shiny
#' @import shinydashboard
#' @import shinyjs
#' @import shinythemes
#' 
#' @return none

# Define the user interface
ui <- fluidPage(
  shinyjs::useShinyjs(),
  theme = shinytheme('cyborg'),
  tags$style(HTML('.some-space { margin-bottom: 50px; }')),
  tags$head(tags$style(HTML('h1 { text-align: center }'))),
  tags$head(tags$style(HTML('.container-fluid { margin-left: 20px; margin-right: 20px; }'))),
  headerPanel(
    div(
      column(style = 'margin-top: 50px; margin-bottom: -230px; ', width = 2, align = 'left', imageOutput('cnio_logo')),
      column(width = 8, align = 'center',
             h1('PDACMolecularOmniClassifier'),
             h2(paste0('v', packageVersion('PDACMOC')))),
      column(style = 'margin-bottom: -180px; ', width = 2, align = 'right', imageOutput('gmeg_logo')))),
  fluidRow(
    style = 'margin-bottom: 150px;',
    column(width = 6, align = 'center', imageOutput('tumor_img')),
    column(width = 6, align = 'center', imageOutput('stroma_img'))
  ),
  fluidRow(
    column(width = 2,
           fileInput('tsvfile', 'Upload a tsv file with genes in rows and samples in columns (max. 1GB)'),
           radioButtons('batch', 'Do you want to apply batch correction?', choices = c('Yes', 'No'), selected = 'Yes'),
           radioButtons('gene_id', 'Select the gene ID type', choices = c('EnsemblID', 'EntrezID', 'GeneSymbol')),
           checkboxGroupInput('tumor_classifiers', 'Select one or more tumor classifiers',
                              choices = list('Collisson' = 'Collisson',
                                             'Moffitt' = 'Moffitt',
                                             'Bailey' = 'Bailey',
                                             'Puleo' = 'Puleo',
                                             'Chan-Seng-Yue' = 'Chan-Seng-Yue',
                                             'PDAConsensus' = 'PDAConsensus'),
                              selected = c('Collisson', 'Moffitt', 'Bailey', 'Puleo', 'Chan-Seng-Yue', 'PDAConsensus')),
           radioButtons('stroma', 'Do you want to include stroma classification?', choices = c('Yes', 'No'), selected = 'Yes'),
           conditionalPanel(
             condition = 'input.stroma == "Yes"',
             checkboxGroupInput('stroma_classifiers', 'Select one or more stroma classifiers',
                                choices = list('Moffitt' = 'Moffitt',
                                               'Maurer' = 'Maurer',
                                               'PDAConsensus' = 'PDAConsensus'),
                                selected = c('Moffitt', 'Maurer', 'PDAConsensus'))
           ),
           actionButton('runButton', 'Run classification', style = 'margin-bottom: 50px;', disabled = TRUE)
    ),
    column(width = 10,
           dropdownMenu(
             type = 'notifications',
             icon = icon('question-circle'),
             badgeStatus = NULL,
             headerText = 'See also: ',
             notificationItem('GitHub', icon = icon('github'),
                              href = 'https://github.com/pavillos/PDACMOC')
           ))
  ),
  h4('Tumor/Stroma proportions'),
  fluidRow(tabPanel('Proportions', DTOutput('proportions')), class = 'some-space'),
  fluidRow(
    column(width = 2, align = 'center', 
           downloadButton('downloadProportions', 'Download tumor/stroma proportions', disabled = TRUE)), class = 'some-space'
  ),
  h4('Tumor classification'),
  fluidRow(tabsetPanel(
             id = 'Tumor',
             tabPanel('Collisson', DTOutput('collisson')),
             tabPanel('Moffitt', DTOutput(('moffitt'))),
             tabPanel('Bailey', DTOutput(('bailey'))),
             tabPanel('Puleo', DTOutput(('puleo'))),
             tabPanel('Chan-Seng-Yue', DTOutput(('chan'))),
             tabPanel('PDAConsensus', DTOutput(('consensus')))
           ),
  ),
  h6('*In red when the probability is below a predefined treshold (70% for Moffitt and PDAConsensus, 50% for Collisson, 33% for Bailey, and 30% for Puleo and Chan-Seng-Yue)', class = 'some-space'),
  fluidRow(
    column(width = 4, align = 'center',
           downloadButton('downloadCollisson', 'Download Collisson classification', disabled = TRUE)),
    column(width = 4, align = 'center',
           downloadButton('downloadMoffitt', 'Download Moffitt classification', disabled = TRUE)),
    column(width = 4, align = 'center',
           downloadButton('downloadBailey', 'Download Bailey classification', disabled = TRUE)), class = 'some-space'
  ),
  fluidRow(
    column(width = 4, align = 'center',
           downloadButton('downloadPuleo', 'Download Puleo classification', disabled = TRUE)),
    column(width = 4, align = 'center',
           downloadButton('downloadChan-Seng-Yue', 'Download Chan-Seng-Yue classification', disabled = TRUE)),
    column(width = 4, align = 'center',
           downloadButton('downloadConsensus', 'Download PDAConsensus classification', disabled = TRUE)), class = 'some-space'
  ),
  h4('Overview tumor classification'),
  fluidRow(
    column(width = 4, align = 'center',
           tabPanel('Overview Collisson', tableOutput('overview_collisson'))),
    column(width = 4, align = 'center',
           tabPanel('Overview Moffitt', tableOutput('overview_moffitt'))),
    column(width = 4, align = 'center',
           tabPanel('Overview Bailey', tableOutput('overview_bailey'))), class = 'some-space'
  ),
  fluidRow(
    column(width = 4, align = 'center',
           tabPanel('Overview Puleo', tableOutput('overview_puleo'))),
    column(width = 4, align = 'center',
           tabPanel('Overview Chan-Sen-Yue', tableOutput('overview_chan'))),
    column(width = 4, align = 'center',
           tabPanel('Overview PDAConsensus', tableOutput('overview_consensus'))), class = 'some-space'
  ),
  h4('Stroma classification'),
  fluidRow(tabsetPanel(
             id = 'Stroma',
             tabPanel('Stroma Moffitt', DTOutput('moffitt_stroma')),
             tabPanel('Stroma Maurer', DTOutput('maurer_stroma')),
             tabPanel('Stroma PDAConsensus', DTOutput('consensus_stroma'))
           )
  ),
  h6('*In red when the probability is below a predefined threshold (70%)', class = 'some-space'),
  fluidRow(
    column(width = 4, align = 'center',
           downloadButton('downloadMoffittstroma', 'Dowload Moffitt stroma classification', disabled = TRUE)),
    column(width = 4, align = 'center',
           downloadButton('downloadMaurerstroma', 'Dowload Maurer stroma classification', disabled = TRUE)),
    column(width = 4, align = 'center',
           downloadButton('downloadConsensusstroma', 'Download PDAConsensus stroma classification', disabled = TRUE)), class = 'some-space'
  ),
  h4('Overview stroma classification'),
  fluidRow(
    column(width = 4, align = 'center',
           tabPanel('Overview stroma Moffit', tableOutput('overview_moffitt_stroma'))),
    column(width = 4, align = 'center',
           tabPanel('Overview stroma Maurer', tableOutput('overview_maurer_stroma'))),
    column(width = 4, align = 'center',
           tabPanel('Overview stroma PDAConsensus', tableOutput('overview_consensus_stroma'))), class = 'some-space'
  ),
  fluidRow(
    style = 'margin-left: 0px; margin-bottom: 30px;',
    actionButton('resetButton', 'Reset app')
  ),
  fluidRow(
    style = 'height: 40 px; margin-left: -50px; margin-bottom: -100px;',
    column(width = 6, align = 'center', imageOutput('pdacmoc_logo'), style = 'margin-top: 20px;'),
    column(width = 6, align = 'center', imageOutput('pdaconsensus_logo'), style = 'margin-top: 110px;')
  )
) 
