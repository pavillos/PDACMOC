# Version shown in the navbar and the browser tab
app_version <- paste0('v', utils::packageVersion('PDACMOC'))

# A result tab: table, note and its download button
result.tab <- function(title, output_id, download_id, download_label, note = NULL) {
  bslib::nav_panel(
    title,
    div(class = 'result-body',
        DTOutput(output_id, height = 'auto'),
        if (!is.null(note)) div(class = 'result-note', note),
        div(class = 'result-actions', downloadButton(download_id, download_label, class = 'btn-primary', disabled = TRUE)))
  )
}

section.title <- function(...) div(class = 'section-title', ...)

# Metrics shown in Help, from the same sources as the README badges
gist_json <- 'https://gist.githubusercontent.com/pavillos/e55b1802a1ed7b3b815189d7e0c0b802/raw/traffic.json'
zenodo_json <- 'https://zenodo.org/api/records/17019896'
usage_json <- 'https://pdacmoc.cnio.es/stats/usage.json'

metric.badge <- function(alt, json, query, label, colour, href) {
  src <- paste0('https://img.shields.io/badge/dynamic/json?url=', utils::URLencode(json, reserved = TRUE),
                '&query=', utils::URLencode(query, reserved = TRUE),
                '&label=', utils::URLencode(label, reserved = TRUE), '&color=', colour)
  tags$a(href = href, target = '_blank', tags$img(alt = alt, src = src))
}

metric.row <- function(name, ...) div(class = 'metric-row', span(class = 'metric-name', name), ...)

threshold.note <- function(threshold) {
  paste0('The probability is how sure the classifier is of the subtype it assigns. Below ',
         threshold, '% the sample is flagged as low confidence: read it as undetermined rather ',
         'than as belonging to that subtype. Subtype colours are the ones used in the paper.')
}

# Define the user interface
#' @title UI
#'
#' @author Villoslada-Blanco, Pablo
#'
#' @description UI for the Shiny app
#'
#' @import DT
#' @import shiny
#' @import shinyjs
#'
#' @return none

ui <- bslib::page_navbar(
  id = 'main_nav',
  title = div(
    class = 'app-title',
    span(class = 'app-name', 'PDACMolecularOmniClassifier'),
    span(class = 'app-version', app_version)
  ),
  window_title = paste('PDACMOC', app_version),
  fillable = FALSE,
  theme = bslib::bs_theme(
    version = 5,
    primary = '#6B2E99',
    secondary = '#185C84',
    base_font = bslib::font_google('IBM Plex Sans'),
    code_font = bslib::font_google('IBM Plex Mono'),
    heading_font = bslib::font_google('Source Serif 4')
  ),
  header = tagList(
    shinyjs::useShinyjs(),
    tags$style(HTML("
      .app-title { display: flex; align-items: baseline; gap: .5rem; flex-wrap: wrap; }
      .app-name { font-size: 1.45rem; font-weight: 600; color: #6B2E99; letter-spacing: -.01em; }
      .app-version { font-family: var(--bs-font-monospace); font-size: .8rem; color: var(--bs-secondary-color); }
      [data-bs-theme='dark'] .app-name { color: #c69bea; }
      .navbar { padding-block: .55rem; }
      .navbar > .container-fluid { align-items: center; }
      .navbar-brand { display: flex; align-items: center; margin-right: 2.5rem !important; padding-block: 0; }
      .navbar-nav { align-items: center; }
      .navbar .nav-link { padding-inline: .85rem; }
      .app-title { align-items: baseline; line-height: 1.1; }

      /* cards and tables take the height of their content, they do not stretch */
      .card, .card-body, .tab-content, .tab-pane, .html-fill-container { flex: 0 0 auto !important; }
      .card-body > .html-fill-item { flex: 0 0 auto !important; }
      .card-body > .datatables { height: auto !important; }
      .card-body { min-height: 0 !important; }

      /* the footer always sits at the bottom of the window */
      body > .container-fluid { display: flex !important; flex-direction: column !important; min-height: calc(100vh - 72px) !important; }
      body > .container-fluid > .tab-content { flex: 0 0 auto; }
      .app-footer { margin-top: auto; }
      .navbar-logos { display: flex; align-items: center; gap: .9rem; margin-right: .75rem; }
      .navbar-logos .logo-chip { background: #ffffff; border-radius: .35rem; padding: .3rem .5rem;
        display: flex; align-items: center; }
      .navbar-logos img { display: block; width: auto; }

      .card { box-shadow: none; }
      .card-header { font-weight: 600; font-size: .95rem; }
      .section-title { font-size: .72rem; font-weight: 600; letter-spacing: .06em;
        text-transform: uppercase; color: var(--bs-secondary-color); margin: 0 0 .35rem; }

      .sidebar .form-group, .sidebar .shiny-input-container { margin-bottom: .35rem; }
      .helper { font-size: .78rem; color: var(--bs-secondary-color); margin: 0 0 1rem; line-height: 1.35; }
      .sidebar .btn { width: 100%; }
      .btn-file { background: var(--bs-primary); border-color: var(--bs-primary); color: #fff; width: auto !important; }
      .btn-file:hover { filter: brightness(.9); color: #fff; }

      .empty-state { display: flex; flex-direction: column; gap: .4rem; padding: .5rem 0 .25rem; }
      .empty-state .step { display: flex; gap: .6rem; align-items: baseline; font-size: .9rem;
        color: var(--bs-secondary-color); }
      .empty-state .num { font-family: var(--bs-font-monospace); font-size: .75rem;
        color: var(--bs-primary); font-weight: 600; }

      .run-head { display: flex; justify-content: space-between; align-items: baseline; gap: .75rem; flex-wrap: wrap; }
      .run-title { font-weight: 600; }
      .run-time { font-family: var(--bs-font-monospace); font-size: .82rem; color: var(--bs-secondary-color); }
      .run-bar { height: .5rem; border-radius: 999px; background: var(--bs-tertiary-bg); overflow: hidden; margin: .75rem 0; }
      .run-bar > div { height: 100%; border-radius: 999px; background: var(--bs-primary); transition: width .6s ease; }
      .run-steps { list-style: none; padding: 0; margin: 0 0 .5rem; display: grid; gap: .3rem; font-size: .88rem; }
      .run-steps li { display: flex; align-items: center; gap: .55rem; color: var(--bs-secondary-color); }
      .run-steps li.now { color: var(--bs-primary); font-weight: 500; }
      .run-steps .dot { width: .45rem; height: .45rem; border-radius: 50%; background: var(--bs-border-color); flex: none; }
      .run-steps li.done .dot { background: #2c6a4a; }
      .run-steps li.now .dot { background: var(--bs-primary); }
      .run-steps li.stopped .dot { background: #9b1c1c; }
      .run-actions { display: flex; justify-content: space-between; align-items: center; gap: .75rem; flex-wrap: wrap; margin-top: .5rem; }
      .run-actions .helper { margin: 0; }
      .run-banner { border-radius: .4rem; padding: .55rem .8rem; margin-bottom: .75rem; font-size: .9rem;
        background: var(--bs-tertiary-bg); display: flex; gap: .75rem; flex-wrap: wrap; align-items: baseline; }
      .run-banner.done { background: #e7f2ec; color: #2c6a4a; }
      .run-banner.error { background: #fdeaea; color: #9b1c1c; }
      [data-bs-theme='dark'] .run-banner.done { background: #1f3329; color: #9ad0b1; }
      [data-bs-theme='dark'] .run-banner.error { background: #4a2a30; color: #f0a9ac; }
      .run-messages { margin-top: 1rem; padding-top: .75rem; border-top: 1px solid var(--bs-border-color); }
      .run-messages-title { font-weight: 600; font-size: .9rem; margin-bottom: .4rem; }
      .run-messages ul { list-style: none; padding: 0; margin: 0; display: grid; gap: .35rem; font-size: .88rem; }
      .run-messages li { display: flex; gap: .55rem; align-items: baseline; color: var(--bs-secondary-color); }
      .run-messages li.error { color: #9b1c1c; }
      [data-bs-theme='dark'] .run-messages li.error { color: #f0a9ac; }
      .run-messages .icon { flex: none; width: 1.1rem; height: 1.1rem; border-radius: 50%; font-size: .72rem;
        font-weight: 700; display: inline-flex; align-items: center; justify-content: center;
        border: 1.5px solid currentColor; font-style: italic; font-family: Georgia, serif; }
      .run-messages li.error .icon { font-style: normal; font-family: inherit; }

      .acc-list { display: grid; gap: .45rem; }
      .acc-row { display: grid; grid-template-columns: 9.5rem 1fr 4rem; gap: .6rem; align-items: center; }
      .acc-name { font-size: .85rem; color: var(--bs-secondary-color); }
      .acc-track { height: .5rem; border-radius: 999px; background: var(--bs-tertiary-bg); overflow: hidden; }
      .acc-fill { height: 100%; border-radius: 999px; background: var(--bs-primary); }
      .acc-value { font-family: var(--bs-font-monospace); font-size: .8rem; text-align: right;
        font-variant-numeric: tabular-nums; }
      @media (max-width: 575px) { .acc-row { grid-template-columns: 1fr 3.5rem; } .acc-track { display: none; } }

      .result-body { padding-top: .35rem; }
      .result-body table.dataTable { width: auto !important; min-width: min(100%, 44rem); margin-left: 0 !important; }
      .result-body table.dataTable th, .result-body table.dataTable td { padding: .55rem 1.4rem .55rem .9rem; }
      .result-body table.dataTable thead th { font-size: .72rem; font-weight: 600; letter-spacing: .06em;
        text-transform: uppercase; color: var(--bs-secondary-color); background: var(--bs-tertiary-bg);
        border-bottom: 1px solid var(--bs-border-color); }
      .result-body table.dataTable tbody td { border-top: 1px solid var(--bs-border-color); }
      .result-body table.dataTable tbody td:first-child { font-family: var(--bs-font-monospace); font-size: .82rem; }
      .result-body table.dataTable tbody td.dt-right { font-family: var(--bs-font-monospace); }
      .result-body .dataTables_wrapper { max-width: 60rem; }
      .summary-empty { font-size: .85rem; color: var(--bs-secondary-color); }
      .result-note { font-size: .82rem; color: var(--bs-secondary-color); margin-top: .75rem; }
      .result-actions { margin-top: .75rem; }

      .subtype { display: inline-flex; align-items: center; gap: .45rem; white-space: nowrap;
        padding: .12rem .6rem .12rem .45rem; border-radius: 999px; font-size: .85rem;
        border: 1px solid color-mix(in srgb, var(--c) 55%, var(--bs-border-color));
        background: color-mix(in srgb, var(--c) 18%, transparent); }
      .subtype i { width: .6rem; height: .6rem; border-radius: 2px; display: inline-block;
        background: var(--c); box-shadow: inset 0 0 0 1px rgba(0,0,0,.18); }
      .confidence { display: inline-block; font-size: .75rem; font-weight: 600; border-radius: 999px;
        padding: .08rem .55rem; }
      .confidence.high { color: #2c6a4a; background: #e7f2ec; }
      .confidence.low { color: #9b1c1c; background: #fdeaea; }
      [data-bs-theme='dark'] .confidence.high { color: #9ad0b1; background: #1f3329; }
      [data-bs-theme='dark'] .confidence.low { color: #f0a9ac; background: #4a2a30; }
      .low-confidence { display: inline-block; font-size: .72rem; font-weight: 600; letter-spacing: .02em;
        color: #9b1c1c; background: #fdeaea; border-radius: 999px; padding: .05rem .5rem; }
      [data-bs-theme='dark'] .low-confidence { color: #f0a9ac; background: #4a2a30; }

      .subtype-summary { display: grid; gap: .5rem; }
      .summary-meta { display: flex; gap: .6rem; align-items: center; font-size: .8rem;
        color: var(--bs-secondary-color); }
      .segbar { display: flex; height: 1.6rem; border-radius: .35rem; overflow: hidden;
        background: var(--bs-tertiary-bg); gap: 2px; }
      .segbar .seg { display: flex; align-items: center; justify-content: center; font-size: .75rem;
        font-family: var(--bs-font-monospace); font-weight: 600; min-width: 2px; }
      .summary-legend { display: flex; flex-wrap: wrap; gap: .3rem .9rem; font-size: .78rem;
        color: var(--bs-secondary-color); }
      .summary-legend .legend-item { display: inline-flex; align-items: center; gap: .35rem; }
      .summary-legend i { width: .65rem; height: .65rem; border-radius: 2px; display: inline-block;
        box-shadow: inset 0 0 0 1px rgba(0,0,0,.18); }

      .help-block h3 { font-size: 1rem; margin-top: 1.1rem; }
      .help-block h3:first-child { margin-top: 0; }
      .help-block p, .help-block li { font-size: .9rem; color: var(--bs-secondary-color); }

      .metric-row { display: flex; align-items: center; gap: .6rem; flex-wrap: wrap; margin-bottom: .55rem; }
      .metric-name { font-size: .85rem; font-weight: 600; min-width: 5.5rem; }
      .metric-row img { display: block; height: 20px; }
      .usage-line { font-size: .85rem; color: var(--bs-secondary-color); margin-top: .8rem; }
      .download-all { display: flex; justify-content: flex-end; align-items: center; gap: .9rem; flex-wrap: wrap; margin: .2rem 0 .8rem; }
      .download-all .helper { margin: 0; }

      .app-footer { border-top: 1px solid var(--bs-border-color); margin-top: auto; padding: 1.1rem 0 1.4rem; }
      .app-footer::before { content: ''; display: block; height: 0; }
      body > .container-fluid > .tab-content { margin-bottom: 1.5rem; }
      .logo-row { display: flex; align-items: center; justify-content: center; gap: 2.5rem; flex-wrap: wrap; }
      .logo-row .shiny-image-output { display: flex; align-items: center; }
      .logo-row img { display: block; width: auto; }

      [data-bs-theme='dark'] {
        --bs-body-bg: #2f2939; --bs-body-color: #efe9f4;
        --bs-secondary-color: #b2a7bd; --bs-border-color: #4b4159;
        --bs-tertiary-bg: #3f374d; --bs-emphasis-color: #ffffff;
      }
      [data-bs-theme='dark'] .card, [data-bs-theme='dark'] .navbar, [data-bs-theme='dark'] .sidebar {
        background-color: #38314595 !important; }
      [data-bs-theme='dark'] .card { border-color: #4b4159; }
    ")),
    tags$script(HTML("$(function(){ $('.acc-list').closest('.card').addClass('acc-card'); });"))
  ),

  bslib::nav_panel(
    'Classify',
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        width = 320,
        title = 'Your samples',
        fileInput('tsvfile', 'Expression file',
                  accept = c('.tsv', '.csv', 'text/tab-separated-values', 'text/csv')),
        div(class = 'helper', 'Raw counts with genes in rows and samples in columns, up to 1 GB.'),
        radioButtons('gene_id', 'Gene ID type', choices = c('EnsemblID', 'EntrezID', 'GeneSymbol')),
        div(class = 'helper', 'The identifiers in the first column of your file.'),
        radioButtons('batch', 'Batch correction', choices = c('Yes', 'No'), selected = 'Yes'),
        div(class = 'helper',
            'Corrects your samples against the 514 training samples. Use the same setting for the whole study.'),
        checkboxGroupInput(
          'tumor_classifiers', 'Tumor classifiers',
          choiceNames = list('Collisson', 'Moffitt', 'Bailey', 'Puleo', 'Chan-Seng-Yue', pdaconsensus.label()),
          choiceValues = c('Collisson', 'Moffitt', 'Bailey', 'Puleo', 'Chan-Seng-Yue', 'PDAConsensus'),
          selected = c('Collisson', 'Moffitt', 'Bailey', 'Puleo', 'Chan-Seng-Yue', 'PDAConsensus')),
        radioButtons('stroma', 'Stroma classification', choices = c('Yes', 'No'), selected = 'Yes'),
        conditionalPanel(
          condition = 'input.stroma == "Yes"',
          checkboxGroupInput(
            'stroma_classifiers', 'Stroma classifiers',
            choiceNames = list('Moffitt', 'Maurer', pdaconsensus.label()),
            choiceValues = c('Moffitt', 'Maurer', 'PDAConsensus'),
            selected = c('Moffitt', 'Maurer', 'PDAConsensus'))),
        div(class = 'helper', 'Stroma needs virtual microdissection, which roughly doubles the time.'),
        actionButton('runButton', 'Run classification', class = 'btn-primary', disabled = TRUE),
        actionButton('resetButton', 'New classification', class = 'btn-outline-primary', disabled = TRUE)
      ),
      bslib::card(
        bslib::card_header('Run'),
        bslib::card_body(uiOutput('run_panel'))
      )
    )
  ),

  bslib::nav_panel(
    'Results',
    div(class = 'download-all',
        span(class = 'helper', 'All the tables and the summary figure in one zip file.'),
        downloadButton('downloadAll', 'Download all (zip)', class = 'btn-primary', disabled = TRUE)),
    bslib::navset_card_tab(
      full_screen = FALSE,
      result.tab('Collisson', 'collisson', 'downloadCollisson',
                 'Download Collisson classification', threshold.note(50)),
      result.tab('Moffitt', 'moffitt', 'downloadMoffitt',
                 'Download Moffitt classification', threshold.note(70)),
      result.tab('Bailey', 'bailey', 'downloadBailey',
                 'Download Bailey classification', threshold.note(33)),
      result.tab('Puleo', 'puleo', 'downloadPuleo',
                 'Download Puleo classification', threshold.note(30)),
      result.tab('Chan-Seng-Yue', 'chan', 'downloadChan-Seng-Yue',
                 'Download Chan-Seng-Yue classification', threshold.note(30)),
      result.tab(pdaconsensus.label(), 'consensus', 'downloadConsensus',
                 tagList('Download ', pdaconsensus.label('classification')), threshold.note(70)),
      bslib::nav_spacer(),
      result.tab('Stroma Moffitt', 'moffitt_stroma', 'downloadMoffittstroma',
                 'Download Moffitt stroma classification', threshold.note(70)),
      result.tab('Stroma Maurer', 'maurer_stroma', 'downloadMaurerstroma',
                 'Download Maurer stroma classification', threshold.note(70)),
      result.tab(tagList('Stroma ', pdaconsensus.label()), 'consensus_stroma', 'downloadConsensusstroma',
                 tagList('Download ', pdaconsensus.label('stroma classification')), threshold.note(70)),
    ),
    bslib::card(
      fill = FALSE,
      bslib::card_header('Tumor and stroma proportions'),
      bslib::card_body(fill = FALSE,
        div(class = 'helper', style = 'margin-bottom:.5rem',
            'Estimated fraction of epithelium (E), stroma (S) and other cells (O) of each sample, from the virtual microdissection. Only filled in when stroma classification is included.'),
        DTOutput('proportions', height = 'auto'),
        div(class = 'result-actions', downloadButton('downloadProportions', 'Download tumor/stroma proportions', class = 'btn-primary', disabled = TRUE))
      )
    )
  ),

  bslib::nav_panel(
    'Summary',
    section.title('Tumor classification'),
    bslib::layout_columns(
      col_widths = c(4, 4, 4), fill = FALSE,
      bslib::card(bslib::card_header('Collisson'), bslib::card_body(uiOutput('overview_collisson'))),
      bslib::card(bslib::card_header('Moffitt'), bslib::card_body(uiOutput('overview_moffitt'))),
      bslib::card(bslib::card_header('Bailey'), bslib::card_body(uiOutput('overview_bailey'))),
      bslib::card(bslib::card_header('Puleo'), bslib::card_body(uiOutput('overview_puleo'))),
      bslib::card(bslib::card_header('Chan-Seng-Yue'), bslib::card_body(uiOutput('overview_chan'))),
      bslib::card(bslib::card_header(pdaconsensus.label()), bslib::card_body(uiOutput('overview_consensus')))
    ),
    section.title('Stroma classification'),
    bslib::layout_columns(
      col_widths = c(4, 4, 4), fill = FALSE,
      bslib::card(bslib::card_header('Moffitt'), bslib::card_body(uiOutput('overview_moffitt_stroma'))),
      bslib::card(bslib::card_header('Maurer'), bslib::card_body(uiOutput('overview_maurer_stroma'))),
      bslib::card(bslib::card_header(pdaconsensus.label()), bslib::card_body(uiOutput('overview_consensus_stroma')))
    ),
    div(class = 'helper', 'Number of samples assigned to each subtype, in the colours of the paper. Filled in when a classification finishes.'),
    div(class = 'result-actions',
        downloadButton('downloadSummary', 'Download summary (PNG)', class = 'btn-primary', disabled = TRUE))
  ),

  bslib::nav_panel(
    'Performance',
    div(class = 'helper',
        paste('Mean balanced accuracy of each classifier, as published (Genome Medicine 2025, Tables 4 and 6),',
              'estimated by repeated 10-times 10-fold cross-validation on the training cohort of 514 samples.',
              'The dashed shape is what a classifier that always predicts the most frequent class would reach.')),
    bslib::layout_columns(
      col_widths = c(6, 6), fill = FALSE,
      bslib::card(
        bslib::card_header('Tumor classifiers'),
        bslib::card_body(plotOutput('tumor_img', height = '420px'))),
      bslib::card(
        bslib::card_header('Stroma classifiers'),
        bslib::card_body(plotOutput('stroma_img', height = '420px')))
    )
  ),

  bslib::nav_panel(
    'Help',
    bslib::layout_columns(
      col_widths = c(6, 6), fill = FALSE,
      bslib::card(bslib::card_header('Input file'), bslib::card_body(
        div(class = 'help-block',
            h3('Format'),
            p('Raw counts in a tab-separated (.tsv) or comma-separated (.csv) file: gene IDs in the first ',
              'column, one column per sample, up to 1 GB.'),
            h3('Gene IDs'),
            p('Ensembl IDs may include the version, which is removed. Rows that end up sharing a gene ID ',
              'are averaged.'),
            h3('Batch correction'),
            p('Your samples are corrected together with the 514 training samples, treating yours as one ',
              'batch. The result depends on which samples are uploaded together, so use the same setting ',
              'for the whole study.'),
            h3('How long it takes'),
            p('About 5 minutes for 5 samples and 20 minutes for 100, with batch correction. Stroma ',
              'classification roughly doubles the time.'))
      )),
      bslib::card(bslib::card_header('Reading the results'), bslib::card_body(
        div(class = 'help-block',
            h3('Probability'),
            p('How sure the classifier is of the subtype it assigns. Below the threshold of its ',
              'classifier the sample is flagged as low confidence and should be read as undetermined.'),
            tags$ul(
              tags$li('70% for Moffitt, the stroma classifiers and ', pdaconsensus.label()),
              tags$li('50% for Collisson'),
              tags$li('33% for Bailey'),
              tags$li('30% for Puleo and Chan-Seng-Yue')),
            h3('Scores'),
            p('NonClassicalScore and ActivatedECMScore are the probability of belonging to the ',
              'non-classical and activated-ECM classes, as continuous values between 0 and 1.'))
      )),
      bslib::card(bslib::card_header('About'), bslib::card_body(
        div(class = 'help-block',
            p(strong('Version: '), app_version),
            p(strong('Cite the method: '),
              'Villoslada-Blanco P, Alonso L, Sabroso-Lasa S, Maquedano M, Estudillo L, Real FX, ',
              'López de Maturana E, Malats N. Development of a consensus molecular classifier for ',
              'pancreatic ductal adenocarcinoma. Genome Medicine 2025;17:142. ',
              tags$a(href = 'https://doi.org/10.1186/s13073-025-01568-9', target = '_blank',
                     '10.1186/s13073-025-01568-9')),
            p(strong('Cite the software: '),
              tags$a(href = 'https://doi.org/10.5281/zenodo.17019896', target = '_blank',
                     '10.5281/zenodo.17019896')),
            p(strong('Code: '),
              tags$a(href = 'https://github.com/pavillos/PDACMOC', target = '_blank',
                     'github.com/pavillos/PDACMOC')),
            p(strong('Contact: '),
              tags$a(href = 'mailto:pvilloslada@cnio.es', 'pvilloslada@cnio.es')))
      )),
      bslib::card(bslib::card_header('Use and citations'), bslib::card_body(
        div(class = 'help-block',
            metric.row('Citations',
                       metric.badge('Citations of the article', gist_json, '$.citations.article',
                                    'Genome Medicine', '6B2E99', 'https://europepmc.org/article/MED/41239365'),
                       metric.badge('Citations of the preprint', gist_json, '$.citations.preprint',
                                    'bioRxiv', 'B8925A', 'https://europepmc.org/article/PPR/PPR987525')),
            metric.row('Downloads',
                       tags$a(href = 'https://github.com/pavillos/PDACMOC/releases', target = '_blank',
                              tags$img(alt = 'Downloads from GitHub', src = paste0(
                                'https://img.shields.io/github/downloads/pavillos/PDACMOC/total',
                                '?label=GitHub&color=24292F'))),
                       metric.badge('Downloads from Zenodo', zenodo_json, '$.stats.downloads',
                                    'Zenodo', '185C84', 'https://doi.org/10.5281/zenodo.17019896')),
            metric.row('This app',
                       metric.badge('Classifications', usage_json, '$.classifications',
                                    'classifications since Sep 2026', '6B2E99', 'https://pdacmoc.cnio.es'),
                       metric.badge('Samples classified', usage_json, '$.samples',
                                    'samples classified since Sep 2026', '6B2E99', 'https://pdacmoc.cnio.es')),
            metric.row('License',
                       tags$a(href = 'https://github.com/pavillos/PDACMOC/blob/main/LICENSE', target = '_blank',
                              tags$img(alt = 'License: CC BY-NC 4.0',
                                       src = 'https://img.shields.io/badge/CC%20BY--NC%204.0-lightgrey'))),
            uiOutput('usage_line'))
      ))
    )
  ),

  bslib::nav_spacer(),
  bslib::nav_item(div(
    class = 'navbar-logos',
    div(class = 'logo-chip', imageOutput('cnio_logo', width = 'auto', height = '24px', inline = TRUE)),
    div(class = 'logo-chip', imageOutput('gmeg_logo', width = 'auto', height = '28px', inline = TRUE))
  )),
  bslib::nav_item(bslib::input_dark_mode(id = 'colour_mode', mode = 'light')),

  footer = div(
    class = 'app-footer',
    div(class = 'logo-row',
        imageOutput('pdacmoc_logo', width = 'auto', height = '54px', inline = TRUE),
        imageOutput('pdaconsensus_logo', width = 'auto', height = '26px', inline = TRUE))
  )
)
