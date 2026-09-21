#' @title Subtype colours
#'
#' @author Villoslada-Blanco, Pablo
#'
#' @description
#' Colours of every subtype, taken from the figures of the paper
#' (annotation_colors_tumor and annotation_colors_stroma). They are set per classifier
#' because some subtype names (e.g. Classical) exist in more than one classifier.
#'
#' @param classifier Classifier name as used in the results (e.g. 'Moffitt',
#' 'Stroma Moffitt', 'PDAConsensus', 'Stroma PDAConsensus')
#'
#' @return Named character vector of hex colours

subtype.colours <- function(classifier) {
  palettes <- list(
    'Collisson' = c('Classical' = '#7570B3', 'QM' = '#E7298A', 'Exocrine-like' = '#66A61E'),
    'Moffitt' = c('Classical' = '#E6AB02', 'Basal-like' = '#A6761D'),
    'Bailey' = c('Squamous' = '#96969c', 'Progenitor' = '#1B9E77', 'Immunogenic' = '#D95F02',
                 'ADEX' = '#7570B3'),
    'Puleo' = c('Pure_classical' = '#E7298A', 'Immune_classical' = '#A6761D',
                'Pure_basal-like' = '#66A61E', 'Desmoplastic' = '#E6AB02',
                'Stroma_activated' = '#c7c7cd'),
    'Chan-Seng-Yue' = c('Classical-A' = '#A6CEE3', 'Classical-B' = '#1F78B4',
                        'Basal-like-A' = '#B2DF8A', 'Basal-like-B' = '#33A02C', 'Hybrid' = '#FB9A99'),
    'PDAConsensus' = c('Consensus_Classical' = '#E6ED17', 'Consensus_Non-classical' = '#626266'),
    'Stroma Moffitt' = c('Normal' = '#7EB5D6', 'Activated' = '#C5E1A5'),
    'Stroma Maurer' = c('ECM-rich' = '#31708E', 'Immune-rich' = '#9C27B0'),
    'Stroma PDAConsensus' = c('Consensus_Normal-immune' = '#0FF7DB',
                              'Consensus_Activated-ECM' = '#765B0D')
  )
  palettes[[classifier]]
}

# Readable label for a subtype level: underscores become spaces
subtype.label <- function(x) gsub('_', ' ', x)

# Dark or light text for a filled background, by relative luminance
text.on <- function(hex) {
  rgb <- grDevices::col2rgb(hex) / 255
  lum <- sum(c(0.2126, 0.7152, 0.0722) * ifelse(rgb <= 0.03928, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4))
  if (lum > 0.4) '#1c1622' else '#ffffff'
}

#' @title Classification table
#'
#' @author Villoslada-Blanco, Pablo
#'
#' @description
#' Table shown in the Results tab for one classifier: the subtype as a coloured label,
#' the probability, a separate low-confidence flag, and the continuous score when the
#' classifier has one. Only the display changes; downloads use the original results.
#'
#' @import DT
#'
#' @param results Data frame returned by a classifier (Predicted.subtype, Probability
#' and, for PDAConsensus, a score column)
#' @param classifier Classifier name, see \code{subtype.colours}
#' @param threshold Probability (in \%) below which a sample is low confidence
#'
#' @return A DT datatable

classification.table <- function(results, classifier, threshold) {

  colours <- subtype.colours(classifier)
  subtype <- as.character(results$Predicted.subtype)
  swatch <- ifelse(subtype %in% names(colours), colours[subtype], '#999999')

  display <- data.frame(
    Sample = htmltools::htmlEscape(rownames(results)),
    Subtype = sprintf('<span class="subtype" style="--c:%s"><i></i>%s</span>',
                      swatch, htmltools::htmlEscape(subtype.label(subtype))),
    Probability = results$Probability,
    Confidence = ifelse(results$Probability < threshold,
                        '<span class="confidence low">Low</span>',
                        '<span class="confidence high">High</span>'),
    check.names = FALSE, stringsAsFactors = FALSE)

  score_col <- intersect(c('NonClassicalScore', 'ActivatedECMScore'), colnames(results))
  if (length(score_col) == 1) {
    display[[score_col]] <- round(results[[score_col]], 4)
  }

  table <- datatable(
    display, rownames = FALSE, escape = FALSE, class = 'compact hover',
    options = list(pageLength = 25, dom = if (nrow(display) > 25) 'tip' else 't',
                   paging = nrow(display) > 25, autoWidth = TRUE,
                   columnDefs = list(list(className = 'dt-right', targets = 2),
                                     list(orderable = FALSE, targets = 3)))) %>%
    formatString('Probability', suffix = '%')

  if (length(score_col) == 1) {
    bar <- if (score_col == 'NonClassicalScore') colours[['Consensus_Non-classical']] else colours[['Consensus_Activated-ECM']]
    table <- table %>%
      formatRound(score_col, 4) %>%
      formatStyle(score_col,
                  background = styleColorBar(c(0, 1), bar),
                  backgroundSize = '100% 4px', backgroundRepeat = 'no-repeat',
                  backgroundPosition = 'right bottom 4px')
  }
  table
}

#' @title Subtype summary bar
#'
#' @author Villoslada-Blanco, Pablo
#'
#' @description
#' Stacked bar with the number of samples assigned to each subtype, in the colours
#' of the paper, followed by a legend with the counts.
#'
#' @param results Data frame returned by a classifier
#' @param classifier Classifier name, see \code{subtype.colours}
#' @param threshold Probability (in \%) below which a sample is low confidence
#'
#' @return An HTML tag

subtype.summary <- function(results, classifier, threshold) {
  colours <- subtype.colours(classifier)
  counts <- table(factor(as.character(results$Predicted.subtype), levels = names(colours)))
  n <- sum(counts)
  low <- sum(results$Probability < threshold)
  shown <- counts[counts > 0]

  segments <- lapply(names(shown), function(s) {
    div(class = 'seg',
        style = sprintf('width:%.4f%%;background:%s;color:%s', 100 * shown[[s]] / n, colours[[s]],
                        text.on(colours[[s]])),
        title = sprintf('%s: %d', subtype.label(s), shown[[s]]),
        if (shown[[s]] / n >= 0.12) shown[[s]])
  })
  legend <- lapply(names(counts), function(s) {
    span(class = 'legend-item',
         tags$i(style = sprintf('background:%s', colours[[s]])),
         sprintf('%s %d', subtype.label(s), counts[[s]]))
  })

  div(class = 'subtype-summary',
      div(class = 'summary-meta',
          span(sprintf('%d samples', n)),
          if (low > 0) span(class = 'low-confidence', sprintf('%d low confidence', low))),
      div(class = 'segbar', segments),
      div(class = 'summary-legend', legend))
}
