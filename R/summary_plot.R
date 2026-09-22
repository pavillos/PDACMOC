#' @title Summary figure
#'
#' @author Villoslada-Blanco, Pablo
#'
#' @description
#' The Summary tab of the Shiny app as a figure for download: one panel per classifier with
#' the number of samples assigned to each subtype, in the colours of the paper, and the
#' number of samples below the confidence threshold of the classifier.
#'
#' @param result Result of \code{omni.classify}
#' @param settings List with the classifiers that were run (tumor_classifiers, stroma,
#' stroma_classifiers)
#'
#' @return A list with the ggplot object and its size in inches (width, height)

summary.plot <- function(result, settings) {
  # classifier, element of the result, palette name, threshold (as in the Summary tab)
  tumor <- list(c('Collisson', 'Collisson', 'Collisson', 50), c('Moffitt', 'Moffitt', 'Moffitt', 70),
                c('Bailey', 'Bailey', 'Bailey', 33), c('Puleo', 'Puleo', 'Puleo', 30),
                c('Chan-Seng-Yue', 'Chan-Seng-Yue', 'Chan-Seng-Yue', 30),
                c('PDAConsensus', 'PDAConsensus', 'PDAConsensus', 70))
  stroma <- list(c('Moffitt', 'Moffitt', 'Stroma Moffitt', 70), c('Maurer', 'Maurer', 'Stroma Maurer', 70),
                 c('PDAConsensus', 'PDAConsensus', 'Stroma PDAConsensus', 70))

  panel <- function(section, spec, results) {
    colours <- subtype.colours(spec[3])
    counts <- table(factor(as.character(results$Predicted.subtype), levels = names(colours)))
    low <- sum(results$Probability < as.numeric(spec[4]))
    title <- sprintf('%s: %s\n%d samples%s', section, spec[1], sum(counts),
                     if (low > 0) sprintf(', %d low confidence', low) else '')
    data.frame(panel = title, subtype = subtype.label(names(counts)), n = as.integer(counts),
               key = paste(spec[3], names(counts)), colour = unname(colours), stringsAsFactors = FALSE)
  }

  parts <- list()
  for (spec in tumor) {
    if (spec[1] %in% settings$tumor_classifiers) parts[[length(parts) + 1]] <- panel('Tumor', spec, result$Tumor[[spec[2]]])
  }
  if (isTRUE(settings$stroma)) {
    for (spec in stroma) {
      if (spec[1] %in% settings$stroma_classifiers) parts[[length(parts) + 1]] <- panel('Stroma', spec, result$Stroma[[spec[2]]])
    }
  }
  data <- do.call(rbind, parts)
  data$panel <- factor(data$panel, levels = unique(data$panel))
  # subtypes from top to bottom in the order of the palette
  data$row <- factor(paste(data$panel, data$subtype), levels = rev(unique(paste(data$panel, data$subtype))))
  fills <- stats::setNames(data$colour, data$key)

  plot <- ggplot(data, aes(x = n, y = row, fill = key)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n), hjust = -0.25, size = 3.4, colour = '#2b2233') +
    facet_wrap(~ panel, ncol = 3, scales = 'free_y') +
    scale_y_discrete(labels = stats::setNames(data$subtype, as.character(data$row))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    scale_fill_manual(values = fills, guide = 'none') +
    labs(x = 'Samples', y = NULL,
         title = 'PDACMOC classification summary',
         caption = sprintf('PDACMOC %s, pdacmoc.cnio.es', utils::packageVersion('PDACMOC'))) +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(hjust = 0, face = 'bold', size = 10, lineheight = 1.1),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(face = 'bold', colour = '#6B2E99'),
          plot.caption = element_text(colour = 'grey45', size = 8),
          plot.background = element_rect(fill = 'white', colour = NA))

  rows <- ceiling(nlevels(data$panel) / 3)
  list(plot = plot, width = 11, height = 1.2 + 2.3 * rows)
}
