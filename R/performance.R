#' @title Published balanced accuracy
#'
#' @author Villoslada-Blanco, Pablo
#'
#' @description
#' Mean balanced accuracy of every classifier as published (Genome Medicine 2025,
#' Tables 4 and 6), and the accuracy a classifier that always predicts the most
#' frequent class would reach, which depends on the number of subtypes.
#'
#' @return A named list with the tumor and stroma values

published.accuracy <- function() {
  list(
    tumor = data.frame(
      classifier = c('Collisson', 'Chan-Seng-Yue', 'PDAConsensus', 'Puleo', 'Bailey', 'Moffitt'),
      accuracy = c(97.44, 98.42, 96.33, 97.95, 94.95, 93.78),
      chance = c(33.33, 20, 50, 20, 25, 50),
      stringsAsFactors = FALSE),
    stroma = data.frame(
      classifier = c('Moffitt', 'PDAConsensus', 'Maurer'),
      accuracy = c(98.86, 98.92, 91.10),
      chance = c(50, 50, 50),
      stringsAsFactors = FALSE)
  )
}

#' @title Accuracy radar
#'
#' @author Villoslada-Blanco, Pablo
#'
#' @description
#' Radar chart of the balanced accuracy of a set of classifiers, with the accuracy
#' of a classifier that always predicts the most frequent class as a reference.
#'
#' @import ggplot2
#'
#' @param values Data frame with columns classifier, accuracy and chance
#' @param dark TRUE to draw the chart for the dark colour mode
#'
#' @return A ggplot object
#' @usage accuracy.radar(values, dark = FALSE)

accuracy.radar <- function(values, dark = FALSE) {

  ink <- if (dark) '#efe9f4' else '#1c1622'
  ink_soft <- if (dark) '#b2a7bd' else '#6a6175'
  line <- if (dark) '#5d5470' else '#ddd5e4'
  brand <- if (dark) '#c69bea' else '#6B2E99'

  n <- nrow(values)
  angle <- pi / 2 - 2 * pi * (seq_len(n) - 1) / n
  ring_at <- c(25, 50, 75, 100)

  # Closed polygons (the first vertex is repeated at the end)
  close_up <- function(r, a = angle) {
    data.frame(x = c(r * cos(a), r[1] * cos(a[1])), y = c(r * sin(a), r[1] * sin(a[1])))
  }
  rings <- do.call(rbind, lapply(ring_at, function(r) cbind(close_up(rep(r, n)), ring = r)))
  spokes <- data.frame(x = 100 * cos(angle), y = 100 * sin(angle))
  labels <- data.frame(
    x = 112 * cos(angle), y = 112 * sin(angle),
    classifier = values$classifier,
    value = sprintf('%.2f%%', values$accuracy),
    hjust = ifelse(abs(cos(angle)) < 0.3, 0.5, ifelse(cos(angle) > 0, 0, 1)),
    stringsAsFactors = FALSE)

  ggplot() +
    geom_path(data = rings, aes(.data$x, .data$y, group = .data$ring), colour = line, linewidth = 0.4) +
    geom_segment(data = spokes, aes(x = 0, y = 0, xend = .data$x, yend = .data$y),
                 colour = line, linewidth = 0.4) +
    geom_polygon(data = close_up(values$chance), aes(.data$x, .data$y),
                 fill = ink_soft, alpha = 0.12, colour = ink_soft, linewidth = 0.5, linetype = '22') +
    geom_polygon(data = close_up(values$accuracy), aes(.data$x, .data$y),
                 fill = brand, alpha = 0.2, colour = brand, linewidth = 0.9) +
    geom_point(data = data.frame(x = values$accuracy * cos(angle), y = values$accuracy * sin(angle)),
               aes(.data$x, .data$y), colour = brand, size = 2) +
    geom_text(data = labels, aes(.data$x, .data$y, label = .data$classifier, hjust = .data$hjust),
              vjust = 0, colour = ink_soft, size = 3.2, fontface = 'plain') +
    geom_text(data = labels, aes(.data$x, .data$y, label = .data$value, hjust = .data$hjust),
              vjust = 1.4, colour = ink, size = 3.2, fontface = 'bold') +
    # ring labels sit on the midpoint of the ring edge between the first two spokes
    geom_label(data = data.frame(x = ring_at * cos(pi / n) * cos(pi / 2 - pi / n),
                                 y = ring_at * cos(pi / n) * sin(pi / 2 - pi / n),
                                 label = paste0(ring_at, '%')),
               aes(.data$x, .data$y, label = .data$label), colour = ink_soft,
               fill = if (dark) '#2f2939' else '#ffffff', label.size = 0,
               label.padding = unit(0.08, 'lines'), size = 2.5) +
    # frame fitted to the shape, so a triangle is not pushed to the top
    coord_fixed(xlim = c(-138, 138), ylim = c(min(100 * sin(angle)) - 24, 118), clip = 'off') +
    theme_void() +
    theme(plot.background = element_rect(fill = NA, colour = NA),
          panel.background = element_rect(fill = NA, colour = NA),
          plot.margin = margin(4, 4, 4, 4))
}
