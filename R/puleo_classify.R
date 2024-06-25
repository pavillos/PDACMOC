#' @title Puleo classifier
#' 
#' @author Villoslada-Blanco, Pablo
#' 
#' @description
#' Classification of tumor according to Puleo et al., 2018
#' 
#' @import reticulate
#'
#' @param new_samples Final expression matrix after batch correction, filtering, 
#' normalization, scaling, and gene ID change
#'
#' @return results_puleo: classification according to Puleo
#' 
#' @example ./inst/examples/example.R

puleo.classify <- function(new_samples) {

  # Get Puleo genes and extract common genes
  file1 <- system.file('training_data', 'all_Puleo.csv', package = 'PDACMOC')
  model_data <- read.csv(file1, header = TRUE, row.names = 1, check.names = FALSE)
  subtype <- model_data[1]
  model_data <- model_data[,-1]
  rm(file1)
  genes <- colnames(model_data)
  common_genes <- intersect(genes, rownames(new_samples))
  genes_imputed <- length(genes) - length(common_genes)
  message('Puleo: the expression of ', genes_imputed, ' genes will be imputed')
  new_samples_filtered <- new_samples[common_genes, ]
  rm(genes, common_genes, genes_imputed)
  
  # Impute the expression of missing genes
  new_samples_imputed <- matrix(nrow = ncol(model_data),
                                ncol = ncol(new_samples_filtered))
  n <- 1
  for (g in colnames(model_data)) {
    if (g %in% rownames(new_samples_filtered)) {
      new_samples_imputed[n,] <- as.numeric(new_samples_filtered[g,])
    }
    else {
      new_samples_imputed[n,] <- rep(mean(aggregate(model_data[,g], subtype, FUN = mean)$x),
                                     ncol(new_samples_filtered))
    }
    n = n + 1
  }
  rownames(new_samples_imputed) <- colnames(model_data)
  colnames(new_samples_imputed) <- colnames(new_samples_filtered)
  rm(g, n, new_samples_filtered)
  data2ML <- t(new_samples_imputed)
  rm(model_data, new_samples_imputed)
  
  # Get labels and class probabilites
  file2 <- system.file('models', 'Puleo.pkl', package = 'PDACMOC')
  trained_model <- py_load_object(file2)
  rm(file2)
  labels <- trained_model$model$predict(data2ML)
  probabilities <- trained_model$model$predict_proba(data2ML)
  rm(trained_model)
  probability <- c()
  for (sample in 1:nrow(probabilities)) {
    max_value <- max(probabilities[sample,])
    max_value <- round(max_value * 100)
    probability <- c(probability, max_value)
  }
  rm(sample, max_value, probabilities)
  results_puleo <- data.frame('Predicted.subtype' = factor(labels, levels = c('Pure_classical', 'Immune_classical', 'Pure_basal-like', 'Desmoplastic', 'Stroma_activated')),
                              'Probability' = probability,
                              row.names = rownames(data2ML))
  rm(data2ML, labels, probability)
  
  return(results_puleo)
  
}
