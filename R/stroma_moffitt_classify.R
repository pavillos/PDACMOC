#' @title Moffitt classifier of stroma
#' 
#' @author Villoslada-Blanco, Pablo
#' 
#' @description
#' Classification of stroma according to Moffitt et al., 2015
#' 
#' @import reticulate
#'
#' @param vm_S Stroma fraction after virtual microdissection. It is produced by
#' virtual.microdissect function
#'
#' @return results_stroma_moffitt: classification according to Moffitt of stroma
#' samples
#' 
#' @example ./inst/examples/example.R

stroma.moffitt.classify <- function(vm_S) {
  
  # Get Moffitt stromal genes and extract common genes
  file1 <- system.file('training_data', 'all_Moffitt_stroma.csv', package = 'PDACMOC')
  model_data <- read.csv(file1, header = TRUE, row.names = 1, check.names = FALSE)
  subtype <- model_data[1]
  model_data <- model_data[,-1]
  rm(file1)
  genes <- colnames(model_data)
  common_genes <- intersect(genes, rownames(vm_S))
  genes_imputed <- length(genes) - length(common_genes)
  message('Stroma Moffitt: the expression of ', genes_imputed, ' genes will be imputed')
  vm_S_filtered <- vm_S[common_genes, ]
  rm(genes, common_genes, genes_imputed)
  
  # Impute the expression of missing genes
  vm_S_imputed <- matrix(nrow = ncol(model_data), ncol = ncol(vm_S_filtered))
  n <- 1
  for (g in colnames(model_data)) {
    if (g %in% rownames(vm_S_filtered)) {
      vm_S_imputed[n,] <- as.numeric(vm_S_filtered[g,])
    }
    else {
      vm_S_imputed[n,] <- rep(mean(aggregate(model_data[,g], subtype, FUN = mean)$x),
                              ncol(vm_S_filtered))
    }
    n = n + 1
  }
  rownames(vm_S_imputed) <- colnames(model_data)
  colnames(vm_S_imputed) <- colnames(vm_S_filtered)
  rm(g, n, subtype, vm_S_filtered)
  vmdata2ML <- t(vm_S_imputed)
  rm(model_data, vm_S_imputed)
  
  # Get labels and class probabilites of virtual microdissected samples
  file2 <- system.file('models', 'Moffitt_stroma.pkl', package = 'PDACMOC')
  trained_model <- py_load_object(file2)
  rm(file2)
  labels <- trained_model$model$predict(vmdata2ML)
  probabilities <- trained_model$model$predict_proba(vmdata2ML)
  rm(trained_model)
  probability <- c()
  for (sample in 1:nrow(probabilities)) {
    max_value <- max(probabilities[sample,])
    max_value <- round(max_value * 100)
    probability <- c(probability, max_value)
  }
  rm(sample, max_value, probabilities)
  results_stroma_moffitt <- data.frame('Predicted.subtype' = factor(labels, levels = c('Normal', 'Activated')),
                                       'Probability' = probability,
                                       row.names = rownames(vmdata2ML))
  rm(vmdata2ML, labels, probability)
  
  return(results_stroma_moffitt)
  
}
