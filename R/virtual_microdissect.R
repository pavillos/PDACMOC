#' @title Virtual microdissection
#' 
#' @author Villoslada-Blanco, Pablo
#' 
#' @description
#' Apply virtual microdissection using ADVOCATE package
#' 
#' @import ADVOCATE
#'
#' @param new_samples Final expression matrix after batch correction, filtering, 
#' normalization, scaling, and gene ID change
#'
#' @return vm_results: a list with three matrixes:
#'  1. new_samples: Final expression matrix after batch correction, filtering, 
#'     normalization, scaling, and gene id change
#'  2. vm_E: tumor fraction after virtual microdissection
#'  3. vm_S: stroma fraction after virtual microdissection
#'  4. proportions: proportions of tumor and stroma
#' 
#' @example ./inst/examples/example.R

virtual.microdissect <- function(new_samples) {
  
  pdf('./tmp.pdf')
  
  # Get train data
  train <- system.file('training_data', 'trainDataADVOCATE.rda', package = 'PDACMOC')
  
  # Compute proportions of epithelium, stroma, and others
  suppressPackageStartupMessages({
    tmp <- capture.output(prop <- predict_bulk_3comp(train, new_samples, epsilon = 0.001))
  })
  rm(tmp)
  
  # Get virtual microdissection
  load(train)
  tmp <- capture.output(vm <- calCellTypeExpression_3comp(expmat, deg, fc, pval,
                                                          sampleInfo,
                                                          new_samples,
                                                          prop, method = 'lcm'))
  vm_E <- subset(vm, select = endsWith(names(vm), '_E'))
  names(vm_E) <- sub('_E$', '', names(vm_E))
  vm_S <- subset(vm, select = endsWith(names(vm), '_S'))
  names(vm_S) <- sub('_S$', '', names(vm_S))
  rm(train, vm)
  
  prop <- round(prop, 3)
  vm_result <- list('new_samples' = new_samples, 'vm_E' = vm_E, 'vm_S' = vm_S, 'proportions' = prop)
  rm(new_samples, vm_E, vm_S, prop)
  
  dev.off()
  file.remove('./tmp.pdf')
  
  return(vm_result)
  
}
