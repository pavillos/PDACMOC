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
  
  # ADVOCATE draws plots: send them to a temporary pdf that is removed at the end
  tmp_pdf <- tempfile(fileext = '.pdf')
  pdf(tmp_pdf)
  pdf_device <- dev.cur()
  on.exit({
    if (pdf_device %in% dev.list()) dev.off(pdf_device)
    unlink(tmp_pdf)
  }, add = TRUE)
  
  # ADVOCATE fails with a single sample, so it is duplicated and only the first copy is kept
  # (ADVOCATE processes each sample independently)
  single_sample <- ncol(new_samples) == 1
  if (single_sample) {
    new_samples <- cbind(new_samples, new_samples)
    colnames(new_samples)[2] <- paste0(colnames(new_samples)[1], '_copy')
  }
  
  # Number of parallel workers used by ADVOCATE (option PDACMOC.cores, default 4)
  old_options <- options(ADVOCATE.cores = getOption('PDACMOC.cores', 4))
  on.exit(options(old_options), add = TRUE)
  
  # Get train data
  train <- system.file('training_data', 'trainDataADVOCATE.rda', package = 'PDACMOC')
  load(train)
  
  run_advocate <- function() {
    # Compute proportions of epithelium, stroma, and others
    suppressPackageStartupMessages({
      tmp <- capture.output(prop <- predict_bulk_3comp(train, new_samples, epsilon = 0.001))
    })
    # A failed parallel task does not stop ADVOCATE, so check the proportions
    if (NROW(prop) != ncol(new_samples) || anyNA(prop)) {
      stop('the proportions could not be estimated for all samples')
    }
    # Get virtual microdissection
    tmp <- capture.output(vm <- calCellTypeExpression_3comp(expmat, deg, fc, pval,
                                                            sampleInfo,
                                                            new_samples,
                                                            prop, method = 'lcm'))
    list(prop = prop, vm = vm)
  }
  
  # Occasional failures of the parallel workers: retry once
  advocate <- tryCatch(run_advocate(), error = function(e) {
    message('Virtual microdissection failed (', conditionMessage(e), '). Retrying once')
    run_advocate()
  })
  prop <- advocate$prop
  vm <- advocate$vm
  rm(advocate)
  
  vm_E <- subset(vm, select = endsWith(names(vm), '_E'))
  names(vm_E) <- sub('_E$', '', names(vm_E))
  vm_S <- subset(vm, select = endsWith(names(vm), '_S'))
  names(vm_S) <- sub('_S$', '', names(vm_S))
  rm(train, vm, expmat, deg, fc, pval, sampleInfo)
  
  if (single_sample) {
    new_samples <- new_samples[, 1, drop = FALSE]
    vm_E <- vm_E[, 1, drop = FALSE]
    vm_S <- vm_S[, 1, drop = FALSE]
    prop <- prop[1, , drop = FALSE]
  }
  
  prop <- round(prop, 3)
  vm_result <- list('new_samples' = new_samples, 'vm_E' = vm_E, 'vm_S' = vm_S, 'proportions' = prop)
  rm(new_samples, vm_E, vm_S, prop)
  
  return(vm_result)
  
}
