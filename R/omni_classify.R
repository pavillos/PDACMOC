#' @title PDACMMolecularOmniClassifier
#' 
#' @author Villoslada-Blanco, Pablo
#' 
#' @description
#' Classification of tumor and stroma according to different classifiers.
#' Classifiers available for tumor fractions:
#'  1. Collisson et al., 2011
#'  2. Moffitt et al., 2015
#'  3. Bailey et al., 2016
#'  4. Puleo et al, 2018
#'  5. Chan-Seng-Yue, 2020
#'  6. PDAConsensus
#' Classifiers available for stroma fractions:
#'  1. Moffitt et al., 2016
#'  2. Maurer et al., 2019
#'  3. PDAConsensus
#'
#' @import shiny
#'
#' @param expmat Matrix of RNAseq raw counts
#' @param batch Boolean to apply or not batch correction (TRUE or FALSE; default = TRUE)
#' @param gene_id Type of gene ID used (EnsemblID, EntrezID, or GeneSymbol; default = EnsemblID)
#' @param classifier Classifier desired to make the classification of tumor fractions
#' ('Collisson', 'Moffitt', 'Bailey', 'Puleo', 'Chan-Seng-Yue', and 'PDAConsensus'; 
#' default = c('Collisson', 'Moffitt', 'Bailey', 'Puleo', 'Chan-Seng-Yue', 'PDAConsensus'))
#' @param stroma Boolean to include or not classification of stroma (TRUE or FALSE; default = FALSE)
#' @param stroma_classifier Classifier desired to make the classification of
#' stroma fractions ('Moffit', 'Maurer', and 'PDAConsensus'; default = '')
#'
#' @return classification: a list containing three list. The first list contains 
#' the proportions of tumor and stroma of the samples. The second list corresponds
#' to the classification of tumor and it contains as much lists as classifiers length.
#' The third list corresponds to the classification of stroma and it contains
#' as much lists as stroma_classifiers length. Second and third list are only
#' included if stroma = TRUE  
#' 
#' @example ./inst/examples/example.R
#' 
#' @export

omni.classify <- function(expmat, batch = TRUE, gene_id = 'EnsemblID',
                          classifier = c('Collisson', 'Moffitt', 'Bailey',
                                         'Puleo', 'Chan-Seng-Yue', 'PDAConsensus'),
                          stroma = FALSE, stroma_classifier = '') {
  if (shiny::isRunning()) {withProgress(message = 'Classification started', value = 0, {Sys.sleep(3)})}
  
  if (any(!classifier %in% c('Collisson', 'Moffitt', 'Bailey',
                             'Puleo', 'Chan-Seng-Yue', 'PDAConsensus'))) {
    stop('The classifier is not correct. Available options are Collisson, Moffitt, Bailey, Puleo, Chan-Seng-Yue, and PDAConsensus')
  }
  
  if(length(stroma) > 1) {
    stop('Stroma parameter only accepts a unique value')
  }
  
  if (!stroma %in% c('TRUE', 'FALSE')) {
    stop('Stroma parameter must be TRUE or FALSE')
  }
  
  if(stroma & any(!stroma_classifier %in% c('Moffitt', 'Maurer', 'PDAConsensus'))) {
    stop('The stroma classifier is not correct. Available options are Moffitt, Maurer, and PDAConsensus')
  }
  
  if(stroma == FALSE & any(stroma_classifier %in% c('Moffitt', 'Maurer', 'PDAConsensus'))) {
    stroma <- TRUE
    message('You have not selected stroma classification but stroma classifiers have been provided. Thus, stroma is forced to TRUE')
  }
  
  # Import data
  if (shiny::isRunning()) {
    withProgress(message = 'Importing samples', value = 0.1, {
      new_samples <- import.and.normalize(expmat, batch, gene_id)
    })
  } else {new_samples <- import.and.normalize(expmat, batch, gene_id)}
  
  # Get virtual microdissection if stroma classification is selected
  if (stroma == TRUE) {
    if (shiny::isRunning()) {
      withProgress(message = 'Virtual microdissection. This may take a while...', value = 0.4, {
        vm_result <- virtual.microdissect(new_samples)
      })
    } else {
      vm_result <- virtual.microdissect(new_samples)
    }
  }
  
  # Classify tumor samples
  if ('Collisson' %in% classifier) {
    if (shiny::isRunning()) {withProgress(message = 'Collisson classification', value = 0.9, {Sys.sleep(3)})}
    results_collisson <- collisson.classify(new_samples)
  }
  
  if ('Moffitt' %in% classifier) {
    if (shiny::isRunning()) {withProgress(message = 'Moffitt classification', value = 0.9, {Sys.sleep(3)})}
    results_moffitt <- moffitt.classify(new_samples)
  }
  
  if ('Bailey' %in% classifier) {
    if (shiny::isRunning()) {withProgress(message = 'Bailey classification', value = 0.9, {Sys.sleep(3)})}
    results_bailey <- bailey.classify(new_samples)
  }
  
  if ('Puleo' %in% classifier) {
    if (shiny::isRunning()) {withProgress(message = 'Puleo classification', value = 0.9, {Sys.sleep(3)})}
    results_puleo <- puleo.classify(new_samples)
  }
  
  if ('Chan-Seng-Yue' %in% classifier) {
    if (shiny::isRunning()) {withProgress(message = 'Chan-Seng-Yue classification', value = 0.9, {Sys.sleep(3)})}
    results_chan <- chan.classify(new_samples)
  }
  
  if ('PDAConsensus' %in% classifier) {
    if (shiny::isRunning()) {withProgress(message = 'PDAConsensus classification', value = 0.9, {Sys.sleep(3)})}
    results_consensus <- PDAConsensus.classify(new_samples)
  }
  rm(new_samples)
  
  results_tumor <- list()
  if ('Collisson' %in% classifier) {
    results_tumor$Collisson <- results_collisson
    rm(results_collisson)
  }
  if ('Moffitt' %in% classifier) {
    results_tumor$Moffitt <- results_moffitt
    rm(results_moffitt)
  }
  if ('Bailey' %in% classifier) {
    results_tumor$Bailey <- results_bailey
    rm(results_bailey)
  }
  if ('Puleo' %in% classifier) {
    results_tumor$Puleo <- results_puleo
    rm(results_puleo)
  }
  if ('Chan-Seng-Yue' %in% classifier) {
    results_tumor$`Chan-Seng-Yue` <- results_chan
    rm(results_chan)
  }
  if ('PDAConsensus' %in% classifier) {
    results_tumor$PDAConsensus <- results_consensus
    rm(results_consensus)
  }
  
  # Classify stromal samples
  if (stroma) {
    classification <- list('Proportions' = vm_result$proportions, 'Tumor' = results_tumor) # we include proportions of tumor and stroma in the output
    rm(results_tumor)
    
    if ('Moffitt' %in% stroma_classifier) {
      if (shiny::isRunning()) {withProgress(message = 'Moffitt stroma classification', value = 0.95, {Sys.sleep(3)})}
      results_stroma_moffitt <- stroma.moffitt.classify(vm_result$vm_S)
    }
    if ('Maurer' %in% stroma_classifier) {
      if (shiny::isRunning()) {withProgress(message = 'Maurer stroma classification', value = 0.95, {Sys.sleep(3)})}
      results_stroma_maurer <- stroma.maurer.classify(vm_result$vm_S)
    }
    if ('PDAConsensus' %in% stroma_classifier) {
      if (shiny::isRunning()) {withProgress(message = 'PDAConsensus stroma classification', value = 0.95, {Sys.sleep(3)})}
      results_stroma_consensus <- stroma.PDAConsensus.classify(vm_result$vm_S)
    }

    results_stroma <- list()
    if ('Moffitt' %in% stroma_classifier) {
      results_stroma$Moffitt <- results_stroma_moffitt
      rm(results_stroma_moffitt)
    }
    if ('Maurer' %in% stroma_classifier) {
      results_stroma$Maurer <- results_stroma_maurer
      rm(results_stroma_maurer)
    }
    if ('PDAConsensus' %in% stroma_classifier) {
      results_stroma$PDAConsensus <- results_stroma_consensus
      rm(results_stroma_consensus)
    }
    classification$Stroma <- results_stroma
    rm(results_stroma, vm_result)
  } else {
    classification <- list('Tumor' = results_tumor)
    rm(results_tumor)
  }
  
  if (shiny::isRunning()) {withProgress(message = 'Classification finished', value =1, {Sys.sleep(3)})}
  
  return(classification)
  
}
