#' @title  Import and normalize data
#' 
#' @author Villoslada-Blanco, Pablo
#' 
#' @description
#' Import matrix of RNAseq raw counts, batch correction, filtering, and normalization using VST
#' (Variance Stabilizing Transformation) and standard scaler.
#' The matrix must have genes in rows and samples in columns
#' 
#' @import DESeq2
#' @import org.Hs.eg.db
#' @import AnnotationDbi
#' @import dplyr
#' @import sva
#'
#' @param expmat Matrix of RNAseq raw counts
#' @param batch Boolean to apply or not batch correction (TRUE or FALSE; default = TRUE)
#' @param gene_id Type of gene ID used (EnsemblID, EntrezID, or GeneSymbol; default = EnsemblID)
#'
#' @return new_samples final expression matrix after batch correction, filtering, 
#' normalization, scaling, and gene ID change
#' 
#' @example ./inst/examples/example.R

import.and.normalize <- function(expmat, batch = TRUE, gene_id = 'EnsemblID') {
  
  # Check gene ID
  if (!gene_id %in% c('EnsemblID', 'EntrezID', 'GeneSymbol')) {
    stop('The gene IDs are not correct')
  }
  
  # Convert gene ID to EnsemblID
  if (gene_id == 'EntrezID') {
    gene_ids <- as.character(suppressMessages(mapIds(org.Hs.eg.db,
                                                     rownames(expmat),
                                                     keytype = 'ENTREZ',
                                                     column = 'ENSEMBL')))
    expmat <- expmat[!is.na(gene_ids), ]
    gene_ids <- gene_ids[!is.na(gene_ids)]
    rownames(expmat) <- gene_ids
    rm(gene_ids)
  } else if (gene_id == 'GeneSymbol') {
    gene_ids <- as.character(suppressMessages(mapIds(org.Hs.eg.db,
                                                     rownames(expmat),
                                                     keytype = 'SYMBOL',
                                                     column = 'ENSEMBL')))
    expmat <- expmat[!is.na(gene_ids), ]
    gene_ids <- gene_ids[!is.na(gene_ids)]
    rownames(expmat) <- gene_ids
    rm(gene_ids)
  } else if (gene_id == 'EnsemblID') {
    # Discard gene version
    if (any(grepl('\\.', rownames(expmat)))) {
      message('The gene IDs include version and is going to be discarded')
      gene_ids <- gsub('\\..*', '', rownames(expmat))
      expmat <- as.data.frame(expmat)
      expmat$GeneID <- gene_ids
      expmat_mean <- expmat %>% group_by(GeneID) %>% summarise_all(mean, na.rm = TRUE)
      expmat <- as.data.frame(expmat_mean[, -1])
      rownames(expmat) <- expmat_mean$GeneID
      rm(gene_ids, expmat_mean)
    }
  }
  
  # Apply batch correction (if there are more than 1 sample)
  if (ncol(expmat) > 1 & batch) {
    all_datasets_file <- system.file('training_data', 'all_datasets.csv', package = 'PDACMOC')
    all_datasets <- read.csv(all_datasets_file, row.names = 1, check.names = FALSE)
    common_genes <- intersect(rownames(expmat), rownames(all_datasets))
    expmat <- expmat[common_genes, ]
    all_datasets <- all_datasets[common_genes, ]
    all_datasets <- cbind(expmat, all_datasets)
    batch_info <- c(rep('new_samples', ncol(expmat)),
                    rep('TCGA', 149),
                    rep('ICGC_v100', 68),
                    rep('ICGC_v84', 190),
                    rep('PanGenEU', 107))
    batch_info <- factor(batch_info, levels = c('new_samples', 'TCGA', 'ICGC_v100', 'ICGC_v84', 'PanGenEU'))
    if (shiny::isRunning()) {
      withProgress(message = 'Applying batch correction', value = 0.2, {
        tmp <- capture.output(all_datasets_corrected <- ComBat_seq(counts = as.matrix(all_datasets), batch = batch_info))
      })
    } else {tmp <- capture.output(all_datasets_corrected <- ComBat_seq(counts = as.matrix(all_datasets), batch = batch_info))}
    rm(tmp)
    expmat <- all_datasets_corrected[, 1:ncol(expmat)]
    rm(all_datasets_file, all_datasets, common_genes, all_datasets_corrected, batch_info)
  }
  
  # Filter out genes with fewer than 5 reads in 50% of samples
  # We do not discard genes that are include in our list of genes used to train
  # ML model, because they would be imputed later and this does not make sense
  num_samples <- ncol(expmat) * 0.5
  file <- system.file('gene_signatures', 'all_genes.csv', package = 'PDACMOC')
  genes <- read.csv(file, header = FALSE)
  genes <- genes$V1
  genes_mod <- as.character(suppressMessages(mapIds(org.Hs.eg.db,
                                                    genes,
                                                    keytype = 'SYMBOL',
                                                    column = 'ENSEMBL')))
  expmat <- expmat[rowSums(expmat >= 5) >= num_samples | rownames(expmat) %in% genes_mod, ]
  rm(genes, genes_mod)
  
  # Apply VST (if there are more than 1 sample) and standard scaler
  expmat <- as.matrix(expmat)
  if (ncol(expmat) > 1) {
    suppressMessages({expmat_vst <- varianceStabilizingTransformation(expmat, blind = FALSE)})
  } else {
    expmat_vst <- expmat
  }
  expmat_scaled <- scale(expmat_vst)
  if (all(dim(expmat) == dim(expmat_vst)) | all(dim(expmat) == dim(expmat_scaled))) {
  } else {
    stop('The dimensions of the normalized expression matrix are not correct')
  }
  rm(expmat, expmat_vst)
  
  # Convert EnsemblID to GeneSymbol
  new_samples <- expmat_scaled
  rownames(new_samples) <- as.character(suppressMessages(mapIds(org.Hs.eg.db,
                                                                rownames(expmat_scaled),
                                                                keytype = 'ENSEMBL',
                                                                column = 'SYMBOL')))
  new_samples <- new_samples[!is.na(rownames(new_samples)), ]
  genes_lost <- nrow(expmat_scaled) - nrow(new_samples) 
  message(genes_lost, ' genes were not been able to convert from ', gene_id,
          ' to GeneSymbol')
  stopifnot(ncol(new_samples) == ncol(expmat_scaled),
            nrow(new_samples) == (nrow(expmat_scaled) - genes_lost))
  rm(expmat_scaled)

  return(new_samples)
    
}
