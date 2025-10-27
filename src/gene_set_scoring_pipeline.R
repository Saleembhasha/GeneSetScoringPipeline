
################################################################################
#             Gene set enrichment Analysis using different methodes            #
################################################################################
suppressPackageStartupMessages({
  # Required libraries
  library(Seurat)
  library(AUCell)
  library(UCell)
  library(GSVA)
  library(singscore)
  library(Matrix)
})


CalculateGeneSetScores <- function(seurat_obj = NULL, expr_matrix = NULL, gene_sets, methods_to_use = c("AddModuleScore", "AUCell", "UCell", "Custom_Average", "singscore", "GSVA")) {
  # Validate input
  if (is.null(seurat_obj) & is.null(expr_matrix)) {
    stop("Provide either a Seurat object or an expression matrix")
  }
  methods_all <- c("AddModuleScore", "AUCell", "UCell", "Custom_Average", "singscore", "GSVA")
  methods_to_use <- intersect(methods_to_use, methods_all)
  if (length(methods_to_use) == 0) stop("No valid methods specified in 'methods_to_use'")
  
  # Extract expression matrix
  if (!is.null(seurat_obj)) {
    assay_name <- DefaultAssay(seurat_obj)
    expr_matrix <- as.matrix(seurat_obj[[assay_name]]@data)
    if (is.null(expr_matrix) || ncol(expr_matrix) == 0) {
      expr_matrix <- as.matrix(seurat_obj[[assay_name]]@counts)
    }
  }
  
  # Filter gene sets for genes present in data
  gene_sets_filtered <- lapply(gene_sets, function(gs) intersect(gs, rownames(expr_matrix)))
  
  results <- list()
  
  # Run AddModuleScore if requested and Seurat object provided
  if ("AddModuleScore" %in% methods_to_use) {
    if (!is.null(seurat_obj)) {
      seurat_obj <- AddModuleScore(seurat_obj, features = gene_sets_filtered, name = "ModuleScore")
      results[["AddModuleScore"]] <- seurat_obj@meta.data[, grep("ModuleScore", colnames(seurat_obj@meta.data)), drop = FALSE]
    } else {
      warning("AddModuleScore requires a Seurat object; skipping.")
    }
  }
  
  # Run AUCell if requested
  if ("AUCell" %in% methods_to_use) {
    cells_rankings <- AUCell_buildRankings(expr_matrix, plotStats=FALSE, nCores=1)
    aucell_scores <- AUCell_calcAUC(gene_sets_filtered, cells_rankings, aucMaxRank = ceiling(0.05 * nrow(expr_matrix)))
    results[["AUCell"]] <- as.data.frame(t(aucell_scores@assays@data))
  }
  
  # Run UCell if requested
  if ("UCell" %in% methods_to_use) {
    results[["UCell"]] <- ScoreSignatures_UCell(expr_matrix, features = gene_sets_filtered)
  }
  
  # Run custom average score if requested
  if ("Custom_Average" %in% methods_to_use) {
    custom_score <- sapply(gene_sets_filtered, function(gs) {
      if (length(gs) == 0) {
        return(rep(NA, ncol(expr_matrix)))
      }
      colMeans(expr_matrix[gs, , drop=FALSE])
    })
    results[["Custom_Average"]] <- as.data.frame(t(custom_score))
  }
  
  # Run singscore if requested
  if ("singscore" %in% methods_to_use) {
    rank_mat <- rankGenes(expr_matrix)
    singscore_list <- lapply(gene_sets_filtered, function(gs) {
      if (length(gs) == 0) {
        return(rep(NA, ncol(expr_matrix)))
      }
      simpleScore(rank_mat, upSet = gs)$TotalScore
    })
    results[["singscore"]] <- as.data.frame(t(do.call(cbind, singscore_list)))
    rownames(results[["singscore"]]) <- names(gene_sets_filtered)
  }
  
  # Run GSVA if requested
  if ("GSVA" %in% methods_to_use) {
    gsva_res <- gsva(expr_matrix, gene_sets_filtered, method = "gsva", mx.diff = TRUE, verbose = FALSE)
    results[["GSVA"]] <- as.data.frame(t(gsva_res))
  }
  
  return(results)
}

# Example:
# gene_sets <- list(GeneSet1 = c("GeneA", "GeneB"), GeneSet2 = c("GeneC", "GeneD"))
# single_method_res <- CalculateGeneSetScores(seurat_obj, gene_sets, methods_to_use = "AUCell")
# multiple_methods_res <- CalculateGeneSetScores(expr_matrix, gene_sets, methods_to_use = c("GSVA", "singscore"))
