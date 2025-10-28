
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
  if (!is.null(seurat_obj)){
    avail <- rownames(GetAssayData(seurat_obj, assay = assay_name, layer = layer))
  }else if (!is.null(expr_matrix)){
    avail <- intersect(gs, rownames(expr_matrix))
  }
  
  gene_sets_filtered <- lapply(gene_sets, function(gs) intersect(unique(gs), avail))
  
  results <- list()
  
  # Run AddModuleScore if requested and Seurat object provided
  if ("AddModuleScore" %in% methods_to_use) {
    if (!is.null(seurat_obj)) {
      seurat_obj <- AddModuleScore(seurat_obj, features = gene_sets_filtered,ctrl = 3,name = NULL)
      results[["AddModuleScore"]] <- seurat_obj@meta.data[, names(gene_sets_filtered), drop = FALSE]
    } else {
      warning("AddModuleScore requires a Seurat object; skipping.")
    }
  }
  
  # Run UCell if requested
  if ("UCell" %in% methods_to_use) {
    seurat_obj <- AddModuleScore_UCell(seurat_obj, features = gene_sets_filtered, name=NULL)
    results[["UCell"]] <- seurat_obj@meta.data[, names(gene_sets_filtered), drop = FALSE]
  }
  
  # Run AUCell if requested
  if ("AUCell" %in% methods_to_use) {
    AUScore<-AUCell_run(expr_matrix, gene_sets_filtered)
    results[["AUCell"]] <- as.data.frame(t(getAUC(AUScore)))
  }
  
  # Run custom average score if requested
  if ("Custom_Average" %in% methods_to_use) {
    custom_score <- sapply(gene_sets_filtered, function(gs) {
      if (length(gs) == 0) {
        return(rep(NA, ncol(expr_matrix)))
      }
      colMeans(expr_matrix[gs, , drop=FALSE])
    })
    results[["Custom_Average"]] <- as.data.frame(custom_score)
  }
  
  # Run singscore if requested
  if ("singscore" %in% methods_to_use) {
    ranked <- rankGenes(expr_matrix)
    singscore <- multiScore(ranked, upSetColc = gene_sets_filtered)
    singscore <-  as.data.frame(t(singscore[["Scores"]]))
    #colnames(singscore)<-names(gene_sets_filtered)
    results[["singscore"]] <- singscore
  }
  
  # Run GSVA if requested
  if ("GSVA" %in% methods_to_use) {
    gsva_params <- gsvaParam(expr_matrix, gene_sets_filtered, maxDiff=FALSE)
    gsva_scores <- gsva(gsva_params)
    results[["GSVA"]] <- as.data.frame(t(gsva_scores))
  }
  
  return(results)
}

# Example:
# gene_sets <- list(GeneSet1 = c("GeneA", "GeneB"), GeneSet2 = c("GeneC", "GeneD"))
# single_method_res <- CalculateGeneSetScores(seurat_obj, gene_sets, methods_to_use = "AUCell")
# multiple_methods_res <- CalculateGeneSetScores(expr_matrix, gene_sets, methods_to_use = c("GSVA", "singscore"))
