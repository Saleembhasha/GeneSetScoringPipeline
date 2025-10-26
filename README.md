# GeneSetScoringPipeline

A flexible R pipeline for scoring gene sets (modules) in single-cell RNA-sequencing data, supporting Seurat, AUCell, UCell, singscore, GSVA, and custom averaging methods. This tool enables you to compute per-cell gene set enrichment using different algorithms in a single unified workflow.

## Features

- Accepts either a Seurat object or expression matrix as input
- Supports the following methods:
  - Seurat AddModuleScore
  - AUCell
  - UCell
  - GSVA
  - singscore
  - Custom average expression
- Flexible: specify one or more methods to run
- Returns a named list of dataframes with scores

## Installation

You need to have R (>=4.0.0) and install the required packages:
install.packages(c("Seurat", "GSVA", "Matrix"))
These may need Bioconductor

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(c("AUCell", "UCell", "singscore"))

## Usage

Place your gene sets in a list format:
gene_sets <- list(
GeneSetA = c("Gene1", "Gene2"),
GeneSetB = c("Gene3", "Gene4")
)

Run the pipeline:
source("src/gene_set_scoring_pipeline.R")

scores <- CalculateGeneSetScores(
seurat_obj = your_seurat_object, # or expr_matrix = your_expression_matrix
gene_sets = gene_sets,
methods_to_use = c("GSVA", "singscore")
)

Scores will be in a list (each element corresponds to method used).

## Example

See `example/demo_example.R` for a demonstration with toy data.

## License

MIT
