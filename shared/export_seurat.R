# Usage: Rscript shared/export_seurat.R input.rds data/my_counts RNA
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2)
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Matrix))
obj <- readRDS(args[1])
assay <- if (length(args) >= 3) args[3] else "RNA"
if (packageVersion("SeuratObject") >= "5.0.0") {
  layers <- Layers(obj[[assay]], search = "^counts")
  if (length(layers) != 1) stop("Expected one counts layer. Explicitly JoinLayers after checking sample identities.")
  counts <- LayerData(obj, assay = assay, layer = layers[1])
} else {
  counts <- GetAssayData(obj, assay = assay, slot = "counts")
}
stopifnot(!anyDuplicated(rownames(counts)), !anyDuplicated(colnames(counts)))
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
writeMM(as(counts, "sparseMatrix"), file.path(args[2], "matrix.mtx"))
write.table(data.frame(rownames(counts), rownames(counts), "Gene Expression"), file.path(args[2], "features.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(colnames(counts), file.path(args[2], "barcodes.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(obj[[]][colnames(counts), , drop = FALSE], file.path(args[2], "metadata.tsv"), sep = "\t", quote = FALSE, col.names = NA)
cat(nrow(counts), "genes;", ncol(counts), "cells exported\n")
