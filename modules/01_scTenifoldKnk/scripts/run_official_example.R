args <- commandArgs(trailingOnly = TRUE)
module <- if (length(args)) args[1] else "modules/01_scTenifoldKnk"
r_lib <- Sys.getenv("KNK_R_LIB", "work/knk/Rlib")
.libPaths(c(r_lib, .libPaths()))
library(scTenifoldKnk)
dir.create(file.path(module, "objects"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(module, "results"), recursive = TRUE, showWarnings = FALSE)
set.seed(1)
X <- matrix(rnbinom(100 * 2000, size = 20, prob = 0.98), ncol = 2000)
rownames(X) <- c(paste0("ng", 1:90), paste0("mt-", 1:10))
colnames(X) <- paste0("cell", seq_len(ncol(X)))
started <- proc.time()
out <- scTenifoldKnk(X, gKO = "ng10", nc_nNet = 10, nc_nCells = 500,
                    td_K = 3, qc_minLibSize = 30, nCores = 1, seed = 1)
saveRDS(X, file.path(module, "objects/official_counts.rds"))
saveRDS(out, file.path(module, "objects/official_example.rds"))
write.csv(out$diffRegulation, file.path(module, "results/official_example_DR.csv"), row.names = FALSE)
write.csv(out$manifoldAlignment, file.path(module, "results/official_example_MA.csv"))
WT <- as.matrix(out$tensorNetworks$WT)
KO <- WT; KO["ng10", ] <- 0
stopifnot(all(KO["ng10", ] == 0), all(WT[-match("ng10", rownames(WT)), ] == KO[-match("ng10", rownames(WT)), ]))
stopifnot(nrow(out$diffRegulation) == nrow(WT), all(is.finite(out$diffRegulation$distance)))
# Official heat mode can screen all 100 simulated genes without unlimited scaling.
screen <- scTenifoldKnk(X, transcriptomeWide = TRUE, nc_nNet = 10,
                       nc_nCells = 500, td_K = 3, qc_minLibSize = 30,
                       nCores = 1, seed = 1)
saveRDS(screen, file.path(module, "objects/official_screen.rds"))
write.csv(screen$perturbationDistances, file.path(module, "results/official_screen_distance.csv"))
write.csv(screen$perturbationDirections, file.path(module, "results/official_screen_direction.csv"))
stopifnot(identical(as.matrix(screen$tensorNetworks$WT), WT))
qcX <- scQC(X,minLibSize=30,removeOutlierCells=TRUE,minPCT=0.05,maxMTratio=0.1)
summary <- list(status = "RUN_OK", input = dim(X), qc_counts = dim(qcX), qc_network_genes = nrow(WT),
                knockout = "ng10", nc_nNet = 10, nc_nCells = 500, td_K = 3,
                ma_nDim = 2, seed = 1, nCores = 1,
                significant_including_target = sum(out$diffRegulation$p.adj < 0.05),
                screen_dimension = dim(screen$perturbationDistances),
                identical_WT_same_seed = TRUE, elapsed_sec = unname((proc.time() - started)[3]),
                scTenifoldKnk = as.character(packageVersion("scTenifoldKnk")),
                scTenifoldNet = as.character(packageVersion("scTenifoldNet")))
jsonlite::write_json(summary, file.path(module, "results/official_run_summary.json"), pretty = TRUE, auto_unbox = TRUE)
capture.output(sessionInfo(), file = file.path(module, "reports/R_sessionInfo.txt"))
print(summary)
