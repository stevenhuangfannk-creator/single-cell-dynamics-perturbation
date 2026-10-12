args <- commandArgs(trailingOnly = TRUE)
module <- if (length(args)) args[1] else "modules/01_scTenifoldKnk"
.libPaths(c(Sys.getenv("KNK_R_LIB", "work/knk/Rlib"), .libPaths()))
library(scTenifoldKnk)
library(Matrix)
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
objects <- file.path(module, "objects"); results <- file.path(module, "results")
config <- jsonlite::fromJSON(file.path(module,"config/pbmc.json"))
X <- as(readMM(file.path(objects, "pbmc_cd14_counts.mtx")), "dgCMatrix")
rownames(X) <- readLines(file.path(objects, "pbmc_cd14_genes.tsv"))
colnames(X) <- readLines(file.path(objects, "pbmc_cd14_cells.tsv"))
targets <- readLines(file.path(objects, "candidate_targets.txt"))
stopifnot(nrow(X) >= 100, ncol(X) >= 200, !anyDuplicated(rownames(X)),
          !anyDuplicated(colnames(X)), all(X@x >= 0), all(is.finite(X@x)), all(targets %in% rownames(X)))
params <- config[intersect(names(config),names(formals(scTenifoldKnk)))]
params$dr_direction <- TRUE
started <- proc.time()
out <- do.call(scTenifoldKnk, c(list(countMatrix = X, gKO = "SPI1"), params))
WT <- as.matrix(out$tensorNetworks$WT)
KO <- WT; KO["SPI1", ] <- 0
stopifnot(identical(rownames(WT), rownames(X)), all(diag(WT) == 0),
          all(is.finite(WT)), sum(WT["SPI1", ] != 0) > 0)
saveRDS(X, file.path(objects, "pbmc_counts.rds"))
saveRDS(out, file.path(objects, "single_SPI1.rds"))
saveRDS(WT, file.path(objects, "WT.rds"))
saveRDS(KO, file.path(objects, "KO.rds"))
write.csv(WT, file.path(objects, "WT.csv"))
write.csv(KO, file.path(objects, "KO.csv"))
write.csv(out$diffRegulation, file.path(results, "single_SPI1.csv"), row.names = FALSE)
write.csv(out$manifoldAlignment, file.path(results, "single_SPI1_MA.csv"))
edgeIdx <- which(WT != 0, arr.ind = TRUE)
edges <- data.frame(regulator = rownames(WT)[edgeIdx[,1]], target = colnames(WT)[edgeIdx[,2]], weight = WT[edgeIdx])
write.csv(edges, file.path(results, "WT_edges.csv"), row.names = FALSE)
koEdges <- edges[edges$regulator != "SPI1", ]
write.csv(koEdges, file.path(results, "KO_edges.csv"), row.names = FALSE)
direction_unadjusted <- knockoutDirection(X, gKO = list("SPI1"), genes = rownames(WT), regressLibSize = FALSE)
direction_adjusted <- knockoutDirection(X, gKO = list("SPI1"), genes = rownames(WT), regressLibSize = TRUE)
write.csv(data.frame(gene = colnames(direction_unadjusted), raw_score = as.numeric(direction_unadjusted), library_adjusted_score = as.numeric(direction_adjusted)), file.path(results, "direction_sensitivity.csv"), row.names = FALSE)
audit <- data.frame(gene = targets, outgoing_edges = rowSums(WT[targets,,drop=FALSE] != 0), positive_outgoing_edges = rowSums(WT[targets,,drop=FALSE] > 0), incoming_edges = colSums(WT[,targets,drop=FALSE] != 0))
write.csv(audit, file.path(results, "target_network_audit.csv"), row.names = FALSE)

# Joint KO uses the exact same denoised WT network, official MA and DR functions.
multi_genes <- c("SPI1", "CEBPB")
multiKO <- WT; multiKO[multi_genes, ] <- 0
set.seed(params$seed)
multiMA <- scTenifoldNet::manifoldAlignment(WT, multiKO, d = params$ma_nDim, nCores = params$nCores)
multiDirection <- knockoutDirection(X, gKO = list(multi_genes), genes = rownames(WT), regressLibSize = TRUE)
multiDR <- dRegulation(multiMA, gKO = multi_genes, direction = setNames(as.numeric(multiDirection), colnames(multiDirection)))
saveRDS(list(tensorNetworks = list(WT = Matrix(WT), KO = Matrix(multiKO)), manifoldAlignment = multiMA, diffRegulation = multiDR), file.path(objects, "multi_SPI1_CEBPB.rds"))
write.csv(multiKO, file.path(objects, "multi_KO.csv"))
write.csv(multiDR, file.path(results, "multi_SPI1_CEBPB.csv"), row.names = FALSE)

# Bounded candidate screen in official transcriptomeWide API. Same seeded WT.
screen <- do.call(scTenifoldKnk, c(list(countMatrix = X, gKO = targets, transcriptomeWide = TRUE), params))
stopifnot(identical(WT, as.matrix(screen$tensorNetworks$WT)))
saveRDS(screen, file.path(objects, "candidate_screen.rds"))
write.csv(screen$perturbationDistances, file.path(results, "candidate_screen_distance.csv"))
write.csv(screen$perturbationDirections, file.path(results, "candidate_screen_direction.csv"))
# Reference is an internal predicted SPI1 heat-KO profile, not a disease signature.
reference <- screen$perturbationDistances["SPI1", ] * screen$perturbationDirections["SPI1", ]
reference["SPI1"] <- 0
pm <- perturbationMap(screen, signature = reference, genes = targets, plot = FALSE)
write.csv(pm, file.path(results, "candidate_perturbation_map.csv"), row.names = FALSE)

# Two independent axes: RNG seed and cells subsampled per network. Preserve all runs.
stability <- list(); records <- list()
referenceDR <- out$diffRegulation[!out$diffRegulation$gene %in% "SPI1", ]
referenceDistance <- setNames(referenceDR$distance, referenceDR$gene)
topReference <- head(referenceDR$gene[order(-referenceDR$distance)], 30)
for (n in config$stability_nCells) for (s in config$stability_seeds) {
  label <- paste0("seed", s, "_cells", n)
  t0 <- proc.time()
  p <- params; p$nc_nCells <- n; p$seed <- s
  run <- if (n == params$nc_nCells && s == params$seed) out else do.call(scTenifoldKnk, c(list(countMatrix = X, gKO = "SPI1"), p))
  dr <- run$diffRegulation
  write.csv(dr, file.path(results, paste0("stability_", label, ".csv")), row.names = FALSE)
  dist <- setNames(dr$distance, dr$gene)[names(referenceDistance)]
  top <- head(names(sort(dist, decreasing=TRUE)), 30)
  records[[label]] <- data.frame(run=label, seed=s, cells=n,
       spearman=cor(dist, referenceDistance, method="spearman"),
       top30_jaccard=length(intersect(top,topReference))/length(union(top,topReference)),
       significant_non_target=sum(dr$p.adj < 0.05 & dr$gene != "SPI1"),
       elapsed_sec=unname((proc.time()-t0)[3]))
  stability[[label]] <- dist
}
write.csv(do.call(rbind,records), file.path(results,"stability_summary.csv"), row.names=FALSE)
write.csv(do.call(cbind,stability), file.path(results,"stability_distances.csv"))
# No-change network is a negative software control; all p values must be one.
set.seed(params$seed)
nullMA <- scTenifoldNet::manifoldAlignment(WT, WT, d=params$ma_nDim, nCores=params$nCores)
nullDR <- dRegulation(nullMA)
write.csv(nullDR, file.path(results,"negative_control_WT_vs_WT.csv"), row.names=FALSE)
stopifnot(all(nullDR$p.value == 1))
saveRDS(nullMA, file.path(objects,"negative_control_MA.rds"))

summary <- list(status="PARTIAL", run_status="RUN_OK", dataset="10x PBMC3k",
                cell_type="CD14+ Monocytes", cells=ncol(X), genes=nrow(X),
                parameters=params, edge_convention="WT[row regulator, column target]; official final transpose retained",
                total_nonzero_edges=sum(WT!=0), single_target_outgoing_edges=sum(WT["SPI1",]!=0),
                significant_non_target=sum(out$diffRegulation$p.adj<0.05 & out$diffRegulation$gene!="SPI1"),
                multi_significant_non_targets=sum(multiDR$p.adj<0.05 & !multiDR$gene %in% multi_genes),
                screen_targets=targets, screen_dimensions=dim(screen$perturbationDistances),
                negative_control_all_p_one=TRUE, tensor_raw_bytes=8*nrow(X)^2*params$nc_nNet,
                elapsed_sec=unname((proc.time()-started)[3]),
                scope="Bounded real-data validation on selected 600-gene network and 9 candidate knockouts; not original paper replication or full-transcriptome biological validation")
jsonlite::write_json(summary,file.path(results,"pbmc_run_summary.json"),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(module,"reports/R_sessionInfo.txt"))
print(summary)
