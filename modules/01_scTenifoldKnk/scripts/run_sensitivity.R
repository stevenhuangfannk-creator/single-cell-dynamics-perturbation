# Pre-registered method sensitivity: 1200 genes, 10 networks, 400 cells, seeds 1/2/3.
args <- commandArgs(trailingOnly=TRUE)
module <- if(length(args)) args[1] else "modules/01_scTenifoldKnk"
.libPaths(c(Sys.getenv("KNK_R_LIB","work/knk/Rlib"),.libPaths()))
library(scTenifoldKnk)
library(Matrix)
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
input <- file.path(module,"objects/sensitivity1200/objects")
output <- file.path(module,"results/sensitivity1200")
dir.create(output,recursive=TRUE,showWarnings=FALSE)
X <- as(readMM(file.path(input,"pbmc_cd14_counts.mtx")),"dgCMatrix")
rownames(X)<-readLines(file.path(input,"pbmc_cd14_genes.tsv"))
colnames(X)<-readLines(file.path(input,"pbmc_cd14_cells.tsv"))
baseline <- read.csv(file.path(module,"results/single_SPI1.csv"))
ref <- setNames(baseline$distance,baseline$gene)
records <- list(); all_distances <- list()
for(s in 1:3){
  free <- ps::ps_system_memory()[["avail"]]
  if(is.finite(free) && free < 0.8*1024^3) {
    message("BLOCKED remaining seeds: free RAM below 0.8 GiB; completed runs retained")
    break
  }
  available <- scTenifoldNet::checkMemory(nrow(X),nNet=10,nConditions=1)
  # checkMemory logs the package estimate; the 0.8 GiB free-RAM gate precedes each run.
  started <- proc.time()
  out <- scTenifoldKnk(X,gKO="SPI1",qc=FALSE,nc_nNet=10,nc_nCells=400,
       nc_nComp=3,nc_q=0.9,td_K=3,ma_nDim=2,nCores=1,seed=s,
       dr_directionRegressLibSize=TRUE)
  saveRDS(out,file.path(module,paste0("objects/sensitivity1200/run_seed",s,".rds")))
  dr <- out$diffRegulation
  write.csv(dr,file.path(output,paste0("seed",s,"_DR.csv")),row.names=FALSE)
  vals <- setNames(dr$distance,dr$gene)
  common <- setdiff(intersect(names(ref),names(vals)),"SPI1")
  top1 <- head(names(sort(ref[common],decreasing=TRUE)),30)
  top2 <- head(names(sort(vals[common],decreasing=TRUE)),30)
  records[[s]] <- data.frame(seed=s,genes=nrow(X),networks=10,cells=400,
       common_non_target_genes=length(common),spearman_to_600=cor(ref[common],vals[common],method="spearman"),
       top30_jaccard_to_600=length(intersect(top1,top2))/length(union(top1,top2)),
       significant_non_target=sum(dr$p.adj<0.05 & dr$gene!="SPI1"),
       elapsed_sec=unname((proc.time()-started)[3]))
  all_distances[[paste0("seed",s)]] <- vals[setdiff(rownames(X),"SPI1")]
  write.csv(do.call(rbind,records),file.path(output,"parameter_sensitivity_summary.csv"),row.names=FALSE)
  rm(out); gc()
}
if(!length(all_distances)) stop("No seed run completed: RAM gate")
d <- do.call(cbind,all_distances)
write.csv(d,file.path(output,"distances.csv"))
write.csv(cor(d,method="spearman"),file.path(output,"within_1200_seed_spearman.csv"))
print(do.call(rbind,records)); print(cor(d,method="spearman"))
