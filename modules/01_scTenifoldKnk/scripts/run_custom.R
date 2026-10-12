# Generic bounded KO for an already audited one-cell-type gene panel.
# Usage: Rscript run_custom.R counts.mtx genes.tsv cells.tsv NEW_OUTPUT target1[,target2] [config.json]
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<5) stop("Need counts.mtx genes.tsv cells.tsv NEW_OUTPUT target[,target] [config.json]")
.libPaths(c(Sys.getenv("KNK_R_LIB","work/knk/Rlib"),.libPaths()))
gene_ids <- readLines(args[2]); cell_ids <- readLines(args[3])
tiny <- length(gene_ids)<=100 && length(cell_ids)<=200
minimum_gib <- if(tiny) 0.25 else 0.8
available_before_packages <- ps::ps_system_memory()[["avail"]]
if(available_before_packages<minimum_gib*1024^3) stop("Need >=",minimum_gib," GiB available RAM before starting this scale of run")
library(scTenifoldKnk)
library(Matrix)
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
available_before_counts <- ps::ps_system_memory()[["avail"]]
message("Available RAM after packages: ",round(available_before_counts/1024^3,3),
        " GiB; minimum for this scale: ",minimum_gib," GiB")
if(available_before_counts<minimum_gib*1024^3) stop("Need >=",minimum_gib," GiB available RAM before loading counts")
if(dir.exists(args[4]) && length(list.files(args[4],all.files=TRUE,no..=TRUE))) {
  stop("Output directory must be new or empty; existing analysis is never overwritten")
}
X <- as(readMM(args[1]),"dgCMatrix")
rownames(X) <- gene_ids; colnames(X) <- cell_ids
targets <- unique(strsplit(args[5],",",fixed=TRUE)[[1]])
stopifnot(!anyDuplicated(rownames(X)),!anyDuplicated(colnames(X)),
          all(is.finite(X@x)),all(X@x>=0),all(X@x==round(X@x)),
          all(colSums(X)>0),all(targets %in% rownames(X)))
if(ncol(X)<100 || nrow(X)<50) stop("Need at least 100 cells and 50 audited genes")
if(nrow(X)>3000) stop("Bounded default permits <=3000 genes; preselect an audited panel before running")
detected <- rowMeans(X[targets,,drop=FALSE]>0)
if(any(detected<0.05)) stop("Target detection below 5%; inspect before model fitting")
params <- list(qc=FALSE,nc_nNet=10,nc_nCells=min(500,ncol(X)),nc_nComp=3,
               nc_q=0.9,nc_lambda=0,td_K=3,td_maxIter=1000,td_maxError=1e-5,
               ma_nDim=2,dr_direction=TRUE,dr_directionRegressLibSize=TRUE,nCores=1,seed=1)
if(length(args)>=6) {
  configured <- jsonlite::fromJSON(args[6])
  keys <- intersect(names(params),names(configured))
  params[keys] <- configured[keys]
}
if(params$nc_nCells>ncol(X)) stop("nc_nCells exceeds available cells")
if(params$nCores>4 || params$nc_nNet>20) stop("Bounded runner permits <=4 cores and <=20 networks")
raw_tensor_bytes <- 8*nrow(X)^2*params$nc_nNet
message("Raw double tensor estimate: ",raw_tensor_bytes," bytes; excludes copies and MA objects")
dir.create(args[4],recursive=TRUE,showWarnings=FALSE)
started <- proc.time()
out <- do.call(scTenifoldKnk,c(list(countMatrix=X,gKO=targets),params))
WT <- as.matrix(out$tensorNetworks$WT); KO <- WT; KO[targets,] <- 0
saveRDS(out,file.path(args[4],"output.rds"))
saveRDS(X,file.path(args[4],"counts.rds"))
write.csv(WT,file.path(args[4],"WT.csv")); write.csv(KO,file.path(args[4],"KO.csv"))
write.csv(out$diffRegulation,file.path(args[4],"diffRegulation.csv"),row.names=FALSE)
if(!is.null(out$manifoldAlignment)) write.csv(out$manifoldAlignment,file.path(args[4],"manifoldAlignment.csv"))
summary <- list(status="RUN_OK",scientific_status="UNVALIDATED",cells=ncol(X),genes=nrow(X),
        targets=targets,detection_rate=setNames(as.numeric(detected),targets),
        target_outgoing_edges=setNames(rowSums(WT[targets,,drop=FALSE]!=0),targets),
        parameters=params,elapsed_sec=unname((proc.time()-started)[3]),
        resource_check=list(minimum_available_gib=minimum_gib,
            available_before_packages_bytes=available_before_packages,
            available_before_counts_bytes=available_before_counts,
            raw_tensor_bytes=raw_tensor_bytes),
        interpretation="Network perturbation only; inspect stability and independent validation before biological conclusions")
jsonlite::write_json(summary,file.path(args[4],"run_summary.json"),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(args[4],"sessionInfo.txt"))
print(summary)
