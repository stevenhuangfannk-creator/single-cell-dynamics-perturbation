# Pinned official R packages; execute from the repository root.
lib <- Sys.getenv("KNK_R_LIB", "work/knk/Rlib")
dir.create(lib, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(lib, .libPaths()))
options(repos=c(CRAN="https://cloud.r-project.org"), timeout=300)
dependencies <- c("Matrix","MASS","RSpectra","RhpcBLASctl","cli","ps",
                  "enrichR","igraph","reshape2","rjson","plyr","jsonlite")
missing <- dependencies[!vapply(dependencies,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) install.packages(missing,lib=lib,type=if(.Platform$OS.type=="windows") "binary" else "source")
packages <- c(scTenifoldNet="1.4.3", scTenifoldKnk="2.0.0")
md5 <- c(scTenifoldNet="78b614b5aac290af675fae71f7f06bb1", scTenifoldKnk="b83586c28adcf57225963df5911071de")
for(pkg in names(packages)) {
  if(pkg %in% rownames(installed.packages()) && as.character(packageVersion(pkg))==packages[[pkg]]) next
  filename <- paste0(pkg,"_",packages[[pkg]],".tar.gz")
  path <- file.path(tempdir(),filename)
  url <- paste0("https://cloud.r-project.org/src/contrib/",filename)
  attempt <- try(download.file(url,path,mode="wb"),silent=TRUE)
  if(inherits(attempt,"try-error")) download.file(paste0("https://cran.r-project.org/src/contrib/Archive/",pkg,"/",filename),path,mode="wb")
  stopifnot(unname(tools::md5sum(path))==md5[[pkg]])
  install.packages(path,repos=NULL,type="source",lib=lib)
  stopifnot(requireNamespace(pkg,quietly=TRUE),as.character(packageVersion(pkg))==packages[[pkg]])
}
capture.output(sessionInfo(),file="environments/knk_install_session.txt")
cat("RUN_OK: scTenifoldNet 1.4.3, scTenifoldKnk 2.0.0\n")
