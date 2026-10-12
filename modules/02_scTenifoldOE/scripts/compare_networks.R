# Compare supplied networks with the same official manifold and DR functions.
read_network_csv <- function(path) {
  header <- read.csv(path, nrows = 0, check.names = FALSE)
  data <- read.csv(path, check.names = FALSE, na.strings = NULL,
                   colClasses = c("character", rep("numeric", ncol(header) - 1L)))
  weights <- as.matrix(data[-1])
  rownames(weights) <- data[[1]]
  weights
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop("WT CSV, perturbed CSV, target, output prefix, R library, dimensions, seed required")
.libPaths(c(args[5], .libPaths()))
stopifnot(as.character(packageVersion("scTenifoldKnk")) == "2.0.0")
stopifnot(as.character(packageVersion("scTenifoldNet")) == "1.4.3")
WT <- read_network_csv(args[1])
perturbed <- read_network_csv(args[2])
stopifnot(identical(dimnames(WT), dimnames(perturbed)), all(is.finite(WT)), all(is.finite(perturbed)))
set.seed(as.integer(args[7]))
MA <- scTenifoldNet::manifoldAlignment(WT, perturbed, d = as.integer(args[6]), nCores = 1)
DR <- scTenifoldKnk::dRegulation(MA, gKO = args[3], empiricalNull = FALSE, direction = NULL)
stopifnot(nrow(DR) == nrow(WT), all(is.finite(DR$distance)), all(is.finite(DR$FC)))
write.csv(MA, paste0(args[4], "_alignment.csv"))
write.csv(DR, paste0(args[4], "_dr.csv"), row.names = FALSE)
saveRDS(list(WT = WT, perturbed = perturbed, alignment = MA, diffRegulation = DR), paste0(args[4], ".rds"))
writeLines(capture.output(sessionInfo()), paste0(args[4], "_session.txt"))
