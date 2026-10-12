# Lightweight real-R check of the production CSV reader only; no packages/MA/DR.
code <- as.list(parse("modules/02_scTenifoldOE/scripts/compare_networks.R"))
reader <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                 identical(x[[2]], as.name("read_network_csv")), code)
stopifnot(length(reader) == 1L)
eval(reader[[1]])
for (ids in list(c("001", "002"), c("001", "002", "NA"))) {
  expected <- matrix(c(0, .5, -1.25, .125, 0, 2.5, -.25, .75, 0), nrow = 3)
  expected <- expected[seq_along(ids), seq_along(ids), drop = FALSE]
  dimnames(expected) <- list(ids, ids)
  path <- tempfile(fileext = ".csv")
  write.csv(expected, path)
  actual <- read_network_csv(path)
  stopifnot(identical(dimnames(actual), dimnames(expected)),
            is.numeric(actual), identical(actual, expected))
  cat("PASS IDs and exact numeric weights:", paste(ids, collapse = ","), "\n")
  unlink(path)
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  # Optional Python-exported fixture uses the same full three-ID matrix above.
  actual <- read_network_csv(args[1])
  stopifnot(is.numeric(actual), identical(actual, expected))
  cat("PASS Python-exported CSV -> R: 001,002,NA and exact numeric weights\n")
}
cat("RUN_OK: production R CSV reader; no perturbation or MA/DR performed\n")
