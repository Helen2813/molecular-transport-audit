options(stringsAsFactors = FALSE)

SCRIPT_VERSION <- "00-install-tcbb-r-packages-v1-no-cli"

cat(strrep("=", 100), "\n", sep = "")
cat("Paper 4 / TCBB - install required R packages into USER library\n")
cat(strrep("=", 100), "\n", sep = "")
cat("Script version: ", SCRIPT_VERSION, "\n", sep = "")

# Use R's standard per-user library. This avoids requiring Administrator rights
# for C:/Program Files/R/R-4.6.1/library.
user_lib <- Sys.getenv("R_LIBS_USER")

if (!nzchar(user_lib)) {
  local_appdata <- Sys.getenv("LOCALAPPDATA")
  r_minor <- strsplit(R.version$minor, "\\.")[[1]][1]
  user_lib <- file.path(
    local_appdata, "R", "win-library",
    paste0(R.version$major, ".", r_minor)
  )
}

dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

cat("R version:    ", R.version.string, "\n", sep = "")
cat("User library: ", user_lib, "\n", sep = "")
cat(".libPaths():\n")
print(.libPaths())

options(repos = c(CRAN = "https://cloud.r-project.org"))

cat("\n[1/4] Installing BiocManager into user library ...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = user_lib, dependencies = TRUE)
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  stop("BiocManager installation failed.")
}

cat("Bioconductor version appropriate for this R: ",
    as.character(BiocManager::version()), "\n", sep = "")

cat("\n[2/4] Installing WGCNA + Bioconductor dependencies ...\n")
BiocManager::install(
  "WGCNA",
  lib = user_lib,
  ask = FALSE,
  update = FALSE
)

if (!requireNamespace("WGCNA", quietly = TRUE)) {
  stop("WGCNA installation failed.")
}

cat("\n[3/4] Installing NetRep ...\n")
install.packages(
  "NetRep",
  lib = user_lib,
  dependencies = TRUE
)

if (!requireNamespace("NetRep", quietly = TRUE)) {
  stop("NetRep installation failed.")
}

cat("\n[4/4] Verification ...\n")
cat("WGCNA version: ", as.character(packageVersion("WGCNA")), "\n", sep = "")
cat("NetRep version: ", as.character(packageVersion("NetRep")), "\n", sep = "")

# Save a reproducibility snapshot.
out_dir <- file.path(
  "D:/paper4_tcbb_data",
  "paper4_tcbb_r_environment_v1"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "R_sessionInfo_after_install_v1.txt")
)

installed <- data.frame(
  package = c("BiocManager", "WGCNA", "NetRep"),
  version = c(
    as.character(packageVersion("BiocManager")),
    as.character(packageVersion("WGCNA")),
    as.character(packageVersion("NetRep"))
  ),
  stringsAsFactors = FALSE
)

write.table(
  installed,
  file = file.path(out_dir, "R_core_package_versions_v1.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("\n", strrep("=", 100), "\n", sep = "")
cat("R PACKAGE INSTALL: PASS\n")
cat(strrep("=", 100), "\n", sep = "")
cat("No R IDE was installed.\n")
cat("Packages were installed into the per-user R library only.\n")
cat("Environment snapshot: ", out_dir, "\n", sep = "")
cat(strrep("=", 100), "\n", sep = "")
