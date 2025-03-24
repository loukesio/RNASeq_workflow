install_and_load <- function(packages, bioc_packages = NULL, github_packages = NULL) {
  new_pkgs <- packages[!packages %in% installed.packages()[,"Package"]]
  if(length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)

  if(!is.null(bioc_packages)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    new_bioc <- bioc_packages[!bioc_packages %in% installed.packages()[,"Package"]]
    if(length(new_bioc)) BiocManager::install(new_bioc)
  }

  if(!is.null(github_packages)) {
    if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
    github_pkg_names <- gsub(".*/(.*)", "\1", github_packages)
    new_github <- github_packages[!github_pkg_names %in% installed.packages()[,"Package"]]
    if(length(new_github)) lapply(new_github, devtools::install_github)
  }

  all_pkgs <- c(packages, bioc_packages, if (!is.null(github_packages)) github_pkg_names)
  invisible(lapply(all_pkgs, library, character.only = TRUE))
}