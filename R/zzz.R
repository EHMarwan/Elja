.onAttach <- function(libname, pkgname) {
  packageStartupMessage("This package has required dplyr, ggplot2, stats and MASS.")
}
