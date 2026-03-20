.onAttach <- function(libname, pkgname) {
  if (interactive()) {
    packageStartupMessage("mirBottleneck loaded. By Mark")
  }
}
