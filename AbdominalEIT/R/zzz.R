.onLoad <- function(libname, pkgname){
  rm(list = ls())
  suppressMessages(FALSE)
  suppressPackageStartupMessages(FALSE)
  packageStartupMessage("abdominal electrical impedance tomography", appendLF = TRUE)
  packageStartupMessage("Type eit() to pop-up a menu", appendLF = TRUE)
}
