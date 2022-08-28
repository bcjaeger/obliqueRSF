
.onAttach <- function(libname, pkgname){
  
  packageStartupMessage(
    "obliqueRSF has been superseded by aorsf. ",
    "I highly recommend you use aorsf to fit oblique ",
    "random survival forests: see https://github.com/bcjaeger/aorsf ",
    "or install from CRAN with install.packages('aorsf')."
  )
  
}