# ::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO HOUSE REQUIRED PACKAGES                   #
# ::::::::::::::::::::::::::::::::::::::::::::::::::: #

# SOURCE THIS SCRIPT AT THE TOP OF KEY SCRIPTS THAT WILL NEED ANY PACKAGE
# ALL PACKAGES AVAILABLE ON CRAN
# USE `install.packages("pkgName")` to install them

suppressPackageStartupMessages({
  library(glmmTMB)  # for efficient fitting of complex GLMs
  library(DHARMa)   # for standardized residual diagnostics
  library(stringr)  # for misc. string manipulations
  library(lme4)     # for bootMer() - provides nice interface for parametric bootstrap
  library(boot)     # for boot.ci() - summarizes output of lme4::bootMer()
  library(scales)   # for transparent colors
  library(snow)     # for parallel computing - makes parametric bootstrap much faster
  library(rlecuyer) # for random number generation in parallel computing
})
