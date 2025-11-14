## * Install packages
# sf and excursions are prerquisites for BayesfMRI
install.packages(c("sf", "excursions", "BayesfMRI"))

## * Install INLA
install.packages("INLA",
                 repos = c(getOption("repos"),
                           INLA = "https://inla.r-inla-download.org/R/stable"),
                 dep = FALSE
                 )

## * Install INLA binaries
INLA::inla.binary.install()

## * Load Libraries
library("ciftiTools")
library("BayesfMRI")
library("INLA")

## * Set Workbench installation location for ciftiTools
ciftiTools.setOption(
  "wb_path","~/Software/workbench/2.1.0/bin_linux64"
)
