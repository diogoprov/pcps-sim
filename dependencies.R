# Packages used indirectly (not called in R/ or analysis/), declared here so
# that renv::snapshot() records them.
library(qs2)        # targets format = "qs"
library(crew)       # targets controller
library(tarchetypes)
library(BiocManager) # renv: Bioconductor dependencies of phylofactor (tests)
library(data.table)  # long table for phylofactor::gpf (tests)
