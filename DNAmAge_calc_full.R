# from 
# https://www.bioconductor.org/packages/devel/bioc/vignettes/methylclock/inst/doc/methylclock.html

install.packages(c("tidyverse", "impute", "Rcpp"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("methylclock")

# using this guide 
# https://bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.html
# calculate beta

library(minfi)
library(minfiData)
library(sva)

baseDir <- "../../RAW"
targets <- read.metharray.sheet(baseDir)
RGSet <- read.metharray.exp(targets = targets)
phenoData <- pData(RGSet)
manifest <- getManifest(RGSet)
MSet <- preprocessRaw(RGSet) 
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
beta <- getBeta(RSet)

# from 
# https://www.bioconductor.org/packages/devel/bioc/vignettes/methylclock/inst/doc/methylclock.html
# calculate DNAmAge

library(methylclockData)
library(methylclock)
library(Biobase)
library(tibble)
library(impute)
library(ggplot2)
library(ggpmisc)
library(GEOquery)

myDNAmAge <- DNAmAge(beta)

# read the metadata
dataOnly = read.csv("../../RAW/Sample_Sheet_Full.csv", skip = 8, header = F)
headers = read.csv("../../RAW/Sample_Sheet_Full.csv", skip = 7, header = F, nrows = 1, as.is = T)
colnames(dataOnly) = headers
metadata_file <- dataOnly

# bind metadata with the DNAmAge
myDNAmAge_with_metadata <- cbind(myDNAmAge, metadata_file)

# remove unneeded columns
drop_cols <- c("PedBE", "Wu", "TL", "BNN")
myDNAmAge_with_metadata_without_pedbe_wu_tl_bnn <- myDNAmAge_with_metadata[ , !(names(myDNAmAge_with_metadata) %in% drop_cols)]
write.csv(myDNAmAge_with_metadata_without_pedbe_wu_tl_bnn, "myDNAmAge_with_metadata_without_pedbe_wu_tl_bnn.csv")