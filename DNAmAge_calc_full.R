# from 
# https://www.bioconductor.org/packages/devel/bioc/vignettes/methylclock/inst/doc/methylclock.html
# calculate DNAmAge

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

#install.packages(c("tidyverse", "impute", "Rcpp"))

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# BiocManager::install("methylclock")

library(methylclockData)
library(methylclock)

library(Biobase)
library(tibble)
library(impute)
library(ggplot2)
library(ggpmisc)
library(GEOquery)

# library(sqldf)
# myData=read.csv.sql("beta1.csv")
myDNAmAge <- DNAmAge(beta)
write.csv(myDNAmAge, "DNAmAge_full.csv")