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

is_full_data <- TRUE # full 96 samples if TRUE, 16 samples for speed if FALSE
if (is_full_data) {
  baseDir <- "RAW"
} else {
  baseDir <- "samples"
}

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

# read the metadata
if (is_full_data) {
  sample_sheet_filename = "/Sample_Sheet_Full.csv"
} else {
  sample_sheet_filename = "/Sample_Sheet_16.csv"
}

path <- paste(baseDir, sample_sheet_filename, sep="")
dataOnly = read.csv(path, skip = 8, header = F)
headers = read.csv(path, skip = 7, header = F, nrows = 1, as.is = T)

colnames(dataOnly) = headers
metadata_file <- dataOnly

# calculate ageAcceleration as specified in 4.2 in https://www.bioconductor.org/packages/devel/bioc/vignettes/methylclock/inst/doc/methylclock.html
myDNAmAge_with_acceleration_age <- DNAmAge(beta, age=metadata_file$Age, cell.count=TRUE)

# bind metadata with the DNAmAge
myDNAmAge_with_acceleration_age_with_metadata <- cbind(myDNAmAge_with_acceleration_age, metadata_file)

# remove unneeded columns
drop_cols <- c("PedBE", "Wu", "TL", "BNN")
myDNAmAge_with_acceleration_with_metadata_without_pedbe_wu_tl_bnn <- myDNAmAge_with_acceleration_age_with_metadata[ , !(names(myDNAmAge_with_acceleration_age_with_metadata) %in% drop_cols)]
write.csv(myDNAmAge_with_acceleration_with_metadata_without_pedbe_wu_tl_bnn, "myDNAmAge_with_acceleration_with_metadata_without_pedbe_wu_tl_bnn.csv")




