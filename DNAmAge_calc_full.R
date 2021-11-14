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

# read real age
# meta <- read.csv("../../RAW/Sample_Sheet_Full.csv", header=TRUE)
dataOnly = read.csv("../../RAW/Sample_Sheet_Full.csv", skip = 8, header = F)
headers = read.csv("../../RAW/Sample_Sheet_Full.csv", skip = 7, header = F, nrows = 1, as.is = T)
colnames(dataOnly) = headers
metadata_file <- dataOnly

#meta_with_data_only <- meta[-c(1:6),]
#real_age <- meta$X.4[8:103]

# add column to myDNAmAge
#myDNAmAge_with_real_age <- myDNAmAge
#myDNAmAge_with_real_age$real_age <- real_age
#write.csv(myDNAmAge_with_real_age, "DNAmAge_full_with_real_age.csv")

# remove PedBE, Wu, TL
#write.csv(myDNAmAge_with_real_age[,c("id", "Horvath", "Hannum", "Levine", "skinHorvath", "real_age")], "DNAmAge_full_with_real_age_without_pedbe_wu_tl.csv")

# add metadata
myDNAmAge_with_metadata <- myDNAmAge
myDNAmAge_with_metadata <- cbind(myDNAmAge_with_metadata, metadata_file)

drop_cols <- c("PedBE", "Wu", "TL")
myDNAmAge_with_metadata_without_pedbe_wu_tl <- myDNAmAge_with_metadata[ , !(names(myDNAmAge_with_metadata) %in% drop_cols)]
write.csv(myDNAmAge_with_metadata_without_pedbe_wu_tl, "myDNAmAge_with_metadata_without_pedbe_wu_tl.csv")