# from 
# https://www.bioconductor.org/packages/devel/bioc/vignettes/methylclock/inst/doc/methylclock.html

install.packages(c("tidyverse", "impute", "Rcpp"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("methylclock")

library(methylclockData)
library(methylclock)

library(Biobase)
library(tibble)
library(impute)
library(ggplot2)
library(ggpmisc)
library(GEOquery)

# misha
library(sqldf)
x=read.csv.sql("beta1.csv")
toBetas=FALSE # x already contains betas
normalize=FALSE # because it's time consuming
fastImp=TRUE

# dealing with missing data we saw with horvath that we're missing 1059 samples
# Get TestDataset data
TestDataset <- get_TestDataset()
cpgs.missing <- checkClocks(TestDataset)
cpgs.missing.GA <- checkClocksGA(TestDataset)
commonClockCpgs(cpgs.missing, "Hannum" )
commonClockCpgs(cpgs.missing.GA, "Bohlin" )

# 4.1Data in Horvath’s format (e.g. csv with CpGs in rows)
library(tidyverse)
MethylationData <- get_MethylationDataExample()
MethylationData
