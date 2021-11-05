if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("minfi")
BiocManager::install("minfiData")
BiocManager::install("GEOquery")
BiocManager::install("ChAMP")
BiocManager::install("DMRcate")

# source("http://www.bioconductor.org/biocLite.R")
# biocLite(c("minfi", "GEOquery"))

# all of the above should be run only once

library(minfi)
library(minfiData)
library(GEOquery)

# following this 
# https://bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.html

# this just displays an example
RGsetEx 

# loads sample data
MsetEx <- preprocessRaw(RGsetEx)
GMsetEx <- mapToGenome(MsetEx)

list.files("RAW")

# Now the real data
# baseDir <- system.file(path = "./RAW", package = "minfiData")
# baseDir <- system.file("/Users/michaelk/Documents/Shiri-R/RAW", package = "minfiData")
# list.files(baseDir)
# list.files(file.path("/Users/michaelk/Documents/Shiri-R/RAW"))
file.ls <- list.files(path="RAW", pattern="GSM3936974")
file.ls # result: [1] "GSM3936974_202184900129_R08C01_Grn.idat" "GSM3936974_202184900129_R08C01_Red.idat"
# read.metharray.sheet("/Users/michaelk/Documents/Shiri-R/RAW")
# file.ls <- list.files(path="../Shiri-R", pattern="GSM3936974")
file.ls <- list.files(path="./")
targets <- read.metharray.sheet("./")
tutorial_csv <-read.csv(file = "TUTORIAL/datMiniAnnotation27k.csv")
View(tutorial_csv)

# rgSet <- read.450k.exp("RAW/idat")
rgSet <- read.metharray.exp("./RAW_SAMPLE")
pData(rgSet)
head(sampleNames(rgSet))
dim(getRed(rgSet))
dim(getGreen(rgSet))
mset <- preprocessIllumina(rgSet)

library(DMRcate)
library(ChAMP)

# targets <- read.metharray.sheet(baseDir, pattern = "Sample_Sheet.csv")

# baseDir <- system.file("extdata", package="minfiData")
targets <- read.450k.sheet(baseDir)
targets <- read.metharray.sheet(baseDir)

targets <- read.metharray.sheet("./RAW_SAMPLE")
RGSet <- read.metharray.exp(base = "./RAW_SAMPLE", targets = targets)

beta <- getBeta(gset.funnorm) #gets matrix
beta <- getBeta(gset.funnorm) #gets matrix
####
setwd("/Users/michaelk/Documents/Shiri-R")
targets <- read.metharray.sheet("./RAW_SAMPLE")
rgSet <- read.metharray.exp("./RAW_SAMPLE")
mset <- preprocessIllumina(rgSet)
phenoData <- pData(rgSet)
mphenoData <- pData(mset)

################################################################################################
# https://www.bioconductor.org/packages/devel/bioc/vignettes/methylclock/inst/doc/methylclock.html
install.packages(c("tidyverse", "impute", "Rcpp"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')
BiocManager::install("methylclock")
BiocManager::install("methylclockData", force = TRUE)

library(methylclockData)
library(methylclock)

library(Biobase)
library(tibble)
library(impute)
library(ggplot2)
library(ggpmisc)
library(GEOquery) 

TestDataset <- get_TestDataset()
cpgs.missing <- checkClocks(TestDataset)
cpgs.missing.GA <- checkClocksGA(TestDataset)