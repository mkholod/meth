# using this guide 
# https://bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.html

library(minfi)
library(minfiData)

################
# SAMPLE
RGsetEx # this is the minfi example - should return values if the libraries are loaded correctly
MsetEx <- preprocessRaw(RGsetEx)
GMsetEx <- mapToGenome(MsetEx)

baseDir <- system.file("extdata", package = "minfiData")
list.files(baseDir)

list.files(file.path(baseDir, "5723646052"))
targets <- read.metharray.sheet(baseDir)

sub(baseDir, "", targets$Basename)
RGset <- read.metharray.exp(targets = targets)
pd <- pData(RGset)
pd[,1:4]

RGset2 <- read.metharray.exp(file.path(baseDir, "5723646052"))
RGset3 <- read.metharray.exp(baseDir, recursive = TRUE)

targets2 <- read.csv(file.path(baseDir, "SampleSheet.csv"), 
                     stringsAsFactors = FALSE, skip = 7)

targets2$Basename <- file.path(baseDir, targets2$Sentrix_ID, 
                               paste0(targets2$Sentrix_ID, 
                                      targets2$Sentrix_Position))

annotation(RGsetEx)

#---------- playing on my own with the example
preProcIll <- preprocessIllumina(RGset)

# preProcIll@assays@data@listData[["Meth"]] holds the following
#                 5723646052_R02C02 5723646052_R04C01 5723646052_R05C02 5723646053_R04C02
#cg00050873            20546.488739      2.031322e+02      1.964028e+04      3.741370e+00
#cg00212031              434.955911      1.870509e+02      9.351586e+01      9.353426e-01
#cg00213748             1320.872604      6.178603e+01      3.545810e+02      0.000000e+00

#--------- trying Tirosh methods from methylation 20190913.r
rgSet <- RGset
gset.funnorm <- preprocessFunnorm(rgSet)

#--------- from here https://www.bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
library(minfi)
library(minfiData)
library(sva)

#in the guide but I didn't use
#source("http://www.bioconductor.org/biocLite.R")
#biocLite(c("minfiData", "sva"))
#
#library(devtools)
#install_github("hansenlab/tutorial.450k")

baseDir <- system.file("extdata", package="minfiData")
#targets <- read.450k.sheet(baseDir) # defunct
targets <- read.metharray.sheet(baseDir)

# RGSet <- read.450k.exp(targets = targets) # defunct
RGSet <- read.metharray.exp(targets = targets)

phenoData <- pData(RGSet)
phenoData[,1:6]

manifest <- getManifest(RGSet)
manifest

head(getProbeInfo(manifest))

MSet <- preprocessRaw(RGSet) 
MSet
head(getMeth(MSet)[,1:3])
head(getUnmeth(MSet)[,1:3])

RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
RSet

beta <- getBeta(RSet)
beta

# results with the following - as expected: values between 0 and 1
#                 5723646052_R02C02 5723646052_R04C01 5723646052_R05C02 5723646053_R04C02 5723646053_R05C02
#cg00050873             0.918911031       0.575905975       0.952967421      0.5422818792       0.444444444
#cg00212031             0.093706873       0.654775604       0.140345269      0.4635838150       0.383292383
#cg00213748             0.808383234       0.477324263       0.705588822      0.3175510204       0.240165631
#cg00214611             0.084430237       0.770388959       0.171671672      0.5501355014       0.463611860
#cg00455876             0.781546990       0.334453782       0.754459115      0.4835355286       0.450693374
#cg01707559             0.091942072       0.390332326       0.079120606      0.3731517510       0.419638495

#------- now trying with our own data

library(minfi)
library(minfiData)
library(sva)

baseDir <- "./samples"
targets <- read.metharray.sheet(baseDir)
RGSet <- read.metharray.exp(targets = targets)
phenoData <- pData(RGSet)
manifest <- getManifest(RGSet)
MSet <- preprocessRaw(RGSet) 
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
beta <- getBeta(RSet)