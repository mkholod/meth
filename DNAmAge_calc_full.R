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

#plot 
plotDNAmAge(myDNAmAge_with_acceleration_with_metadata_without_pedbe_wu_tl_bnn$Horvath, metadata_file$Age)
plotDNAmAge(myDNAmAge_with_acceleration_with_metadata_without_pedbe_wu_tl_bnn$Hannum, metadata_file$Age, tit="Hannum")
plotDNAmAge(myDNAmAge_with_acceleration_with_metadata_without_pedbe_wu_tl_bnn$Levine, metadata_file$Age, tit="Levine")
plotDNAmAge(myDNAmAge_with_acceleration_with_metadata_without_pedbe_wu_tl_bnn$skinHorvath, metadata_file$Age, tit="skinHorvath")


# make for each clock by sample group
myData <- myDNAmAge_with_acceleration_with_metadata_without_pedbe_wu_tl_bnn
filter_sporadic <- subset(myData, Sample_Group == "Sporadic")
filter_MEN1 <- subset(myData, Sample_Group == "MEN1")
filter_VHL <- subset(myData, Sample_Group == "VHL")
rownames(filter_sporadic) <- NULL

plotDNAmAge(filter_sporadic$Horvath, filter_sporadic$Age, tit="Hovarth (sporadic)")
plotDNAmAge(filter_MEN1$Horvath, filter_MEN1$Age, tit="Hovarth (MEN1)")
plotDNAmAge(filter_VHL$Horvath, filter_VHL$Age, tit="Hovarth (VHL)")

plotDNAmAge(filter_sporadic$Hannum, filter_sporadic$Age, tit="Hannum (sporadic)")
plotDNAmAge(filter_MEN1$Hannum, filter_MEN1$Age, tit="Hannum (MEN1)")
plotDNAmAge(filter_VHL$Hannum, filter_VHL$Age, tit="Hannum (VHL)")

plotDNAmAge(filter_sporadic$Levine, filter_sporadic$Age, tit="Levine (sporadic)")
plotDNAmAge(filter_MEN1$Levine, filter_MEN1$Age, tit="Levine (MEN1)")
plotDNAmAge(filter_VHL$Levine, filter_VHL$Age, tit="Levine (VHL)")

plotDNAmAge(filter_sporadic$skinHorvath, filter_sporadic$Age, tit="skinHorvath (sporadic)")
plotDNAmAge(filter_MEN1$skinHorvath, filter_MEN1$Age, tit="skinHorvath (MEN1)")
plotDNAmAge(filter_VHL$skinHorvath, filter_VHL$Age, tit="skinHorvath (VHL)")

# for each clock calculate
# 1. mean deviation 2. average and present the value for AnyMet (yes/no)
filter_AnyMetYes <- subset(myData, AnyMet == "Yes")
filter_AnyMetNo <- subset(myData, AnyMet == "No")

avg_by_clock <- c()
for (clock_type in c("Horvath", "Hannum", "Levine", "skinHorvath")) {
    avg_yes <- mean(filter_AnyMetYes[[clock_type]], na.rm=TRUE)
    avg_no <- mean(filter_AnyMetNo[[clock_type]], na.rm=TRUE)
    avg_by_clock <- rbind(avg_by_clock, c(clock_type, avg_yes, avg_no))
}
colnames(avg_by_clock) <- c("clock_type", "Yes_avg", "No_avg")
write.csv(avg_by_clock, "any_met_avg_yes_no_by_clocks.csv")


# TODO metastasis predicition accel vs features
# TODO does it divide to sample groups unsupervised
# TODO VHL+MEN1 vs Sporadic
# TODO in MEN1 check for only under the age of 50
