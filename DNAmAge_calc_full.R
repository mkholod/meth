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
cpgs.missing <- checkClocks(myDNAmAge_with_acceleration_age)
cpgs.missing.GA <- checkClocksGA(myDNAmAge_with_acceleration_age)

# bind metadata with the DNAmAge
myDNAmAge_with_acceleration_age_with_metadata <- cbind(myDNAmAge_with_acceleration_age, metadata_file)

write.csv(myDNAmAge_with_acceleration_age_with_metadata, "csv/myDNAmAge_with_acceleration_age_with_metadata.csv")

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

is_filter_negative_calculated_age = TRUE
out_file_avg = "any_met_avg_yes_no_by_clocks.csv"
out_file_sd = "any_met_sd_yes_no_by_clocks.csv"

if (is_filter_negative_calculated_age) {
  filter_AnyMetYes <- subset(filter_AnyMetYes, Levine > 0)
  filter_AnyMetNo <- subset(filter_AnyMetNo, Levine > 0)
  out_file_avg = "any_met_avg_yes_no_by_clocks_without_negative_levine.csv"
  out_file_sd = "any_met_sd_yes_no_by_clocks_without_negative_levine.csv"
}

avg_by_clock <- c()
sd_by_clock <- c()
for (clock_type in c("ageAcc.Horvath", "ageAcc2.Horvath", "ageAcc3.Horvath", "ageAcc.Hannum", "ageAcc2.Hannum", "ageAcc3.Hannum", "ageAcc.Levine", "ageAcc2.Levine", "ageAcc3.Levine", "ageAcc.Hovarth2", "ageAcc2.Hovarth2", "ageAcc3.Hovarth2")) {
    avg_yes <- mean(filter_AnyMetYes[[clock_type]], na.rm=TRUE)
    avg_no <- mean(filter_AnyMetNo[[clock_type]], na.rm=TRUE)
    avg_by_clock <- rbind(avg_by_clock, c(clock_type, avg_yes, avg_no))
    
    sd_yes <- sd(filter_AnyMetYes[[clock_type]], na.rm=TRUE)
    sd_no <- sd(filter_AnyMetNo[[clock_type]], na.rm=TRUE)
    sd_by_clock <- rbind(sd_by_clock, c(clock_type, sd_yes, sd_no))
}
colnames(avg_by_clock) <- c("clock_type", "Yes_avg", "No_avg")
write.csv(avg_by_clock, out_file_avg)

colnames(sd_by_clock) <- c("clock_type", "Yes_sd", "No_sd")
write.csv(sd_by_clock, out_file_sd)

# To perform a TTEST for acceleration1 for Hovarth and Levine
# And check if there is a difference between Sporadic and Vhl
# using https://www.youtube.com/watch?v=RlhnNbPZC0A&ab_channel=MarinStatsLectures-RProgramming%26Statistics
all_data = myDNAmAge_with_acceleration_age_with_metadata
boxplot(myDNAmAge_with_acceleration_age_with_metadata$ageAcc.Levine ~ myDNAmAge_with_acceleration_age_with_metadata$Sample_Group)
boxplot(myDNAmAge_with_acceleration_age_with_metadata$ageAcc.Horvath ~ myDNAmAge_with_acceleration_age_with_metadata$Sample_Group)

t.test(all_data$ageAcc.Horvath[all_data$Sample_Group == "Sporadic"], all_data$ageAcc.Horvath[all_data$Sample_Group == "VHL"])

#Welch Two Sample t-test

#data:  all_data$ageAcc.Horvath[all_data$Sample_Group == "Sporadic"] and all_data$ageAcc.Horvath[all_data$Sample_Group == "VHL"]
#t = -0.92099, df = 39.473, p-value = 0.3626
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -13.464380   5.036927
#sample estimates:
#  mean of x mean of y 
#8.956855 13.170581 

t.test(all_data$ageAcc.Levine[all_data$Sample_Group == "Sporadic"], all_data$ageAcc.Levine[all_data$Sample_Group == "VHL"])

# Welch Two Sample t-test

#data:  all_data$ageAcc.Levine[all_data$Sample_Group == "Sporadic"] and all_data$ageAcc.Levine[all_data$Sample_Group == "VHL"]
#t = -0.95366, df = 38.685, p-value = 0.3462
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -12.430501   4.466156
#sample estimates:
#  mean of x mean of y 
#-15.76965 -11.78748 

t.test(all_data$ageAcc2.Horvath[all_data$Sample_Group == "Sporadic"], all_data$ageAcc2.Horvath[all_data$Sample_Group == "VHL"])

# Welch Two Sample t-test

# data:  all_data$ageAcc2.Horvath[all_data$Sample_Group == "Sporadic"] and all_data$ageAcc2.Horvath[all_data$Sample_Group == "VHL"]
# t = 0.24433, df = 19.742, p-value = 0.8095
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -9.04976 11.44871
#sample estimates:
#  mean of x  mean of y 
#0.9136663 -0.2858094 

t.test(all_data$ageAcc2.Levine[all_data$Sample_Group == "Sporadic"], all_data$ageAcc2.Levine[all_data$Sample_Group == "VHL"])


#Welch Two Sample t-test

#data:  all_data$ageAcc2.Levine[all_data$Sample_Group == "Sporadic"] and all_data$ageAcc2.Levine[all_data$Sample_Group == "VHL"]
#t = 0.064188, df = 28.41, p-value = 0.9493
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -8.272391  8.807964
#sample estimates:
#  mean of x mean of y 
#3.897955  3.630168 

t.test(all_data$ageAcc.Horvath[all_data$Sample_Group != "VHL"], all_data$ageAcc.Horvath[all_data$Sample_Group == "VHL"])

# Welch Two Sample t-test

#data:  all_data$ageAcc.Horvath[all_data$Sample_Group != "VHL"] and all_data$ageAcc.Horvath[all_data$Sample_Group == "VHL"]
#t = -1.2458, df = 21.332, p-value = 0.2263
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -11.822740   2.959044
#sample estimates:
#  mean of x mean of y 
#8.738733 13.170581 

t.test(all_data$ageAcc.Levine[all_data$Sample_Group != "VHL"], all_data$ageAcc.Levine[all_data$Sample_Group == "VHL"])



# TODO metastasis predicition accel vs features
# TODO does it divide to sample groups unsupervised
# TODO VHL+MEN1 vs Sporadic
# TODO in MEN1 check for only under the age of 50