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
write.csv(beta, "./beta.csv")

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

###
load_DNAm_Clocks_data()
horvath_cgs <- coefHorvath[-1,1]
levine_cgs <- coefLevine[-1,1]
row_names_beta <- row.names(beta)

library(stringr)
row_names_beta_sorted <- str_sort(row_names_beta)

horvath_cgs_sorted <- str_sort(horvath_cgs)
missing <- horvath_cgs_sorted %in% row_names_beta_sorted
missing_horvath_cgs <- horvath_cgs_sorted[!missing]
# > missing_horvath_cgs
# [1] "cg02654291" "cg02972551" "cg04431054" "cg05590257" "cg06117855" "cg09785172" "cg09869858" "cg13682722" "cg14329157"
# [10] "cg16494477" "cg17408647" "cg19046959" "cg19167673" "cg19273182" "cg19569684" "cg19945840" "cg24471894" "cg27016307"
# [19] "cg27319898"

horvath_beta_row_names <- row_names_beta_sorted %in% horvath_cgs_sorted
horvath_beta <- beta[row_names_beta_sorted %in% horvath_cgs_sorted,]
write.csv(horvath_beta, "horvath_beta.csv")

levine_cgs_sorted <- str_sort(levine_cgs)
missing_places_levine <- levine_cgs_sorted %in% row_names_beta_sorted
missing_levine_cgs <- levine_cgs_sorted[!missing_places_levine]

levine_beta_row_names <- row_names_beta_sorted %in% levine_cgs_sorted
levine_beta <- beta[row_names_beta_sorted %in% levine_cgs_sorted,]
write.csv(levine_beta, "levine_beta.csv")
# cpg.names <- getInputCpgNames(myDNAmAge_with_acceleration_age)
# if (!all(c("coefHorvath", "coefHannum", "coefLevine", "coefSkin", 
#            "coefPedBE", "coefWu", "coefTL") %in% ls(.GlobalEnv))) {
#   load_DNAm_Clocks_data()
# }
# checkHorvath <- coefHorvath$CpGmarker[-1][!coefHorvath$CpGmarker[-1] %in% cpg.names]
# checkHannum <- coefHannum$CpGmarker[!coefHannum$CpGmarker %in% 
#                                       cpg.names]
# checkLevine <- coefLevine$CpGmarker[-1][!coefLevine$CpGmarker[-1] %in% 
#                                           cpg.names]
# checkHorvath <- coefHorvath$CpGmarker[-1][!coefHorvath$CpGmarker[-1] %in% cpg.names]
###

# bind metadata with the DNAmAge
myDNAmAge_with_acceleration_age_with_metadata <- cbind(myDNAmAge_with_acceleration_age, metadata_file)

write.csv(myDNAmAge_with_acceleration_age_with_metadata, "csv/myDNAmAge_with_acceleration_age_with_metadata.csv")

# remove unneeded columns
drop_cols <- c("PedBE", "Wu", "TL", "BNN")
myDNAmAge_with_acceleration_with_metadata_without_pedbe_wu_tl_bnn <- myDNAmAge_with_acceleration_age_with_metadata[ , !(names(myDNAmAge_with_acceleration_age_with_metadata) %in% drop_cols)]
write.csv(myDNAmAge_with_acceleration_with_metadata_without_pedbe_wu_tl_bnn, "myDNAmAge_with_acceleration_with_metadata_without_pedbe_wu_tl_bnn.csv")


#plot 
# plotDNAmAge(myDNAmAge_with_acceleration_with_metadata_without_pedbe_wu_tl_bnn$Horvath, metadata_file$Age)
# plotDNAmAge(myDNAmAge_with_acceleration_with_metadata_without_pedbe_wu_tl_bnn$Hannum, metadata_file$Age, tit="Hannum")
# plotDNAmAge(myDNAmAge_with_acceleration_with_metadata_without_pedbe_wu_tl_bnn$Levine, metadata_file$Age, tit="Levine")
# plotDNAmAge(myDNAmAge_with_acceleration_with_metadata_without_pedbe_wu_tl_bnn$skinHorvath, metadata_file$Age, tit="skinHorvath")


# make for each clock by sample group
myData <- myDNAmAge_with_acceleration_with_metadata_without_pedbe_wu_tl_bnn
filter_sporadic <- subset(myData, Sample_Group == "Sporadic")
filter_MEN1 <- subset(myData, Sample_Group == "MEN1")
filter_VHL <- subset(myData, Sample_Group == "VHL")
rownames(filter_sporadic) <- NULL

# plotDNAmAge(filter_sporadic$Horvath, filter_sporadic$Age, tit="Hovarth (sporadic)")
# plotDNAmAge(filter_MEN1$Horvath, filter_MEN1$Age, tit="Hovarth (MEN1)")
# plotDNAmAge(filter_VHL$Horvath, filter_VHL$Age, tit="Hovarth (VHL)")

# plotDNAmAge(filter_sporadic$Hannum, filter_sporadic$Age, tit="Hannum (sporadic)")
# plotDNAmAge(filter_MEN1$Hannum, filter_MEN1$Age, tit="Hannum (MEN1)")
# plotDNAmAge(filter_VHL$Hannum, filter_VHL$Age, tit="Hannum (VHL)")

# plotDNAmAge(filter_sporadic$Levine, filter_sporadic$Age, tit="Levine (sporadic)")
# plotDNAmAge(filter_MEN1$Levine, filter_MEN1$Age, tit="Levine (MEN1)")
# plotDNAmAge(filter_VHL$Levine, filter_VHL$Age, tit="Levine (VHL)")

# plotDNAmAge(filter_sporadic$skinHorvath, filter_sporadic$Age, tit="skinHorvath (sporadic)")
# plotDNAmAge(filter_MEN1$skinHorvath, filter_MEN1$Age, tit="skinHorvath (MEN1)")
# plotDNAmAge(filter_VHL$skinHorvath, filter_VHL$Age, tit="skinHorvath (VHL)")

# for each clock calculate
# 1. mean deviation 2. average and present the value for AnyMet (yes/no)
# filter_AnyMetYes <- subset(myData, AnyMet == "Yes")
# filter_AnyMetNo <- subset(myData, AnyMet == "No")

is_filter_negative_calculated_age = TRUE
# out_file_avg = "any_met_avg_yes_no_by_clocks.csv"
# out_file_sd = "any_met_sd_yes_no_by_clocks.csv"

# if (is_filter_negative_calculated_age) {
#   filter_AnyMetYes <- subset(filter_AnyMetYes, Levine > 0)
#   filter_AnyMetNo <- subset(filter_AnyMetNo, Levine > 0)
#   out_file_avg = "any_met_avg_yes_no_by_clocks_without_negative_levine.csv"
#   out_file_sd = "any_met_sd_yes_no_by_clocks_without_negative_levine.csv"
# }

# avg_by_clock <- c()
# sd_by_clock <- c()
# for (clock_type in c("ageAcc.Horvath", "ageAcc2.Horvath", "ageAcc3.Horvath", "ageAcc.Hannum", "ageAcc2.Hannum", "ageAcc3.Hannum", "ageAcc.Levine", "ageAcc2.Levine", "ageAcc3.Levine", "ageAcc.Hovarth2", "ageAcc2.Hovarth2", "ageAcc3.Hovarth2")) {
#     avg_yes <- mean(filter_AnyMetYes[[clock_type]], na.rm=TRUE)
#     avg_no <- mean(filter_AnyMetNo[[clock_type]], na.rm=TRUE)
#     avg_by_clock <- rbind(avg_by_clock, c(clock_type, avg_yes, avg_no))
#     
#     sd_yes <- sd(filter_AnyMetYes[[clock_type]], na.rm=TRUE)
#     sd_no <- sd(filter_AnyMetNo[[clock_type]], na.rm=TRUE)
#     sd_by_clock <- rbind(sd_by_clock, c(clock_type, sd_yes, sd_no))
# }
# colnames(avg_by_clock) <- c("clock_type", "Yes_avg", "No_avg")
# write.csv(avg_by_clock, out_file_avg)
# 
# colnames(sd_by_clock) <- c("clock_type", "Yes_sd", "No_sd")
# write.csv(sd_by_clock, out_file_sd)

# To perform a TTEST for acceleration1 for Hovarth and Levine
# And check if there is a difference between Sporadic and Vhl
# using https://www.youtube.com/watch?v=RlhnNbPZC0A&ab_channel=MarinStatsLectures-RProgramming%26Statistics
all_data = myDNAmAge_with_acceleration_age_with_metadata
# boxplot(myDNAmAge_with_acceleration_age_with_metadata$ageAcc.Levine ~ myDNAmAge_with_acceleration_age_with_metadata$Sample_Group)
# boxplot(myDNAmAge_with_acceleration_age_with_metadata$ageAcc.Horvath ~ myDNAmAge_with_acceleration_age_with_metadata$Sample_Group)

t.test(all_data$Horvath[all_data$Sample_Group == "Sporadic"], all_data$Horvath[all_data$Sample_Group == "VHL"])

# data:  all_data$Horvath[all_data$Sample_Group == "Sporadic"] and all_data$Horvath[all_data$Sample_Group == "VHL"]
# t = 0.81329, df = 15.477, p-value = 0.4284
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -7.247031 16.228669
# sample estimates:
#   mean of x mean of y 
# 63.66140  59.17058 

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

################ LEVINE ##############

t.test(all_data$Levine[all_data$Sample_Group == "Sporadic"], all_data$Levine[all_data$Sample_Group == "VHL"])

# data:  all_data$Levine[all_data$Sample_Group == "Sporadic"] and all_data$Levine[all_data$Sample_Group == "VHL"]
# t = 0.97601, df = 21.57, p-value = 0.3399
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -5.323502 14.768248
# sample estimates:
#   mean of x mean of y 
# 38.93489  34.21252 

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

t.test(all_data$ageAcc3.Levine[all_data$Sample_Group == "Sporadic"], all_data$ageAcc3.Levine[all_data$Sample_Group == "VHL"])

# data:  all_data$ageAcc3.Levine[all_data$Sample_Group == "Sporadic"] and all_data$ageAcc3.Levine[all_data$Sample_Group == "VHL"]
# t = 2.352, df = 32.698, p-value = 0.02485
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   1.143115 15.833131
# sample estimates:
#   mean of x mean of y 
# 5.073437 -3.414686 

########## LEVINE NO NEGATIVE ###############

no_negative_levine <- read.csv("csv/myDNAmAge_with_acceleration_with_metadata_without_pedbe_wu_tl_bnn_no_negative_levine.csv")
t.test(no_negative_levine$ageAcc.Levine[no_negative_levine$Sample_Group == "Sporadic"], no_negative_levine$ageAcc.Levine[no_negative_levine$Sample_Group == "VHL"])

# Welch Two Sample t-test

# data:  no_negative_levine$ageAcc.Levine[no_negative_levine$Sample_Group == "Sporadic"] and no_negative_levine$ageAcc.Levine[no_negative_levine$Sample_Group == "VHL"]
# t = -0.66194, df = 38.345, p-value = 0.512
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -11.260729   5.709995
# sample estimates:
#   mean of x mean of y 
# -14.56285 -11.78748 

t.test(no_negative_levine$ageAcc2.Levine[no_negative_levine$Sample_Group == "Sporadic"], no_negative_levine$ageAcc2.Levine[no_negative_levine$Sample_Group == "VHL"])

# Welch Two Sample t-test
# 
# data:  no_negative_levine$ageAcc2.Levine[no_negative_levine$Sample_Group == "Sporadic"] and no_negative_levine$ageAcc2.Levine[no_negative_levine$Sample_Group == "VHL"]
# t = 0.49885, df = 26.175, p-value = 0.6221
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -6.318595 10.370033
# sample estimates:
#   mean of x mean of y 
# 5.655887  3.630168 

t.test(no_negative_levine$ageAcc3.Levine[no_negative_levine$Sample_Group == "Sporadic"], no_negative_levine$ageAcc3.Levine[no_negative_levine$Sample_Group == "VHL"])

# data:  no_negative_levine$ageAcc3.Levine[no_negative_levine$Sample_Group == "Sporadic"] and no_negative_levine$ageAcc3.Levine[no_negative_levine$Sample_Group == "VHL"]
# t = 2.7853, df = 31.115, p-value = 0.009023
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   2.648601 17.126450
# sample estimates:
#   mean of x mean of y 
# 6.472839 -3.414686 


############ WU #####

t.test(all_data$Wu[all_data$Sample_Group == "Sporadic"], all_data$Wu[all_data$Sample_Group == "VHL"])

# data:  all_data$Wu[all_data$Sample_Group == "Sporadic"] and all_data$Wu[all_data$Sample_Group == "VHL"]
# t = -1.2821, df = 15.569, p-value = 0.2186
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -3.1842357  0.7875924
# sample estimates:
#   mean of x mean of y 
# 5.974587  7.172909 

t.test(all_data$ageAcc.Wu[all_data$Sample_Group == "Sporadic"], all_data$ageAcc.Wu[all_data$Sample_Group == "VHL"])

# data:  all_data$ageAcc.Wu[all_data$Sample_Group == "Sporadic"] and all_data$ageAcc.Wu[all_data$Sample_Group == "VHL"]
# t = -3.3153, df = 35.516, p-value = 0.002117
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -15.963663  -3.842071
# sample estimates:
#   mean of x mean of y 
# -48.72996 -38.82709 

t.test(all_data$ageAcc2.Wu[all_data$Sample_Group == "Sporadic"], all_data$ageAcc2.Wu[all_data$Sample_Group == "VHL"])

# data:  all_data$ageAcc2.Wu[all_data$Sample_Group == "Sporadic"] and all_data$ageAcc2.Wu[all_data$Sample_Group == "VHL"]
# t = -1.9217, df = 17.287, p-value = 0.07129
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -3.4122120  0.1570926
# sample estimates:
#   mean of x  mean of y 
# -0.6333620  0.9941978 

t.test(all_data$ageAcc3.Wu[all_data$Sample_Group == "Sporadic"], all_data$ageAcc3.Wu[all_data$Sample_Group == "VHL"])

# data:  all_data$ageAcc3.Wu[all_data$Sample_Group == "Sporadic"] and all_data$ageAcc3.Wu[all_data$Sample_Group == "VHL"]
# t = -1.2944, df = 18.689, p-value = 0.2113
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -2.7222215  0.6432136
# sample estimates:
#   mean of x  mean of y 
# -0.3746660  0.6648379 

# TODO metastasis predicition accel vs features
# TODO does it divide to sample groups unsupervised
# TODO VHL+MEN1 vs Sporadic
# TODO in MEN1 check for only under the age of 50

# https://www.youtube.com/watch?v=BCJ7gLQ9MZM&ab_channel=CrackEconomicsandStatistics
# test age for MAN1, Sporadic, VHL groups for normal distribution
curr_data <- myDNAmAge_with_acceleration_age_with_metadata
horvath_age_sporadic <- curr_data$Horvath[curr_data$Sample_Group=='Sporadic']
shapiro.test(horvath_age_sporadic)
horvath_age_vhl <- curr_data$Horvath[curr_data$Sample_Group=='VHL']
shapiro.test(horvath_age_vhl)
vhl <- 'VHL'
# data:  horvath_age_vhl
# W = 0.94214, p-value = 0.5771
shapiro.test(curr_data$ageAcc.Horvath[curr_data$Sample_Group==vhl])
# data:  curr_data$ageAcc.Horvath[curr_data$Sample_Group == "VHL"]
# W = 0.87703, p-value = 0.1206
shapiro.test(curr_data$ageAcc2.Horvath[curr_data$Sample_Group=='VHL'])
# data:  curr_data$ageAcc2.Horvath[curr_data$Sample_Group == "VHL"]
# W = 0.93607, p-value = 0.5101
shapiro.test(curr_data$ageAcc3.Horvath[curr_data$Sample_Group=='VHL'])
# data:  curr_data$ageAcc3.Horvath[curr_data$Sample_Group == "VHL"]
# W = 0.94136, p-value = 0.5682
shapiro.test(curr_data$Levine[curr_data$Sample_Group=='VHL'])
# data:  curr_data$Levine[curr_data$Sample_Group == "VHL"]
# W = 0.94839, p-value = 0.6495
shapiro.test(curr_data$ageAcc.Levine[curr_data$Sample_Group=='VHL'])
# data:  curr_data$ageAcc.Levine[curr_data$Sample_Group == "VHL"]
# W = 0.90279, p-value = 0.235
shapiro.test(curr_data$ageAcc2.Levine[curr_data$Sample_Group=='VHL'])
# data:  curr_data$ageAcc2.Levine[curr_data$Sample_Group == "VHL"]
# W = 0.95345, p-value = 0.7094
shapiro.test(curr_data$ageAcc3.Levine[curr_data$Sample_Group=='VHL'])
# data:  curr_data$ageAcc3.Levine[curr_data$Sample_Group == "VHL"]
# W = 0.85219, p-value = 0.06168

sporadic <- 'Sporadic'
shapiro.test(curr_data$Horvath[curr_data$Sample_Group==sporadic])
shapiro.test(curr_data$ageAcc.Horvath[curr_data$Sample_Group==sporadic])
shapiro.test(curr_data$ageAcc2.Horvath[curr_data$Sample_Group==sporadic])
shapiro.test(curr_data$ageAcc3.Horvath[curr_data$Sample_Group==sporadic])
shapiro.test(curr_data$Levine[curr_data$Sample_Group==sporadic])
shapiro.test(curr_data$ageAcc.Levine[curr_data$Sample_Group==sporadic])
shapiro.test(curr_data$ageAcc2.Levine[curr_data$Sample_Group==sporadic])
shapiro.test(curr_data$ageAcc3.Levine[curr_data$Sample_Group==sporadic])

# > shapiro.test(curr_data$Horvath[curr_data$Sample_Group==sporadic])
# Shapiro-Wilk normality test
# data:  curr_data$Horvath[curr_data$Sample_Group == sporadic]
# W = 0.95538, p-value = 0.08736
# 
# > shapiro.test(curr_data$ageAcc.Horvath[curr_data$Sample_Group==sporadic])
# Shapiro-Wilk normality test
# data:  curr_data$ageAcc.Horvath[curr_data$Sample_Group == sporadic]
# W = 0.97847, p-value = 0.5736
# 
# > shapiro.test(curr_data$ageAcc2.Horvath[curr_data$Sample_Group==sporadic])
# Shapiro-Wilk normality test
# data:  curr_data$ageAcc2.Horvath[curr_data$Sample_Group == sporadic]
# W = 0.96049, p-value = 0.1357
# 
# > shapiro.test(curr_data$ageAcc3.Horvath[curr_data$Sample_Group==sporadic])
# Shapiro-Wilk normality tes
# data:  curr_data$ageAcc3.Horvath[curr_data$Sample_Group == sporadic]
# W = 0.95599, p-value = 0.09213
# 
# > shapiro.test(curr_data$Levine[curr_data$Sample_Group==sporadic])
# Shapiro-Wilk normality test
# data:  curr_data$Levine[curr_data$Sample_Group == sporadic]
# W = 0.9706, p-value = 0.3177
# 
# > shapiro.test(curr_data$ageAcc.Levine[curr_data$Sample_Group==sporadic])
# Shapiro-Wilk normality test
# data:  curr_data$ageAcc.Levine[curr_data$Sample_Group == sporadic]
# W = 0.94559, p-value = 0.03781
# 
# > shapiro.test(curr_data$ageAcc2.Levine[curr_data$Sample_Group==sporadic])
# Shapiro-Wilk normality test
# data:  curr_data$ageAcc2.Levine[curr_data$Sample_Group == sporadic]
# W = 0.96554, p-value = 0.2089
# 
# > shapiro.test(curr_data$ageAcc3.Levine[curr_data$Sample_Group==sporadic])
# Shapiro-Wilk normality test
# data:  curr_data$ageAcc3.Levine[curr_data$Sample_Group == sporadic]
# W = 0.98108, p-value = 0.677

men1 <- 'MEN1'
shapiro.test(curr_data$Horvath[curr_data$Sample_Group==men1])
shapiro.test(curr_data$ageAcc.Horvath[curr_data$Sample_Group==men1])
shapiro.test(curr_data$ageAcc2.Horvath[curr_data$Sample_Group==men1])
shapiro.test(curr_data$ageAcc3.Horvath[curr_data$Sample_Group==men1])
shapiro.test(curr_data$Levine[curr_data$Sample_Group==men1])
shapiro.test(curr_data$ageAcc.Levine[curr_data$Sample_Group==men1])
shapiro.test(curr_data$ageAcc2.Levine[curr_data$Sample_Group==men1])
shapiro.test(curr_data$ageAcc3.Levine[curr_data$Sample_Group==men1])

# > shapiro.test(curr_data$Horvath[curr_data$Sample_Group==men1])
# Shapiro-Wilk normality test
# data:  curr_data$Horvath[curr_data$Sample_Group == men1]
# W = 0.91676, p-value = 0.004756
# 
# > shapiro.test(curr_data$ageAcc.Horvath[curr_data$Sample_Group==men1])
# Shapiro-Wilk normality test
# data:  curr_data$ageAcc.Horvath[curr_data$Sample_Group == men1]
# W = 0.94279, p-value = 0.03581
# 
# > shapiro.test(curr_data$ageAcc2.Horvath[curr_data$Sample_Group==men1])
# Shapiro-Wilk normality test
# data:  curr_data$ageAcc2.Horvath[curr_data$Sample_Group == men1]
# W = 0.89891, p-value = 0.001336
# 
# > shapiro.test(curr_data$ageAcc3.Horvath[curr_data$Sample_Group==men1])
# Shapiro-Wilk normality tes
# data:  curr_data$ageAcc3.Horvath[curr_data$Sample_Group == men1]
# W = 0.88053, p-value = 0.0003955
# 
# > shapiro.test(curr_data$Levine[curr_data$Sample_Group==men1])
# Shapiro-Wilk normality tes
# data:  curr_data$Levine[curr_data$Sample_Group == men1]
# W = 0.93218, p-value = 0.01535
# 
# > shapiro.test(curr_data$ageAcc.Levine[curr_data$Sample_Group==men1])
# Shapiro-Wilk normality test
# data:  curr_data$ageAcc.Levine[curr_data$Sample_Group == men1]
# W = 0.9552, p-value = 0.09931
# 
# > shapiro.test(curr_data$ageAcc2.Levine[curr_data$Sample_Group==men1])
# Shapiro-Wilk normality test
# data:  curr_data$ageAcc2.Levine[curr_data$Sample_Group == men1]
# W = 0.92167, p-value = 0.006854
# 
# > shapiro.test(curr_data$ageAcc3.Levine[curr_data$Sample_Group==men1])
# Shapiro-Wilk normality test
# data:  curr_data$ageAcc3.Levine[curr_data$Sample_Group == men1]
# W = 0.93548, p-value = 0.01992

length(curr_data$Horvath[curr_data$Sample_Group==men1])

