# random forest
# https://www.youtube.com/watch?v=6EXPYzbfLCE&t=70s

install.packages("stats")
install.packages("dplyr")
install.packages("randomForest")
install.packages("doSNOW")

library(stats)
library(dplyr)
library(randomForest)

# Normalize the values to be between 0 and 1
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x)) # can it be that max - min = 0? then we have division by 0
}

baseDir <- getwd()
sample_sheet_filename = "/Sample_Sheet_Full.csv"
path <- paste(baseDir, sample_sheet_filename, sep="")
dataOnly = read.csv(path, skip = 8, header = F)
headers = read.csv(path, skip = 7, header = F, nrows = 1, as.is = T)
colnames(dataOnly) = headers
metadata_file <- dataOnly
# metadata_file_pnet_only <- metadata_file[metadata_file$Site == 'PNET', ]

metastasis <- as.factor(metadata_file$AnyMet)

csv_dir <- file.path(getwd(), "")
beta_csv_path <- file.path(csv_dir, "beta.csv") 

##### change from here

beta <- read.csv(file=beta_csv_path, row.names = 1) # nrows=100000
beta_without_na <- na.omit(beta) # remove cgX which have NA values since KNN needs to have full data

# transpose so we have a row of a single patient instead of a column
data_frame <- data.frame(t(beta_without_na))

# normalize to have distributed values between 0 and 1
data_frame_norm <- as.data.frame(lapply(data_frame, min_max_norm))
which(is.na(data_frame_norm)) # make sure there are no NA

model <- randomForest(data_frame_norm, metastasis, proximity = TRUE)

### for 100,000
# Call:
#   randomForest(x = data_frame_norm, y = metastasis, proximity = TRUE) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 316
# 
# OOB estimate of  error rate: 35.42%
# Confusion matrix:
#     No Yes class.error
# No   2  26   0.9285714
# Yes  8  60   0.1176471
# 

### for all
# Call:
#   randomForest(x = data_frame_norm, y = metastasis, proximity = TRUE) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 929
# 
# OOB estimate of  error rate: 34.38%
# Confusion matrix:
#     No Yes class.error
# No   1  27  0.96428571
# Yes  6  62  0.08823529

model_1000 <- randomForest(data_frame_norm, metastasis, ntree=1000, proximity = TRUE)

# 
# Call:
#   randomForest(x = data_frame_norm, y = metastasis, ntree = 1000,      proximity = TRUE) 
# Type of random forest: classification
# Number of trees: 1000
# No. of variables tried at each split: 929
# 
# OOB estimate of  error rate: 34.38%
# Confusion matrix:
#      No Yes class.error
# No   0  28  1.00000000
# Yes  5  63  0.07352941

oob.error.data <- data.frame(
  Trees=rep(1:nrow(model_1000$err.rate), times=3),
  Type=rep(c("OOB", "Yes", "No"), each=nrow(model_1000$err.rate)),
  Error=c(model_1000$err.rate[,"OOB"],
          model_1000$err.rate[,"Yes"],
          model_1000$err.rate[,"No"]))

obb_error_data_125 <- oob.error.data[oob.error.data$Trees < 125,]

library(ggplot2)
ggplot(data=oob.error.data, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))
ggplot(data=obb_error_data_125, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))
oob.error.data[oob.error.data["Trees"]==50,]

# > oob.error.data[oob.error.data["Trees"]==50,] - 50 trees is the optimum
# Trees Type     Error
# 50      50  OOB 0.3020833
# 1050    50  Yes 0.1323529
# 2050    50   No 0.7142857

##########################

model_500_100 <- randomForest(data_frame_norm, metastasis, ntree=500, mtry=100, proximity = TRUE)

oob.error.data_500_100 <- data.frame(
  Trees=rep(1:nrow(model_500_100$err.rate), times=3),
  Type=rep(c("OOB", "Yes", "No"), each=nrow(model_500_100$err.rate)),
  Error=c(model_500_100$err.rate[,"OOB"],
          model_500_100$err.rate[,"Yes"],
          model_500_100$err.rate[,"No"]))

ggplot(data=oob.error.data_500_100, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))

obb_error_data_50_100 <- oob.error.data[oob.error.data$Trees < 50,]
ggplot(data=obb_error_data_50_100, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))


# Call:
#   randomForest(x = data_frame_norm, y = metastasis, ntree = 50,      mtry = 100, proximity = TRUE) 
# Type of random forest: classification
# Number of trees: 50
# No. of variables tried at each split: 100
# 
# OOB estimate of  error rate: 37.5%
# Confusion matrix:
#   No Yes class.error
# No   3  25   0.8928571
# Yes 11  57   0.1617647

# Call:
#   randomForest(x = data_frame_norm, y = metastasis, ntree = 500,      mtry = 100, proximity = TRUE) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 100
# 
# OOB estimate of  error rate: 33.33%
# Confusion matrix:
#   No Yes class.error
# No   1  27  0.96428571
# Yes  5  63  0.07352941

# 18 is the optimum
obb_error_data_50_100[obb_error_data_50_100["Trees"]==18,]

# Trees Type     Error
# 18      18  OOB 0.3125000
# 1018    18  Yes 0.1323529
# 2018    18   No 0.7500000

##########################

model_50_3 <- randomForest(data_frame_norm, metastasis, ntree=50, mtry=3, proximity = TRUE)

oob.error.data_50_3 <- data.frame(
  Trees=rep(1:nrow(model_50_3$err.rate), times=3),
  Type=rep(c("OOB", "Yes", "No"), each=nrow(model_50_3$err.rate)),
  Error=c(model_50_3$err.rate[,"OOB"],
          model_50_3$err.rate[,"Yes"],
          model_50_3$err.rate[,"No"]))

ggplot(data=oob.error.data_50_3, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))
oob.error.data_50_3[oob.error.data_50_3["Trees"]==43,]

########################## HORVATH RANDOM FOREST

horvath_levine_csv_dir <- file.path(getwd(), "csv")
horvath_beta_csv_path <- file.path(horvath_levine_csv_dir, "horvath_beta.csv") 
horvath_beta <- read.csv(file=horvath_beta_csv_path, row.names = 1) 
horvath_beta_without_na <- na.omit(horvath_beta) # remove cgX which have NA values since KNN needs to have full data

# transpose so we have a row of a single patient instead of a column
horvath_data_frame <- data.frame(t(horvath_beta_without_na))

# normalize to have distributed values between 0 and 1
horvath_data_frame_norm <- as.data.frame(lapply(horvath_data_frame, min_max_norm))
which(is.na(horvath_data_frame_norm)) # make sure there are no NA

horvath_random_forest <- randomForest(horvath_data_frame_norm, metastasis, proximity = TRUE)

oob.error.data_horvath <- data.frame(
  Trees=rep(1:nrow(horvath_random_forest$err.rate), times=3),
  Type=rep(c("OOB", "Yes", "No"), each=nrow(horvath_random_forest$err.rate)),
  Error=c(horvath_random_forest$err.rate[,"OOB"],
          horvath_random_forest$err.rate[,"Yes"],
          horvath_random_forest$err.rate[,"No"]))

ggplot(data=oob.error.data_horvath, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))

oob.error.data_horvath_50 <- oob.error.data_horvath[oob.error.data_horvath$Trees < 50,]
ggplot(data=oob.error.data_horvath_50, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))

horvath_random_forest_38 <- randomForest(horvath_data_frame_norm, ntree=38, metastasis, proximity = TRUE)
horvath_random_forest_38

# > horvath_random_forest_38
# 
# Call:
#   randomForest(x = horvath_data_frame_norm, y = metastasis, ntree = 38,      proximity = TRUE) 
# Type of random forest: classification
# Number of trees: 38
# No. of variables tried at each split: 18
# 
# OOB estimate of  error rate: 29.17%
# Confusion matrix:
#     No Yes  class.error
# No   5  23  0.82142857
# Yes  5  63  0.07352941

########################## LEVINE RANDOM FOREST

horvath_levine_csv_dir <- file.path(getwd(), "csv")
levine_beta_csv_path <- file.path(horvath_levine_csv_dir, "levine_beta.csv") 
levine_beta <- read.csv(file=levine_beta_csv_path, row.names = 1) 
levine_beta_without_na <- na.omit(levine_beta) # remove cgX which have NA values since KNN needs to have full data

# transpose so we have a row of a single patient instead of a column
levine_data_frame <- data.frame(t(levine_beta_without_na))

# normalize to have distributed values between 0 and 1
levine_data_frame_norm <- as.data.frame(lapply(levine_data_frame, min_max_norm))
which(is.na(horvath_data_frame_norm)) # make sure there are no NA

set.seed(43) 
levine_random_forest <- randomForest(levine_data_frame_norm, metastasis, proximity = TRUE)

oob.error.data_levine <- data.frame(
  Trees=rep(1:nrow(levine_random_forest$err.rate), times=3),
  Type=rep(c("OOB", "Yes", "No"), each=nrow(levine_random_forest$err.rate)),
  Error=c(levine_random_forest$err.rate[,"OOB"],
          levine_random_forest$err.rate[,"Yes"],
          levine_random_forest$err.rate[,"No"]))

ggplot(data=oob.error.data_levine, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))

# oob.error.data_horvath_50 <- oob.error.data_horvath[oob.error.data_horvath$Trees < 50,]
# ggplot(data=oob.error.data_horvath_50, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))
# 
# horvath_random_forest_38 <- randomForest(horvath_data_frame_norm, ntree=38, metastasis, proximity = TRUE)
# horvath_random_forest_38

# levine_random_forest$err.rate 
#  [23,] 0.2812500 0.8214286 0.05882353

levine_random_forest_23 <- randomForest(levine_data_frame_norm, metastasis, ntree=23, proximity = TRUE)
levine_random_forest_23

# > levine_random_forest_23
# 
# Call:
#   randomForest(x = levine_data_frame_norm, y = metastasis, ntree = 23,      proximity = TRUE) 
# Type of random forest: classification
# Number of trees: 23
# No. of variables tried at each split: 22
# 
# OOB estimate of  error rate: 30.21%
# Confusion matrix:
#       No Yes class.error
# No   6  22   0.7857143
# Yes  7  61   0.1029412

########################## LEVINE SPORADIC ONLY

levine_data_frame_norm$Sample_Group <- metadata_file$Sample_Group
levine_data_frame_norm$AnyMet <- as.factor(metadata_file$AnyMet)
levine_sporadic <- levine_data_frame_norm[levine_data_frame_norm$Sample_Group == 'Sporadic',]
levine_sporadic_rf <- randomForest(levine_sporadic, levine_sporadic$AnyMet, ntree=17, proximity = TRUE)

oob.error.levine_sporadic_rf <- data.frame(
  Trees=rep(1:nrow(levine_sporadic_rf$err.rate), times=3),
  Type=rep(c("OOB", "Yes", "No"), each=nrow(levine_sporadic_rf$err.rate)),
  Error=c(levine_sporadic_rf$err.rate[,"OOB"],
          levine_sporadic_rf$err.rate[,"Yes"],
          levine_sporadic_rf$err.rate[,"No"]))

ggplot(data=oob.error.levine_sporadic_rf , aes(x=Trees, y=Error)) + geom_line(aes(color=Type))

# Call:
#   randomForest(x = levine_sporadic, y = levine_sporadic$AnyMet,      ntree = 17, proximity = TRUE) 
# Type of random forest: classification
# Number of trees: 17
# No. of variables tried at each split: 22
# 
# OOB estimate of  error rate: 11.36%
# Confusion matrix:
#     No Yes class.error
# No   9   4  0.30769231
# Yes  1  30  0.03225806

########################## HORVATH SPORADIC ONLY

horvath_data_frame_norm$Sample_Group <- metadata_file$Sample_Group
horvath_data_frame_norm$AnyMet <- as.factor(metadata_file$AnyMet)
horvath_sporadic <- horvath_data_frame_norm[horvath_data_frame_norm$Sample_Group == 'Sporadic',]
horvath_sporadic_rf <- randomForest(horvath_sporadic, horvath_sporadic$AnyMet, ntree=17, proximity = TRUE)

oob.error.horvath_sporadic_rf <- data.frame(
  Trees=rep(1:nrow(horvath_sporadic_rf$err.rate), times=3),
  Type=rep(c("OOB", "Yes", "No"), each=nrow(horvath_sporadic_rf$err.rate)),
  Error=c(horvath_sporadic_rf$err.rate[,"OOB"],
          horvath_sporadic_rf$err.rate[,"Yes"],
          horvath_sporadic_rf$err.rate[,"No"]))

ggplot(data=oob.error.horvath_sporadic_rf , aes(x=Trees, y=Error)) + geom_line(aes(color=Type))

# Call:
#   randomForest(x = horvath_sporadic, y = horvath_sporadic$AnyMet,      ntree = 17, proximity = TRUE) 
# Type of random forest: classification
# Number of trees: 17
# No. of variables tried at each split: 18
# 
# OOB estimate of  error rate: 22.73%
# Confusion matrix:
#   No Yes class.error
# No   7   6   0.4615385
# Yes  4  27   0.1290323

########################## HORVATH meN1 ONLY

# horvath_data_frame_norm$Sample_Group <- metadata_file$Sample_Group
# horvath_data_frame_norm$AnyMet <- as.factor(metadata_file$AnyMet)
horvath_men1 <- horvath_data_frame_norm[horvath_data_frame_norm$Sample_Group == 'MEN1',]
horvath_men1_rf <- randomForest(horvath_men1, horvath_men1$AnyMet, ntree=15, proximity = TRUE)

oob.error.horvath_men1_rf <- data.frame(
  Trees=rep(1:nrow(horvath_men1_rf$err.rate), times=3),
  Type=rep(c("OOB", "Yes", "No"), each=nrow(horvath_men1_rf$err.rate)),
  Error=c(horvath_men1_rf$err.rate[,"OOB"],
          horvath_men1_rf$err.rate[,"Yes"],
          horvath_men1_rf$err.rate[,"No"]))

ggplot(data=oob.error.horvath_men1_rf, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))

# Call:
#   randomForest(x = horvath_men1, y = horvath_men1$AnyMet, ntree = 15,      proximity = TRUE) 
# Type of random forest: classification
# Number of trees: 15
# No. of variables tried at each split: 18
# 
# OOB estimate of  error rate: 26.19%
# Confusion matrix:
#   No Yes class.error
# No   2   7   0.7777778
# Yes  4  29   0.1212121

########################## HORVATH VHL ONLY

horvath_vhl <- horvath_data_frame_norm[horvath_data_frame_norm$Sample_Group == 'VHL',]
horvath_vhl_rf <- randomForest(horvath_vhl, horvath_vhl$AnyMet, ntree=8, proximity = TRUE)

oob.error.horvath_vhl_rf <- data.frame(
  Trees=rep(1:nrow(horvath_vhl_rf$err.rate), times=3),
  Type=rep(c("OOB", "Yes", "No"), each=nrow(horvath_vhl_rf$err.rate)),
  Error=c(horvath_vhl_rf$err.rate[,"OOB"],
          horvath_vhl_rf$err.rate[,"Yes"],
          horvath_vhl_rf$err.rate[,"No"]))

ggplot(data=oob.error.horvath_vhl_rf, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))

# Call:
#   randomForest(x = horvath_vhl, y = horvath_vhl$AnyMet, ntree = 8,      proximity = TRUE) 
# Type of random forest: classification
# Number of trees: 8
# No. of variables tried at each split: 18
# 
# OOB estimate of  error rate: 30%
# Confusion matrix:
#   No Yes class.error
# No   4   2   0.3333333
# Yes  1   3   0.2500000

########################## HORVATH SPORADIC + MEN1 ONLY

horvath_sporadic_men1 <- horvath_data_frame_norm[horvath_data_frame_norm$Sample_Group != 'VHL',]
horvath_sporadic_men1_rf <- randomForest(horvath_sporadic_men1, horvath_sporadic_men1$AnyMet, ntree=16, proximity = TRUE)

oob.error.horvath_sporadic_men1_rf <- data.frame(
  Trees=rep(1:nrow(horvath_sporadic_men1_rf$err.rate), times=3),
  Type=rep(c("OOB", "Yes", "No"), each=nrow(horvath_sporadic_men1_rf$err.rate)),
  Error=c(horvath_sporadic_men1_rf$err.rate[,"OOB"],
          horvath_sporadic_men1_rf$err.rate[,"Yes"],
          horvath_sporadic_men1_rf$err.rate[,"No"]))

ggplot(data=oob.error.horvath_sporadic_men1_rf, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))

# Call:
#   randomForest(x = horvath_sporadic_men1, y = horvath_sporadic_men1$AnyMet,      ntree = 16, proximity = TRUE) 
# Type of random forest: classification
# Number of trees: 16
# No. of variables tried at each split: 18
# 
# OOB estimate of  error rate: 22.09%
# Confusion matrix:
#   No Yes class.error
# No  12  10   0.4545455
# Yes  9  55   0.1406250
####################

########################## HORVATH SPORADIC + VHL ONLY

horvath_sporadic_vhl <- horvath_data_frame_norm[horvath_data_frame_norm$Sample_Group != 'Sporadic',]
horvath_sporadic_vhl_rf <- randomForest(horvath_sporadic_vhl, horvath_sporadic_vhl$AnyMet, ntree=28, proximity = TRUE)

oob.error.horvath_sporadic_vhl_rf <- data.frame(
  Trees=rep(1:nrow(horvath_sporadic_vhl_rf$err.rate), times=3),
  Type=rep(c("OOB", "Yes", "No"), each=nrow(horvath_sporadic_vhl_rf$err.rate)),
  Error=c(horvath_sporadic_vhl_rf$err.rate[,"OOB"],
          horvath_sporadic_vhl_rf$err.rate[,"Yes"],
          horvath_sporadic_vhl_rf$err.rate[,"No"]))

ggplot(data=oob.error.horvath_sporadic_vhl_rf, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))

# Call:
#   randomForest(x = horvath_sporadic_vhl, y = horvath_sporadic_vhl$AnyMet,      ntree = 37, proximity = TRUE) 
# Type of random forest: classification
# Number of trees: 37
# No. of variables tried at each split: 18
# 
# OOB estimate of  error rate: 24.07%
# Confusion matrix:
#   No Yes class.error
# No  10   9   0.4736842
# Yes  4  31   0.1142857

### HORVATH MEN1 VHL

# Type of random forest: classification
# Number of trees: 28
# No. of variables tried at each split: 18
# 
# OOB estimate of  error rate: 21.15%
# Confusion matrix:
#   No Yes class.error
# No   8   7   0.4666667
# Yes  4  33   0.1081081

########################## LEVINE SPORADIC + MEN1 ONLY

levine_data_frame_norm$Sample_Group <- metadata_file$Sample_Group
levine_data_frame_norm$AnyMet <- as.factor(metadata_file$AnyMet)
levine_sporadic_men1 <- levine_data_frame_norm[levine_data_frame_norm$Sample_Group != 'Sporadic',]
levine_sporadic_men1_rf <- randomForest(levine_sporadic_men1, levine_sporadic_men1$AnyMet, ntree=23, proximity = TRUE)

oob.error.levine_sporadic_men1_rf <- data.frame(
  Trees=rep(1:nrow(levine_sporadic_men1_rf$err.rate), times=3),
  Type=rep(c("OOB", "Yes", "No"), each=nrow(levine_sporadic_men1_rf$err.rate)),
  Error=c(levine_sporadic_men1_rf$err.rate[,"OOB"],
          levine_sporadic_men1_rf$err.rate[,"Yes"],
          levine_sporadic_men1_rf$err.rate[,"No"]))

ggplot(data=oob.error.levine_sporadic_men1_rf, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))

# Call:
#   randomForest(x = levine_sporadic_men1, y = levine_sporadic_men1$AnyMet,      ntree = 21, proximity = TRUE) 
# Type of random forest: classification
# Number of trees: 21
# No. of variables tried at each split: 22
# 
# OOB estimate of  error rate: 22.09%
# Confusion matrix:
#   No Yes class.error
# No   6  16   0.7272727
# Yes  3  61   0.0468750

########################## LEVINE SPORADIC VHL

# Type of random forest: classification
# Number of trees: 48
# No. of variables tried at each split: 22
# 
# OOB estimate of  error rate: 25.93%
# Confusion matrix:
#   No Yes class.error
# No   8  11  0.57894737
# Yes  3  32  0.08571429

######## LEVINE MEN1 VHL

# Type of random forest: classification
# Number of trees: 23
# No. of variables tried at each split: 22
# 
# OOB estimate of  error rate: 15.38%
# Confusion matrix:
#   No Yes class.error
# No  10   5  0.33333333
# Yes  3  34  0.08108108

levine_men1 <- levine_data_frame_norm[levine_data_frame_norm$Sample_Group == 'MEN1',]
levine_men1_rf <- randomForest(levine_men1, levine_men1$AnyMet, ntree=30, proximity = TRUE)

oob.error.levine_men1_rf <- data.frame(
  Trees=rep(1:nrow(levine_men1_rf$err.rate), times=3),
  Type=rep(c("OOB", "Yes", "No"), each=nrow(levine_men1_rf$err.rate)),
  Error=c(levine_men1_rf$err.rate[,"OOB"],
          levine_men1_rf$err.rate[,"Yes"],
          levine_men1_rf$err.rate[,"No"]))

ggplot(data=oob.error.levine_men1_rf, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))

# Call:
#   randomForest(x = levine_men1, y = levine_men1$AnyMet, ntree = 30,      proximity = TRUE) 
# Type of random forest: classification
# Number of trees: 30
# No. of variables tried at each split: 22
# 
# OOB estimate of  error rate: 21.43%
# Confusion matrix:
#     No Yes class.error
# No   0   9           1
# Yes  0  33           0

levine_vhl <- levine_data_frame_norm[levine_data_frame_norm$Sample_Group == 'VHL',]
levine_vhl_rf <- randomForest(levine_vhl, levine_vhl$AnyMet, ntree=100, proximity = TRUE)

oob.error.levine_vhl_rf <- data.frame(
  Trees=rep(1:nrow(levine_vhl_rf$err.rate), times=3),
  Type=rep(c("OOB", "Yes", "No"), each=nrow(levine_vhl_rf$err.rate)),
  Error=c(levine_vhl_rf$err.rate[,"OOB"],
          levine_vhl_rf$err.rate[,"Yes"],
          levine_vhl_rf$err.rate[,"No"]))

ggplot(data=oob.error.levine_vhl_rf, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))

# Call:
#   randomForest(x = levine_vhl, y = levine_vhl$AnyMet, ntree = 100,      proximity = TRUE) 
# Type of random forest: classification
# Number of trees: 100
# No. of variables tried at each split: 22
# 
# OOB estimate of  error rate: 30%
# Confusion matrix:
#     No Yes class.error
# No   6   0        0.00
# Yes  3   1        0.75

##########################
# trying to use 6 cores
# from https://stackoverflow.com/questions/7830255/suggestions-for-speeding-up-random-forests

# library("foreach")
# library("doSNOW")
# registerDoSNOW(makeCluster(6, type="SOCK"))
# 
# # x <- matrix(runif(500), 100)
# # y <- gl(2, 50)
# 
# model_1000_mtx <- foreach(ntree = rep(250, 4), .combine = combine, .packages = "randomForest") %dopar% randomForest(x, y, ntree = ntree)
# model_1000_mtx
