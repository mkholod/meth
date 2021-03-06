# from https://www.youtube.com/watch?v=GtgJEVxl7DY and
# https://www.youtube.com/watch?v=DkLNb0CXw84 (part 2)

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

metadata_sample_group <- metadata_file$AnyMet
metadata_sample_group_len <- length(metadata_sample_group)
# testset contains around 10% of the original set
train_coef <- 0.9
train_len <- as.integer(metadata_sample_group_len * train_coef)
train_target <- as.integer(as.factor(metadata_sample_group[1:train_len]))
test_target <- as.integer(as.factor(metadata_sample_group[(train_len + 1):metadata_sample_group_len]))

csv_dir <- file.path(getwd(), "")
beta_csv_path <- file.path(csv_dir, "beta.csv") 

# lines_to_keep <- 100000
# total_lines <- R.utils::countLines(beta_csv_path) # 866092
# total_lines <- 866092
# beta <- read.csv(file=beta_csv_path, row.names = 1, skip=total_lines-lines_to_keep)
# beta <- read.csv(file=beta_csv_path, row.names = 1, nrows=300000) # remove nrows to have all the betas

beta <- read.csv(file=beta_csv_path, row.names = 1, nrows=10000)
beta_orig <- beta
# random_idx <- sample(1:866092, 100000, replace=FALSE)
# beta_idx <- cbind(seq_len(nrow(beta_orig)), beta_orig)
# beta_rand <- beta_orig[beta_idx$`seq_len(nrow(beta_orig))` %in%  random_idx, ]
# beta <- beta_rand

beta_without_na <- na.omit(beta) # remove cgX which have NA values since KNN needs to have full data

# transpose so we have a row of a single patient instead of a column
data_frame <- data.frame(t(beta_without_na))

# data_frame_pnet_only <- data_frame[metadata_file$Site == 'PNET', ] # use data_frame only for PNET site

# normalize to have distributed values between 0 and 1
data_frame_norm <- as.data.frame(lapply(data_frame, min_max_norm))
which(is.na(data_frame_norm)) # make sure there are no NA

train_set <- data_frame_norm[1:train_len, ]
test_set <- data_frame_norm[(train_len + 1):metadata_sample_group_len, ]

require(class)
classifier_knn <- knn(train=train_set, test=test_set, cl=train_target, k=3)
m1_table <- table(test_target, classifier_knn)
misClassError <- mean(classifier_knn != test_target)
print(paste('Accuracy =', 1-misClassError))

# [1] "Accuracy = 0.6"
# > m1_table
# classifier_knn
# test_target 1 2
#           1 2 3
#           2 1 4

####################### KNN on Horvath beta

csv_dir <- file.path(getwd(), "csv")
beta_csv_path <- file.path(csv_dir, "horvath_beta.csv") 

beta_horvath <- read.csv(file=beta_csv_path, row.names = 1)
beta_horvath_without_na <- na.omit(beta_horvath) # remove cgX which have NA values since KNN needs to have full data

# transpose so we have a row of a single patient instead of a column
horvath_data_frame <- data.frame(t(beta_horvath_without_na))

# data_frame_pnet_only <- data_frame[metadata_file$Site == 'PNET', ] # use data_frame only for PNET site

# normalize to have distributed values between 0 and 1
horvath_data_frame_norm <- as.data.frame(lapply(horvath_data_frame, min_max_norm))
which(is.na(horvath_data_frame_norm)) # make sure there are no NA

train_set <- horvath_data_frame_norm[1:train_len, ]
test_set <- horvath_data_frame_norm[(train_len + 1):metadata_sample_group_len, ]

require(class)
classifier_knn <- knn(train=train_set, test=test_set, cl=train_target, k=3)
m1_table <- table(test_target, classifier_knn)
misClassError <- mean(classifier_knn != test_target)
print(paste('Accuracy =', 1-misClassError))

# accuracy 0.6

####################### KNN on Levine beta

beta_csv_path <- file.path(csv_dir, "levine_beta.csv") 

beta_levine <- read.csv(file=beta_csv_path, row.names = 1)
beta_levine_without_na <- na.omit(beta_levine) # remove cgX which have NA values since KNN needs to have full data

# transpose so we have a row of a single patient instead of a column
levine_data_frame <- data.frame(t(beta_levine_without_na))

# data_frame_pnet_only <- data_frame[metadata_file$Site == 'PNET', ] # use data_frame only for PNET site

# normalize to have distributed values between 0 and 1
levine_data_frame_norm <- as.data.frame(lapply(levine_data_frame, min_max_norm))
which(is.na(levine_data_frame_norm)) # make sure there are no NA

train_set <- levine_data_frame_norm[1:train_len, ]
test_set <- levine_data_frame_norm[(train_len + 1):metadata_sample_group_len, ]

require(class)
classifier_knn <- knn(train=train_set, test=test_set, cl=train_target, k=3)
m1_table <- table(test_target, classifier_knn)
misClassError <- mean(classifier_knn != test_target)
print(paste('Accuracy =', 1-misClassError))

# accuracy 0.5
