# from https://www.youtube.com/watch?v=GtgJEVxl7DY and
# https://www.youtube.com/watch?v=DkLNb0CXw84 (part 2)

# Normalize the values to be between 0 and 1
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

baseDir <- getwd()
sample_sheet_filename = "/Sample_Sheet_Full.csv"
path <- paste(baseDir, sample_sheet_filename, sep="")
dataOnly = read.csv(path, skip = 8, header = F)
headers = read.csv(path, skip = 7, header = F, nrows = 1, as.is = T)
colnames(dataOnly) = headers
metadata_file <- dataOnly
metadata_file_pnet_only <- metadata_file[metadata_file$Site == 'PNET', ]

metadata_pnet_sample_group <- metadata_file_pnet_only$Sample_Group
metadata_sample_group_len <- length(metadata_pnet_sample_group)
# testset contains around 10% of the original set
train_len <- as.integer(metadata_sample_group_len * 0.89)
train_target <- as.integer(as.factor(metadata_pnet_sample_group[1:train_len]))
test_target <- as.integer(as.factor(metadata_pnet_sample_group[(train_len + 1):metadata_sample_group_len]))

csv_dir <- file.path(getwd(), "")
beta_csv_path <- file.path(csv_dir, "beta.csv") 
# beta <- read.csv(file=beta_csv_path, row.names = 1, nrows=100000) # remove nrows to have all the betas
beta <- read.csv(file=beta_csv_path, row.names = 1) 
beta_without_na <- na.omit(beta) # remove cgX which have NA values since KNN needs to have full data

# transpose so we have a row of a single patient instead of a column
data_frame <- data.frame(t(beta_without_na))

# use data_frame only for PNET site
data_frame_pnet_only <- data_frame[metadata_file$Site == 'PNET', ]

# normalize to have distributed values between 0 and 1
data_frame_pnet_only_norm <- as.data.frame(lapply(data_frame_pnet_only, min_max_norm))

train_set <- data_frame_pnet_only_norm[1:train_len, ]
test_set <- data_frame_pnet_only_norm[(train_len + 1):metadata_sample_group_len, ]

require(class)
classifier_knn <- knn(train=train_set, test=test_set, cl=train_target, k=3)
m1_table <- table(test_target, classifier_knn)
misClassError <- mean(classifier_knn != test_target)
print(paste('Accuracy =', 1-misClassError))

# k=1 "Accuracy = 0.272727272727273"
# k=3,5,7, 13 acc=0
# k=9, 11, 17, 19, 23 "Accuracy = 0.0909090909090909"
# k=15 "Accuracy = 0.181818181818182"

# > table(test_target, classifier_knn)
# classifier_knn
# test_target 26 28 32 34 36 37 39 41 42 44 47 48 50 51 53 54 55 56 58 59 60 62 63 65 66 67 68 69 70 72 73 78 80
#          28  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0
#          37  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#          42  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#          43  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#          45  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0
#          53  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#          55  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#          56  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#          63  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  2  0  0  0  0  0  0  0
#          69  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0

# add the actual age as another column
# sample_sheet_filename_path <- file.path(getwd(), "Sample_Sheet_Full.csv")
# metadata_df <- read.csv(sample_sheet_filename_path, skip = 7)
# data_frame$Age <- metadata_df$Age

#############################################################
# or with the DNA age and predict if it has metastasis or not

