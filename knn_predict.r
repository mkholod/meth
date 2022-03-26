# from https://www.youtube.com/watch?v=GtgJEVxl7DY and
# https://www.youtube.com/watch?v=DkLNb0CXw84 (part 2)

# Predict patient age out of beta

csv_dir <- file.path(getwd(), "csv")
horvath_beta_csv_path <- file.path(csv_dir, "horvath_beta.csv")
horvath_beta <- read.csv(file=horvath_beta_csv_path, row.names = 1)

# transpose so we have a row of a single patient
data_frame <- data.frame(t(horvath_beta))

# Normalize the values to be between 0 and 1
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

data_frame_norm <- as.data.frame(lapply(data_frame, min_max_norm))

# testset contains around 10% of the original set
train_set <- data_frame_norm[1:85, ]
test_set <- data_frame_norm[86:96, ]
train_target <- metadata_df$Age[1:85]
test_target <- metadata_df$Age[86:96]

require(class)
m1 <- knn(train=train_set, test=test_set, cl=train_target, k=11)
m1_table <- table(test_target, m1)

# > table(test_target, m1)
# m1
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

