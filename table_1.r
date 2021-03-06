# https://cran.r-project.org/web/packages/table1/index.html
# https://cran.r-project.org/web/packages/table1/table1.pdf
# https://github.com/benjaminrich/table1

install.packages("table1")
require(table1)
require(survival)

# sample_sheet_filename = "/Sample_Sheet_Full.csv"
# 
# path <- paste(baseDir, sample_sheet_filename, sep="")
# dataOnly = read.csv(path, skip = 8, header = F)
# headers = read.csv(path, skip = 7, header = F, nrows = 1, as.is = T)
# 
# colnames(dataOnly) = headers
# t1 <- dataOnly
# t1$Metastasis <- factor(t1$AnyMet, 
#                                     levels=c('Yes','No'),
#                                     labels=c('Yes', 'No' # Reference
#                                              ))
# 
# table1(~ Age + Gender + Metastasis | Sample_Group , data=t1)

##### Now with Horvath and Levine ages

data_with_clock <- read.csv("csv/myDNAmAge_with_acceleration_age_normalized_with_metadata.csv")
names(data_with_clock)[names(data_with_clock) == "AnyMet"] <- "Metastasis"
table1(~ Age + Gender + Metastasis + Horvath + Levine | Sample_Group , data=data_with_clock)
