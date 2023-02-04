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

rndr <- function(x, ...) {
  y <- render.default(x, ...)
  if (is.logical(x)) y[2] else y
}

my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=3), c("",
                                                           "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}

setwd("c:/zz - private/meth/meth/meth")
data_with_clock <- read.csv("myDNAmAge_with_metadata_without_pedbe_wu_tl_bnn.csv")
names(data_with_clock)[names(data_with_clock) == "AnyMet"] <- "Metastasis"
names(data_with_clock)[names(data_with_clock) == "Gender"] <- "Males"
names(data_with_clock)[names(data_with_clock) == "Horvath"] <- "HMC"
names(data_with_clock)[names(data_with_clock) == "Levine"] <- "LMC"

data_with_clock$Metastasis <- data_with_clock$Metastasis == "Yes"
data_with_clock$Males      <- data_with_clock$Males == "m"

pnet_data_with_clock <- data_with_clock[data_with_clock$Site == "PNET",]
table1_pnet <- table1(~ Age + HMC + LMC + Males + Metastasis  | 
                   Sample_Group, 
                 render = rndr,
                 render.continuous=my.render.cont,
                 data=pnet_data_with_clock)
table1_pnet
write.csv(table1_pnet, "export_PNET/table1_pnet.csv")

# For all the samples
## Add calculated age by Horvath and Levine (only mean?)
## Display only mean
## Display only male
## Display only metastasis yes

setwd("c:/zz - private/meth/meth/meth")
data_with_clock <- read.csv("myDNAmAge_with_metadata_without_pedbe_wu_tl_bnn.csv")
names(data_with_clock)[names(data_with_clock) == "AnyMet"] <- "Metastasis"
names(data_with_clock)[names(data_with_clock) == "Gender"] <- "Males"
names(data_with_clock)[names(data_with_clock) == "Horvath"] <- "HMC"
names(data_with_clock)[names(data_with_clock) == "Levine"] <- "LMC"


data_with_clock$Metastasis <- data_with_clock$Metastasis == "Yes"
data_with_clock$Males      <- data_with_clock$Males == "m"

table1 <- table1(~ Age + HMC + LMC + Males + Metastasis  | 
         Sample_Group, 
       render = rndr,
       render.continuous=my.render.cont,
       data=data_with_clock)
table1
write.csv(table1, "export_PNET/table1_all.csv")

# dis column - empty values are no,
dis_yes <- data_with_clock[data_with_clock$Dis == "Yes",]
dis_no <- data_with_clock[data_with_clock$Dis != "Yes",]

# percentage of yes vs no
print(nrow(dis_yes) / nrow(data_with_clock))
## 29/96 = 0.3020833

# what is the average age of dis=yes
print(result.mean_age <- mean(dis_yes$Age))
# 57

# what is the average hovarth and levine for dis=yes
print(result.mean_age <- mean(dis_yes$HMC)) 
# 64.90531

print(result.mean_age <- mean(dis_yes$LMC)) 
# 44.36717

# what is the average age of dis=yes
print(result.mean_age <- mean(dis_no$Age))
# 51.04478

# what is the average hovarth and levine for dis=yes
print(result.mean_age <- mean(dis_no$HMC)) 
# 60.80572

print(result.mean_age <- mean(dis_no$LMC)) 
# 29.63403


####
# dis column - empty values are no,
ln_yes <- data_with_clock[data_with_clock$LN == "Yes",]
ln_no <- data_with_clock[data_with_clock$LN != "Yes",]

# percentage of yes vs no
print(nrow(ln_yes) / nrow(data_with_clock))
## 0.53125

# what is the average age of dis=yes
print(result.mean_age <- mean(ln_yes$Age))
# 54.88235

# what is the average hovarth and levine for dis=yes
print(result.mean_age <- mean(ln_yes$HMC)) 
# 62.1113

print(result.mean_age <- mean(ln_yes$LMC)) 
# 31.28166

# what is the average age of dis=yes
print(result.mean_age <- mean(ln_no$Age))
# 50.53333

# what is the average hovarth and levine for dis=yes
print(result.mean_age <- mean(ln_no$HMC)) 
# 61.96802

print(result.mean_age <- mean(ln_no$LMC)) 
# 37.26141

## ++ Export table1 to csv
## ++ Export table1_pnet to csv
## ++ Not round the values of table1 - 1 digit after the point
## Correlation graphs should have pvalue and n
# ++ Hovarth should be HMC in graph, Levine LMC, Y should be chronological age
# Box plot should be on PNET and on everything  

