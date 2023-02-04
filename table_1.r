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

# setwd("c:/zz - private/meth/meth/meth")
# data_with_clock <- read.csv("myDNAmAge_with_metadata_without_pedbe_wu_tl_bnn.csv")
# names(data_with_clock)[names(data_with_clock) == "AnyMet"] <- "Metastasis"
# names(data_with_clock)[names(data_with_clock) == "Gender"] <- "Males"
# names(data_with_clock)[names(data_with_clock) == "Horvath"] <- "HMC"
# names(data_with_clock)[names(data_with_clock) == "Levine"] <- "LMC"
# 
# 
# data_with_clock$Metastasis <- data_with_clock$Metastasis == "Yes"
# data_with_clock$Males      <- data_with_clock$Males == "m"

data_with_clock$Dis[data_with_clock$Dis == ''] <- 'No'
dis_table1 <- table1(~ Age + HMC + LMC + Males + Metastasis  | 
                   Sample_Group + Dis, 
                 render = rndr,
                 render.continuous=my.render.cont,
                 data=data_with_clock)
dis_table1
write.csv(table1, "export_PNET/table1_dis.csv")


table1 <- table1(~ Age + HMC + LMC + Males + Metastasis  | 
         Sample_Group, 
       render = rndr,
       render.continuous=my.render.cont,
       data=data_with_clock)
table1
write.csv(table1, "export_PNET/table1_all.csv")



####
data_with_clock$LN[data_with_clock$LN == ''] <- 'No'
ln_table1 <- table1(~ Age + HMC + LMC + Males + Metastasis  | 
                       Sample_Group + LN, 
                     render = rndr,
                     render.continuous=my.render.cont,
                     data=data_with_clock)
ln_table1
write.csv(table1, "export_PNET/table1_ln.csv")

#### chi test
# https://www.r-tutor.com/elementary-statistics/goodness-fit/chi-squared-test-independence
library(MASS)
dis_tbl = table(data_with_clock$Age, data_with_clock$Dis)
chisq.test(dis_tbl) 

# Pearson's Chi-squared test
# 
# data:  dis_tbl
# X-squared = 69.28, df = 34, p-value = 0.000332
# 
# Warning message:
# In chisq.test(dis_tbl) : Chi-squared approximation may be incorrect

dis_tbl_hmc = table(data_with_clock$HMC, data_with_clock$Dis)
chisq.test(dis_tbl_hmc) 

# Pearson's Chi-squared test
# 
# data:  dis_tbl_hmc
# X-squared = 96, df = 95, p-value = 0.452

dis_tbl_lmc = table(data_with_clock$LMC, data_with_clock$Dis)
chisq.test(dis_tbl_lmc) 

# Pearson's Chi-squared test
# 
# data:  dis_tbl_lmc
# X-squared = 96, df = 95, p-value = 0.452

data_with_clock_men1 <- data_with_clock[data_with_clock$Sample_Group == "MEN1",]
dis_tbl_men1 = table(data_with_clock_men1$Age, data_with_clock_men1$Dis)
chisq.test(dis_tbl_men1) 

# Pearson's Chi-squared test
# 
# data:  dis_tbl_men1
# X-squared = 42, df = 12, p-value = 3.331e-05

dis_tbl_men1_hmc = table(data_with_clock_men1$HMC, data_with_clock_men1$Dis)
chisq.test(dis_tbl_men1_hmc) 

# Pearson's Chi-squared test
# 
# data:  dis_tbl_men1_hmc
# X-squared = 42, df = 41, p-value = 0.4274

dis_tbl_men1_lmc = table(data_with_clock_men1$LMC, data_with_clock_men1$Dis)
chisq.test(dis_tbl_men1_lmc) 

# Pearson's Chi-squared test
# 
# data:  dis_tbl_men1_lmc
# X-squared = 42, df = 41, p-value = 0.4274

data_with_clock_sporadic <- data_with_clock[data_with_clock$Sample_Group == "Sporadic",]
dis_tbl_sporadic = table(data_with_clock_sporadic$Age, data_with_clock_sporadic$Dis)
chisq.test(dis_tbl_sporadic) 

# Pearson's Chi-squared test
# 
# data:  dis_tbl_sporadic
# X-squared = 35.933, df = 28, p-value = 0.1443

dis_tbl_sporadic_hmc = table(data_with_clock_sporadic$HMC, data_with_clock_sporadic$Dis)
chisq.test(dis_tbl_sporadic_hmc) 

# Pearson's Chi-squared test
# 
# data:  dis_tbl_sporadic_hmc
# X-squared = 44, df = 43, p-value = 0.429

dis_tbl_sporadic_lmc = table(data_with_clock_sporadic$LMC, data_with_clock_sporadic$Dis)
chisq.test(dis_tbl_sporadic_lmc) 

# Pearson's Chi-squared test
# 
# data:  dis_tbl_sporadic_lmc
# X-squared = 44, df = 43, p-value = 0.429


## ++ Export table1 to csv
## ++ Export table1_pnet to csv
## ++ Not round the values of table1 - 1 digit after the point
## Correlation graphs should have pvalue and n
# ++ Hovarth should be HMC in graph, Levine LMC, Y should be chronological age
# Box plot should be on PNET and on everything  

