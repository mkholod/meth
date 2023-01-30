# Probably some packages and imports need to be taken from DNAmAge_calc_full.R
library(methylclockData)
library(methylclock)
library(Biobase)
library(tibble)
library(impute)
library(ggplot2)
library(ggpmisc)
library(GEOquery)

setwd("c:/zz - private/meth/meth/meth/")

data <- read.csv("myDNAmAge_with_acceleration_with_metadata_without_pedbe_wu_tl_bnn.csv")
pnet_data <- data[data$Site == "PNET",]

plot.new()

# All tumors
p1 <- plotDNAmAge(pnet_data$Horvath, pnet_data$Age)
p1 + theme(
  plot.title = element_text(color="black", size=30, face="bold.italic", text("title")),
  axis.title.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=14, face="bold"),
) + labs(title="Horvath", x ="Horvath DNA methylation age", y = "Measured age")

p2 <- plotDNAmAge(pnet_data$Levine, pnet_data$Age, tit="Levine")
p2 + theme(
  plot.title = element_text(color="black", size=30, face="bold.italic", text("title")),
  axis.title.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=14, face="bold"),
) + labs(title="Levine", x ="Levine DNA methylation age", y = "Measured age")

filter_sporadic <- subset(pnet_data, Sample_Group == "Sporadic")
filter_MEN1 <- subset(pnet_data, Sample_Group == "MEN1")
filter_VHL <- subset(pnet_data, Sample_Group == "VHL")
rownames(filter_sporadic) <- NULL

p <- plotDNAmAge(filter_sporadic$Horvath, filter_sporadic$Age, tit="Horvarth sporadic")
p + theme(
  plot.title = element_text(color="black", size=30, face="bold.italic"),
  axis.title.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=14, face="bold"),
) + labs(title="Horvath sporadic", x ="Horvath DNA methylation age", y = "Measured age")


p3 <- plotDNAmAge(filter_MEN1$Horvath, filter_MEN1$Age, tit="Hovarth (MEN1)")
p3 + theme(
  plot.title = element_text(color="black", size=30, face="bold.italic"),
  axis.title.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=14, face="bold"),
) + labs(title="Horvath MEN1", x ="Horvath DNA methylation age", y = "Measured age")

p4 <- plotDNAmAge(filter_VHL$Horvath, filter_VHL$Age, tit="Hovarth (VHL)")
p4 + theme(
  plot.title = element_text(color="black", size=30, face="bold.italic"),
  axis.title.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=14, face="bold"),
) + labs(title="Horvath VHL", x ="Horvath DNA methylation age", y = "Measured age")

p5 <- plotDNAmAge(filter_sporadic$Levine, filter_sporadic$Age, tit="Levine (sporadic)")
p5 + theme(
  plot.title = element_text(color="black", size=30, face="bold.italic"),
  axis.title.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=14, face="bold"),
) + labs(title="Levine sporadic", x ="Levine DNA methylation age", y = "Measured age")

p6 <- plotDNAmAge(filter_MEN1$Levine, filter_MEN1$Age, tit="Levine (MEN1)")
p6 + theme(
  plot.title = element_text(color="black", size=30, face="bold.italic"),
  axis.title.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=14, face="bold"),
) + labs(title="Levine MEN1", x ="Levine DNA methylation age", y = "Measured age")

p7 <- plotDNAmAge(filter_VHL$Levine, filter_VHL$Age, tit="Levine (VHL)")
p7 + theme(
  plot.title = element_text(color="black", size=30, face="bold.italic"),
  axis.title.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=14, face="bold"),
) + labs(title="Levine VHL", x ="Levine DNA methylation age", y = "Measured age")

##############
# 
# res <- boxplot(pnet_data$ageAcc.Levine ~ pnet_data$Sample_Group)
# dat1 <- as.numeric(formatC(res$stats[,1], digits = 2, format = "f"))
# text(y=fivenum(dat1),labels=dat1,x=1.5)
# dat2 <- as.numeric(formatC(res$stats[,2], digits = 2, format = "f"))
# text(y=fivenum(dat2),labels=dat2,x=2.5)
# dat3 <- as.numeric(formatC(res$stats[,3], digits = 2, format = "f"))
# text(fivenum(dat3),labels=dat3,x=3.5)

boxplot(pnet_data$ageAcc.Levine ~ pnet_data$Sample_Group, xlab="", ylab="Age Acc Levine", notch = FALSE)
boxplot(pnet_data$ageAcc.Horvath ~ pnet_data$Sample_Group, xlab="", ylab="Age Acc Horvath", notch = FALSE)
