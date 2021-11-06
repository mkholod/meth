library(dplyr)
library(minfi)
library(RColorBrewer)
library(reshape2)
library(plotly)
library(webshot)
library(ggplot2)
library("org.Hs.eg.db")
library("annotate")
library("biomaRt")
library(doParallel)
library(annotatr)
library(xlsx)
library(Gviz)
library(Vennerable) #venn 
library(ConsensusClusterPlus)
library("clusterProfiler")
library("ReactomePA")
library("gplots")
library(heatmap.plus)
library(genefilter)
library(outliers)
library(ggpubr)
library(pathfindR) #Pathway analysis

#memory.limit(size=10000)

#BiocManager::install("genefilter")

setwd("C:/Users/tiros/Documents/Work/Projects/Methylation/raw data")

#######################
# minfi  data analaysis
#######################

baseDir <- ("C:/Users/tiros/Documents/Work/Projects/Methylation/raw data")
targets <- read.metharray.sheet(baseDir, pattern = "Sample_Sheet.csv")
#saveRDS(targets, "targets_96.rds")

#subset VHL
#targets <- targets[c(5, 27, 35, 43, 52, 60, 75, 83, 84, 92),]

#subset MEN1
#targets <- targets[c(2,  3,  4,  6, 10, 11, 12, 13, 14, 19, 20, 21, 22, 26, 28, 29, 30, 34, 36, 37, 38, 42, 44, 45, 46, 50, 51, 53, 59, 61, 67, 68, 69, 74, 76, 77, 82, 85, 89, 90, 91, 93),]
#
#subset PNET
#targets <- targets[targets$Site=="PNET",]
#targets <- droplevels(targets)

#subset Sporadic
#targets <- targets[c(1,  7 , 8 , 9 ,15 ,16 ,17 ,18 ,23 ,24, 25, 31, 32, 33, 39, 40, 41, 47, 48, 49, 54, 55, 56, 57, 58, 62, 63, 64, 65, 66, 70, 71, 72, 73, 78, 79, 80, 81, 86, 87, 88, 94, 95, 96),]
#targets <- droplevels(targets)
#targets <- droplevels(targets)

rgSet <- read.metharray.exp(targets = targets, force=TRUE)  #rgSet unprocessed format
#saveRDS(rgSet, "rgSet_96.rds")
#saveRDS(rgSet, "rgSet_MEN1.rds")
#saveRDS(rgSet, "rgSet_PNET.rds")

#rgSet <- readRDS("rgSet_96.rds")

gset.funnorm <- preprocessFunnorm(rgSet)

#saveRDS(gset.funnorm, "gsetFunnorm_96.rds")
#saveRDS(gset.funnorm, "gsetFunnorm_MEN1.rds")
#saveRDS(gset.funnorm, "gsetFunnorm_NI.rds")
#saveRDS(gset.funnorm, "gsetFunnorm_PNET.rds")

#gset.funnorm <- readRDS("gsetFunnorm_96.rds")
#gset.funnorm <- readRDS("gsetFunnorm_MEN1.rds")
#gset.funnorm <- readRDS("gsetFunnorm_NI.rds")

#plotSex(getSex(gset.funnorm, cutoff = -2))

#QC plots
#densityPlot(rgSet, sampGroups = targets$Sample_Group)
#densityBeanPlot(rgSet, sampGroups = targets$Sample_Group)

#SNPs
#Identify
#snps <- getSnpInfo(gset.funnorm)
#add data to the ratio values matrix
#gset.funnorm <- addSnpInfo(gset.funnorm)
#Remove SNPs
gset.funnorm <- dropLociWithSnps(gset.funnorm, snps=c("SBE","CpG"))

#saveRDS(gset.funnorm, "gsetFunnorm_96_DropLoci.rds")
#saveRDS(gset.funnorm, "gsetFunnorm_MEN1_DropLoci.rds")
#saveRDS(gset.funnorm, "gsetFunnorm_NI_DropLoci.rds")

gset.funnorm <- readRDS("gsetFunnorm_96_DropLoci.rds")

#probe locations as a genomic ranges
gr <- granges(gset.funnorm)

#DMP finder
Island_data <- getAnnotation(gset.funnorm)

#saveRDS(Island_data, "island_data_96.rds") # not saved! (saved gsetfunnorm instead)
#saveRDS(Island_data, "island_data_MEN1.rds")  # not saved! (saved gsetfunnorm instead)
#saveRDS(Island_data, "island_data_NI.rds")  # not saved! (saved gsetfunnorm instead)
#saveRDS(Island_data, "island_data_PNET.rds")  # ok

Island_data <- readRDS("island_data_PNET.rds")
View(as.data.frame(Island_data))

beta <- getBeta(gset.funnorm) #gets matrix

#saveRDS(beta, "beta_96.rds")
#saveRDS(beta, "beta_MEN1.rds")
#saveRDS(beta, "beta_NI.rds")
#saveRDS(beta, "beta_PNET.rds")

#beta_samples <- readRDS("beta_96.rds")
#beta_rotated <- t(beta_samples)
#write.csv(beta_rotated, "all_beta.csv")

#setwd("C:/Users/tiros/Documents/Work/Projects/Methylation/raw data NI")
#beta_NI <- readRDS("beta_NI.rds")
#setwd("C:/Users/tiros/Documents/Work/Projects/Methylation/raw data")

#beta_100 <- cbind(beta_samples, beta_NI)
#saveRDS(beta_100, "beta_100.rds")

#dd<- beta_samples[c(1:100),]
#outlierTest(dd)
beta_clean <- rm.outlier(beta_samples, fill = TRUE, median = TRUE, opposite = FALSE)

#M <- M[startsWith(rownames(M), "cg"), ]

#statistical comparisons
b_Sporadic_mean <- 0 
b_MEN1_mean <- 0
b_VHL_mean <- 0

mat <- matrix(, ncol = 8, nrow = nrow(beta_samples))

mat[,1] <- rownames(beta_samples)
beta_samples <- beta_clean
for(i in 1:nrow(beta_samples)){
  b_Sporadic <- beta_samples[i, c(1,  7 , 8 , 9 ,15 ,16 ,17 ,18 ,23 ,24, 25, 31, 32, 33, 39, 40, 41, 47, 48, 49, 54, 55, 56, 57, 58, 62, 63, 64, 65, 66, 70, 71, 72, 73, 78, 79, 80, 81, 86, 87, 88, 94, 95, 96)]
  b_MEN1 <- beta_samples[i, c(2,  3,  4,  6, 10, 11, 12, 13, 14, 19, 20, 21, 22, 26, 28, 29, 30, 34, 36, 37, 38, 42, 44, 45, 46, 50, 51, 53, 59, 61, 67, 68, 69, 74, 76, 77, 82, 85, 89, 90, 91, 93)]
  b_VHL <- beta_samples[i, c(5, 27, 35, 43, 52, 60, 75, 83, 84, 92)]
 #PNETs:
  #b_Sporadic <- beta_samples[i, c(4, 7, 9, 13, 18, 23, 28, 31, 34, 35, 39, 40, 44, 45, 49, 50)]
  #b_MEN1 <- beta_samples[i, c(1, 2, 3, 5, 6, 8, 10, 12, 14, 16, 17, 19, 31, 32, 33, 36, 38, 41, 46, 47)]
  #b_VHL <- beta_samples[i, c(3, 11, 15, 20, 26, 30, 37, 42, 43, 48)]
  #b_NI <- beta_NI[i, 1:4]                         
  #Means calculations
  mat[i,2] <- mean(as.numeric(b_Sporadic))
  mat[i,3] <- mean(as.numeric(b_MEN1))
  mat[i,4] <- mean(as.numeric(b_VHL))
}

design

data("mtcars")
mtcars %>% 
  group_by_(.dots=c("mpg","hp","wt")) %>% 
  summarize(x=mean(gear))
View(mtcars)


saveRDS(mat, "mat_all,rds")
  #Comparisons
  #Sporadic vs. NI
  #MEN1 vs. NI
  #VHL vs. NI
  mat[i,6] <- t.test(x = b_Sporadic,
         y = b_NI,
         alternative = "two.sided",
         paired = FALSE,
         var.equal = FALSE)$p.value
  mat[i,7] <- t.test(x = b_MEN1,
                            y = b_NI,
                            alternative = "two.sided",
                            paired = FALSE,
                            var.equal = FALSE)$p.value
  mat[i,8] <- t.test(x = b_VHL,
                            y = b_NI,
                            alternative = "two.sided",
                            paired = FALSE,
                            var.equal = FALSE)$p.value
}

#saveRDS(mat, "mat_100.rds")
#saveRDS(mat, "mat_PNET+NI.rds")
#setwd("C:/Users/tiros/Documents/Work/Projects/Methylation/raw data")

#mat <- readRDS("mat_PNET+NI.rds")
#p_val_spo <- mat[,6]
#p_val_MEN1 <- mat[,7]
#p_val_VHL <- mat[,8]

#FDR_spo <- p.adjust(p_val_spo, method = "fdr", n = 835424)
#FDR_MEN1 <- p.adjust(p_val_MEN1, method = "fdr", n = 835424)
#FDR_VHL <- p.adjust(p_val_VHL, method = "fdr", n = 835424)

#mat <- cbind(mat, p_val_spo, p_val_MEN1, p_val_VHL)
ncol(mat)
mat1 <- as.data.frame(mat, stringsAsFactors = FALSE)
colnames(mat) <- c("CpG", "b_Sporadic", "b_MEN1", "b_VHL", "Spo_hyper_1y2n", "Spo_hypo_1y2n", "MEN1_hyper_1y2n", "MEN1_hypo_1y2n") #, "VHL_hyper_1y2n", "VHL_hypo_1y2n", "empty7")

#saveRDS(mat, "mat_9132019.rds")
mat <- readRDS("mat_9132019.rds") # final matrix - all samples beta values

## Hypermethylation
mat_spo_hyper <- mat[levels(mat$b_Sporadic) > 0.8,] # Sporadic
mat_MEN1_hyper <- mat[levels(mat$b_MEN1) > 0.8,] # MEN1
mat_VHL_hyper <- mat[levels(mat$b_VHL) > 0.8,] # VHL

## Hypomethylation
mat_spo_hypo <- mat[levels(mat$b_Sporadic) < 0.2 ,] # Sporadic
mat_MEN1_hypo <- mat[levels(mat$b_MEN1)  < 0.2 ,] # MEN1
mat_VHL_hypo <- mat[levels(mat$b_VHL) < 0.2 ,] # VHL

CG_spo_hyper <- mat_spo_hyper[,1] %>% as.character() 
CG_MEN1_hyper <- mat_MEN1_hyper[,1]  %>% as.character()
CG_VHL_hyper <- mat_VHL_hyper[,1] %>% as.character()
CG_spo_hypo <- mat_spo_hypo[,1] %>% as.character()
CG_MEN1_hypo <- mat_MEN1_hypo[,1] %>% as.character()
CG_VHL_hypo <- mat_VHL_hypo[,1] %>% as.character()

spohyper <- length(CG_spo_hyper) #41882
spohypo <- length(CG_spo_hypo)  #4781
MEN1hyper <- length(CG_MEN1_hyper)#36402
MEN1hypo <- length(CG_MEN1_hypo) #4502
VHLhyper <- length(CG_VHL_hyper) #32168
VHLhypo <- length(CG_VHL_hypo)  #7789

n_cpg <- 835424 #nrow(beta_samples)


#hyper spo MEN1
spo_men1_hyper <-matrix(c(spohyper, n_cpg-spohyper, MEN1hyper, n_cpg - MEN1hyper),ncol=2,byrow=TRUE)
spo_men1_hypo <-matrix(c(spohypo, n_cpg-spohypo, MEN1hypo, n_cpg - MEN1hypo),ncol=2,byrow=TRUE)
spo_VHL_hyper <-matrix(c(spohyper, n_cpg-spohyper, VHLhyper, n_cpg - VHLhyper),ncol=2,byrow=TRUE)
spo_VHL_hypo <-matrix(c(spohypo, n_cpg-spohypo, VHLhypo, n_cpg - VHLhypo),ncol=2,byrow=TRUE)
VHL_men1_hyper <-matrix(c(VHLhyper, n_cpg-VHLhyper, MEN1hyper, n_cpg - MEN1hyper),ncol=2,byrow=TRUE)
VHL_men1_hypo <-matrix(c(VHLhypo, n_cpg-VHLhypo, MEN1hypo, n_cpg - MEN1hypo),ncol=2,byrow=TRUE)
fisher.test(spo_men1_hyper,alternative="two.sided")$p.value * 12
fisher.test(spo_men1_hypo,alternative="two.sided")$p.value * 12
fisher.test(spo_VHL_hyper,alternative="two.sided")$p.value * 12
fisher.test(spo_VHL_hypo,alternative="two.sided")$p.value * 12
fisher.test(VHL_men1_hyper,alternative="two.sided")$p.value * 12
fisher.test(VHL_men1_hypo,alternative="two.sided")$p.value * 12


## Venn figure

CG_spo_hyper <- CG_spo_hyper1
CG_MEN1_hyper <- CG_MEN1_hyper1
CG_VHL_hyper <- CG_VHL_hyper1
CG_spo_hypo <- CG_spo_hypo1
CG_MEN1_hypo <- CG_MEN1_hypo1
CG_VHL_hypo <- CG_VHL_hypo1


hyper <- list(Sporadic = CG_spo_hyper, MEN1 = CG_MEN1_hyper, VHL = CG_VHL_hyper)

hypo <- list(Sporadic = CG_spo_hypo, MEN1 = CG_MEN1_hypo, VHL = CG_VHL_hypo)
#saveRDS(hyper, "hyper.rds")
#saveRDS(hypo, "hypo.rds")

v1 <- Venn(hyper)
saveRDS(v1, "v1_revision.rds")
v1
v2 <- Venn(hypo)
saveRDS(v2, "v2_revision.rds")
v2

dev.off()
plot(v1, doWeights = FALSE, type = "circles", show = list(SetLabels = TRUE, Faces = FALSE))
plot(v2, doWeights = FALSE, type = "circles", show = list(SetLabels = TRUE, Faces = FALSE))

dev.copy(png, "venn_diagram.png")
dev.off()


#############################
### CpGs number by tumor type
#############################
#which(targets$Site == "PNET")
#3  4  5  8 11 12 16 19 24 26 27 30 32 34 35 37 38 40 42 43 45 46 48 50 51 52 53 56 59 60 64 67 68 71 72 74 75 76 79 80 82 83 84 87 88 90 91 92 95 96

#which(targets$Site == "SINET")
#1  7  9 15 17 18 23 25 31 33 39 41 49 54 57 62 63 66 73 78 81 86 94

#which(targets$Site == "DNET")
#2  6 10 13 14 20 21 22 28 29 36 44 47 55 77 85 89 93

#statistical comparisons
b_PNET_mean <- 0 
b_SINET_mean <- 0
b_DNET_mean <- 0

b_ttest_PNET <- 0
b_ttest_SINET <- 0
b_ttest_DNET <- 0

mat_site <- matrix(, ncol = 8, nrow = nrow(beta_samples))

mat_site[,1] <- rownames(beta_samples)
for(i in 1:nrow(beta_samples)){
  b_PNET <- beta_samples[i, c(3,  4,  5,  8, 11, 12, 16, 19, 24, 26, 27, 30, 32, 34, 35, 37, 38, 40, 42, 43, 45, 46, 48, 50, 51, 52, 53, 56, 59, 60, 64, 67, 68, 71, 72, 74, 75, 76, 79, 80, 82, 83, 84, 87, 88, 90, 91, 92, 95, 96)]
  b_SINET <- beta_samples[i, c(1,  7,  9, 15, 17, 18, 23, 25, 31, 33, 39, 41, 49, 54, 57, 62, 63, 66, 73, 78, 81, 86, 94)]
  b_DNET <- beta_samples[i, c(2,  6, 10, 13, 14, 20, 21, 22, 28, 29, 36, 44, 47, 55, 77, 85, 89, 93)]
  b_NI <- beta_NI[i, 1:4]  
  #Means calculations
  mat_site[i,2] <- mean(as.numeric(b_PNET))/mean(as.numeric(b_NI))
  mat_site[i,3] <- mean(as.numeric(b_SINET))/mean(as.numeric(b_NI))
  mat_site[i,4] <- mean(as.numeric(b_DNET))/mean(as.numeric(b_NI))
  mat_site[i,5] <- mean(as.numeric(b_NI))
  
  #Comparisons
  #PNET vs. NI
  #SINET vs. NI
  #DNET vs. NI
  mat_site[i,6] <- t.test(x = b_PNET,
                     y = b_NI,
                     alternative = "two.sided",
                     paired = FALSE,
                     var.equal = FALSE)$p.value
  mat_site[i,7] <- t.test(x = b_SINET,
                     y = b_NI,
                     alternative = "two.sided",
                     paired = FALSE,
                     var.equal = FALSE)$p.value
  mat_site[i,8] <- t.test(x = b_DNET,
                     y = b_NI,
                     alternative = "two.sided",
                     paired = FALSE,
                     var.equal = FALSE)$p.value
}

#saveRDS(mat, "mat_100_sites.rds")
saveRDS(mat_site, "mat_100_site.rds")

setwd("C:/Users/tiros/Documents/Work/Projects/Methylation/raw data")
#mat <- readRDS("mat_100.rds") #comparisons by groups: sporadic, MEN1 and VHL
#mat_site <- readRDS("mat_100_sites.rds")
p_val_PNET <- mat_site[,6]
p_val_SINET <- mat_site[,7]
p_val_DNET <- mat_site[,8]

FDR_PNET <- p.adjust(p_val_PNET, method = "fdr", n = 835424)
FDR_SINET <- p.adjust(p_val_SINET, method = "fdr", n = 835424)
FDR_DNET <- p.adjust(p_val_DNET, method = "fdr", n = 835424)

mat_site <- cbind(mat_site, p_val_PNET, p_val_SINET, p_val_DNET)
ncol(mat_site)

colnames(mat_site) <- c("CpG", "b_PNET/NI", "b_SINET/NI", "b_DNET/NI", "b_NI", "b_ttest_PNET_NI", "b_ttest_SINET_NI", "b_ttest_DNET_NI", "FDR_PNET", "FDR_SINET", "FDR_DNET")

## Hypermethylation
mat_PNET_hyper <- mat_site[mat_site[,2] > 2 & mat_site[,9] <0.05,] # PNET
mat_SINET_hyper <- mat_site[mat_site[,3] > 2 & mat_site[,10] <0.05,] # SINET
mat_DNET_hyper <- mat_site[mat_site[,4] > 2 & mat_site[,11] <0.05,] # DNET

## Hypomethylation
mat_PNET_hypo <- mat_site[mat_site[,2] < 0.5 & mat_site[,9] <0.05,] # PNET
mat_SINET_hypo <- mat_site[mat_site[,3] < 0.5 & mat_site[,10] <0.05,] # SINET
mat_DNET_hypo <- mat_site[mat_site[,4] < 0.5 & mat_site[,11] <0.05,] # DNET

CG_PNET_hyper <- mat_PNET_hyper[,1]
CG_SINET_hyper <- mat_SINET_hyper[,1]
CG_DNET_hyper <- mat_DNET_hyper[,1]
CG_PNET_hypo <- mat_PNET_hypo[,1]
CG_SINET_hypo <- mat_SINET_hypo[,1]
CG_DNET_hypo <- mat_DNET_hypo[,1]

PNEThyper <- length(CG_PNET_hyper) #32256
PNEThypo <- length(CG_PNET_hypo)  #4705
SINEThyper <- length(CG_SINET_hyper)#46959
SINEThypo <- length(CG_SINET_hypo) #5890
DNEThyper <- length(CG_DNET_hyper) #40731
DNEThypo <- length(CG_DNET_hypo)  #5388

n_cpg <- nrow(beta_samples)

#Calculating n CpGs and comparing by site
PNET_SINET_hyper <-matrix(c(SINEThyper, n_cpg-SINEThyper, PNEThyper, n_cpg-PNEThyper),ncol=2,byrow=TRUE)
PNET_SINET_hypo <-matrix(c(PNEThypo, n_cpg-PNEThypo, SINEThypo, n_cpg - SINEThypo),ncol=2,byrow=TRUE)
PNET_DNET_hyper <-matrix(c(PNEThyper, n_cpg-PNEThyper, DNEThyper, n_cpg - DNEThyper),ncol=2,byrow=TRUE)
PNET_DNET_hypo <-matrix(c(PNEThypo, n_cpg-PNEThypo, DNEThypo, n_cpg - DNEThypo),ncol=2,byrow=TRUE)
DNET_SINET_hyper <-matrix(c(DNEThyper, n_cpg-DNEThyper, SINEThyper, n_cpg - SINEThyper),ncol=2,byrow=TRUE)
DNET_SINET_hypo <-matrix(c(DNEThypo, n_cpg-DNEThypo, SINEThypo, n_cpg - SINEThypo),ncol=2,byrow=TRUE)
fisher.test(PNET_SINET_hyper,alternative="two.sided")$p.value * 12
fisher.test(PNET_SINET_hypo,alternative="two.sided")$p.value * 12
fisher.test(PNET_DNET_hyper,alternative="two.sided")$p.value * 12
fisher.test(PNET_DNET_hypo,alternative="two.sided")$p.value * 12
fisher.test(DNET_SINET_hyper,alternative="two.sided")$p.value * 12
fisher.test(DNET_SINET_hypo,alternative="two.sided")$p.value * 12














DMR_MEN1_hyper <- Island_data[Island_data$Name %in% CG_MEN1_hyper & Island_data$DMR == "DMR", ]
View(as.data.frame(DMR_MEN1_hyper))

DMR_VHL_hypo <- Island_data[Island_data$Name %in% CG_VHL_hypo & Island_data$DMR == "DMR", ]
DMP_VHL_hypo <- Island_data[Island_data$Name %in% CG_VHL_hypo, ]
DMP_VHL_hypo_df <- as.data.frame(DMP_VHL_hypo)
DMP_VHL_hypo_df_genes <- gsub(";.*", "", DMP_VHL_hypo_df$UCSC_RefGene_Name)
write.xlsx(DMP_VHL_hypo_df_genes, "DMP_VHL_hypo_genes.xlsx")  

DMR_VHL_hypo_df <- as.data.frame(DMR_VHL_hypo)
DMR_VHL_hypo_df_genes <- gsub(";.*", "", DMR_VHL_hypo_df$UCSC_RefGene_Name)

DMR_VHL_hypo_df_promoter <- DMR_VHL_hypo_df[grep("Promoter", DMR_VHL_hypo_df$Regulatory_Feature_Group),]
DMR_VHL_hypo_df_promoter_genes <- gsub(";.*", "", DMR_VHL_hypo_df_promoter$UCSC_RefGene_Name)


write.xlsx(DMR_VHL_hypo_df_genes, "VHL_hypo_genes.xlsx")  
write.xlsx(DMR_VHL_hypo_df_promoter_genes, "VHL_hypo_promoter_genes.xlsx")  

## Venn figure

hyper <- list(VHL = levels(VHL_hyper), MEN1 = levels(MEN1_hyper), Sporadic = levels(spo_hyper))
hypo <- list(VHL = levels(VHL_hypo), MEN1 = levels(MEN1_hypo), Sporadic = levels(spo_hypo))

v1 <- Venn(hyper)
saveRDS(v1, "v1.rds")
v1
v2 <- Venn(hypo)
saveRDS(v2, "v2.rds")
v2

dev.off()
plot(v1, doWeights = FALSE, type = "circles", show = list(SetLabels = TRUE, Faces = FALSE))
plot(v2, doWeights = FALSE, type = "circles", show = list(SetLabels = TRUE, Faces = FALSE))
dev.copy(png, "venn_diagram.png")
dev.off()


y$group[y$Sporadic > 2 & y$MEN1 > 2 & y$VHL > 2] <- 2 #All hyper
y$group[y$Sporadic < -2 & y$MEN1 < -2 & y$VHL < -2] <- 3 # All hypo
y$group[y$Sporadic > -2 & y$MEN1 > -2 & y$VHL < -2] <- 4 # Only VHL hypo
y$group[y$Sporadic > -2 & y$MEN1 < -2 & y$VHL > -2] <- 5 # Only MEN1 hypo
y$group[y$Sporadic < -2 & y$MEN1 > -2 & y$VHL > -2] <- 6 # Only Sporadic hypo
y$group[y$Sporadic < 2 & y$MEN1 < 2 & y$VHL > 2] <- 7 # Only VHL hyper
y$group[y$Sporadic < 2 & y$MEN1 > 2 & y$VHL < 2] <- 8 # Only MEN1 hyper
y$group[y$Sporadic > 2 & y$MEN1 < 2 & y$VHL < 2] <- 9 # Only Sporadic hyper


x1 <- melt(y,  measure.vars=c("Sporadic", "MEN1", "VHL"))
x1$value <- as.numeric(as.character(x1$value))
fit <- aov(x1$value ~ x1$variable, data=x1) #for the ANOVA p-value
summary.aov(fit)
colnames(x1) <- c("C", "", "", "", "", "", "", "Group", "M value")
xx <- TukeyHSD(fit)
plot(xx)


#j1 <- droplevels(x1[x1$value > 3 | x1$value < -3, ])
interaction.plot(j1$cg, j1$variable, j1$value)

#boxplot(x1$value[x1$value < -2] ~ x1$variable[x1$value < -2], main="M value per NFPanNET group (mean is black dot)", xlab="Group", ylab="M Value", col=rainbow(3))



p <- ggplot(x1, aes(x1$variable, x1$value)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.001)
p
q <- ggplot(x1, aes(x1$variable, x1$value)) + 
  geom_violin(scale = "area", stat = "ydensity")
q


group <- pData(gset.funnorm)$Sample_Group #<- c(rep("NET", 29), rep("NI", 4))
dmp <- dmpFinder(M, pheno = group , type = "categorical", qCutoff = 0.05, shrinkVar = FALSE)
dmps <- rownames(dmp)
dmp$Name <- rownames(dmp)

# subsetting for large heatmap # used before, now all DMPs and used for heatmap
dmp_new <- dmp #[dmp$qval < 0.001,]# & (dmp$intercept >5 | dmp$intercept < (-5)) ,]
dmps <- rownames(dmp_new)
M_new <- as.data.frame(M[dmps,])

#Select rows for clustering based on their variances (from https://bioconductor.org/packages/release/bioc/vignettes/ConsensusClusterPlus/inst/doc/ConsensusClusterPlus.pdf)

#Calculating variance (mean absolute deviation of M values)
mad_M=apply(M,1,mad)
MM=M[rev(order(mad_M))[1:5000],]

#Normalizing matrix for clustering:
MM = sweep(MM,1, apply(MM,1,median,na.rm=T))


#for(i in 1:nrow(M_new)){
#  M_new$ratio_sq[i] <- (mean(as.numeric(M_new[i,1:29])) - mean(as.numeric(M_new[i,30:33])))^2
#}

MM <- M_new #[M_new$ratio_sq > 49,]
MM[,34] <- NULL


#####################################
####### Added 7/26/2019  START PCA+heatmap - working
#####################################
#####################################

#Principal components plot shows additional but rough clustering of samples

#calculating variances per gene for all samples

#Total variance
beta <- beta_clean #after removing outliers
#beta <- beta_samples #for MEN1 - without removing outliers
#rv <- rowVars(beta)

#Between group variance - using ANOVA F-test
#syndromes

condition <- targets$Sample_Group
condition <- targets$Site
condition <- condition_site_MEN1 <- c("DNET", "PNET", "PNET", "DNET", "DNET", "PNET", "PNET",  "DNET", "DNET", "PNET", "DNET", "DNET", "DNET", "PNET", "DNET", "DNET", "PNET", "PNET", "DNET", "PNET", "PNET", "PNET", "DNET", "PNET", "PNET", "PNET", "PNET", "PNET", "PNET", "GNET", "PNET", "PNET", "GNET", "PNET", "PNET", "DNET", "PNET", "DNET", "DNET", "PNET", "PNET", "DNET")
z <- rowFtests(beta, as.factor(condition), var.equal = TRUE)
select <- order(z$statistic, decreasing=T)[seq_len(min(2200,length(z$statistic)))]

#select_sample_group <- select
#select_site <- select
#select_men1 <- select
#saveRDS(select_sample_group, "select_sample_group.rds")
#saveRDS(select_site, "select_site.rds")
#saveRDS(select_men1, "select_men1.rds")

# calculates PCA based on genes with the highest variances, from the database after log transform and variance stabilization (vsd)
pc <- prcomp(t(beta[select,]))

# set condition
color_syndrome <-  c("Grey", "Black", "Blue")
color_site <-  c("Red", "Pink", "Green", "Purple", "Yellow")

scores <- data.frame(pc$x, condition)

dev.off()

(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
  + geom_point(size = 4)
  + stat_ellipse(level = 0.95)
  + ggtitle("Principal Components Analysis by Syndrome - MEN1")
  #+ scale_colour_brewer(name = " ", palette = "Set1")
  + scale_color_manual(values = color_site)
  + theme(
    plot.title = element_text(face = 'bold'),
    legend.position = c(.9,.85),
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_blank(),
    axis.text.y = element_text(colour = "Black"),
    axis.text.x = element_text(colour = "Black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(color = 'black',fill = NA) 
  ))

ggsave(pcaplot,file="PCA_plot_site_MEN1.jpg")

# scatter plot of rlog transformations between Sample conditions
# nice way to compare control and experimental samples
#head(assay(rld))
#plot(log2(1+counts(dds,normalized=T)[,1:2]),col='black',pch=20,cex=0.3, main='Log2 transformed')
#plot(assay(rld)[,c(1:3,4:6)],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
#plot(assay(rld)[,2:4],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
#plot(assay(rld)[,5:6],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
#plot(assay(rld)[,1:6],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")

###########################
##### Heatmap - by syndrome
###########################
####
#### Added to remove outliers from heamap
####

#breaks for the core of the distribution
breaks=seq(-4, 3, by=0.2) #41 values
#now add outliers
breaks=append(breaks, 10)
breaks=append(breaks, -10, 0)
#create colour panel with length(breaks)-1 colours
mycol <- colorpanel(n=length(breaks)-1,low="blue",high="yellow")

####
#### Added to remove outliers from heamap
####

my_palette <- colorRampPalette(c("blue",'yellow'))(n=4000)

x <- beta[select_sample_group[1:2000],]
dev.off()

#Col labels colors
condition <- targets$Sample_Group
cc_syndrome <- gsub("Sporadic", "Black", condition)
cc_syndrome <- gsub("MEN1", "Grey", cc_syndrome)
cc_syndrome <- gsub("VHL", "Blue", cc_syndrome)

heatmap.2(x, col=mycol,
          scale="row", 
          breaks=breaks,
          key=T, 
          keysize=1, 
          symkey=F,
          density.info="none", 
          trace="none",
          #cexCol=0.5,  
          labRow=F, #rownames(x), #rownames(x), 
          margins = c(0,10),
          main="2000 Most Variably Methylated CpGs",
          distfun  = function(x) dist(x, method="euclidean"),
          hclustfun= function(x) hclust(x, method="complete"),
          ColSideColors=cc_syndrome)
par(lend = 1)
legend("topright", 
       legend = c("Sporadic", "MEN1", "VHL"),
       col = c("Black", "Grey", "Blue"),
       lwd = 10,
       lty = 1,
       cex=0.6)

###################
### Heatmap by site
###################

condition <- targets$Site
cc_site <- gsub("GNET", "Pink", condition)
cc_site <- gsub("DNET", "Red", cc_site)
cc_site <- gsub("SINET", "Purple", cc_site)
cc_site <- gsub("Unknown primary", "Yellow", cc_site)
cc_site <- gsub("PNET", "Green", cc_site)

####
#### Added to remove outliers from heamap
####

#breaks for the core of the distribution
breaks=seq(-4, 3, by=0.1)
#now add outliers
breaks=append(breaks, 5)
breaks=append(breaks, -5, 0)
#create colour panel with length(breaks)-1 colours
mycol <- colorpanel(n=length(breaks)-1,low="blue",high="yellow")

####
#### Added to remove outliers from heamap
####

cc <- cc_site
dev.off()
x <- beta[select_site[1:1000],]

heatmap.2(x, col=mycol,
          scale="row", 
          key=T, 
          breaks=breaks,
          keysize=1, 
          symkey=T,
          density.info="none", 
          trace="none",
          #cexCol=0.5,  
          labRow=F, #rownames(x), #rownames(x), 
          margins = c(0,15),
          #main="1000 Most Variably Methylated CpGs",
          distfun  = function(x) dist(x, method="euclidean"),
          hclustfun= function(x) hclust(x, method="complete"),
          ColSideColors=cc_site)
par(lend = 1)
legend("topright", 
       legend = c("GNET", "DNET", "SINET", "PNET", "Unknown primary"),
       col = c("Pink", "Red", "Purple", "Green", "Yellow"),
       lwd = 10,
       lty = 1,
       cex=0.6)

#############
### MEN1 only
#############

x <- beta[select_men1[1:1000],]
dev.off()
condition <- condition_site_MEN1
cc_men1 <- gsub("GNET", "Pink", condition)
cc_men1 <- gsub("DNET", "Red", cc_men1)
cc_men1 <- gsub("SINET", "Purple", cc_men1)
cc_men1 <- gsub("Unknown primary", "Yellow", cc_men1)
cc_men1 <- gsub("PNET", "Green", cc_men1)

#breaks for the core of the distribution
breaks=seq(-4, 3, by=0.2)
#now add outliers
breaks=append(breaks, 5)
breaks=append(breaks, -5, 0)
#create colour panel with length(breaks)-1 colours
mycol <- colorpanel(n=length(breaks)-1,low="blue",high="yellow")

heatmap.2(x, col=mycol,
          breaks=breaks,
          scale="row", 
          key=T, 
          keysize=1, 
          symkey=T,
          density.info="none", 
          trace="none",
          #cexCol=0.5,  
          labRow=F, #rownames(x), #rownames(x), 
          margins = c(0,15),
          main="1000 Most Variably Methylated CpGs - MEN1",
          distfun  = function(x) dist(x, method="euclidean"),
          hclustfun= function(x) hclust(x, method="complete"),
          ColSideColors=cc_men1)
par(lend = 1)
legend("topright", 
       legend = c("GNET", "DNET", "PNET"),
       col = c("Pink", "Red", "Green"),
       lwd = 10,
       lty = 1,
       cex=0.6)

## Heatmap with multiple rows
x <- beta[select_site[1:1500],] #select for site comparisons

syn <- c("Sporadic", "MEN1", "VHL")
site <- c("GNET", "DNET", "SINET", "PNET", "Unknown primary")

dd <- cbind(cc_site, cc_syndrome)

heatmap.plus(x, col=mycol,
          scale="row", 
          key=T, 
          breaks=breaks,
          keysize=1, 
          symkey=T,
          density.info="none", 
          trace="none",
          #cexCol=0.5,  
          labRow=F, #rownames(x), #rownames(x), 
          margins = c(0,15),
          main="1500 Most Variably Methylated CpGs",
          distfun  = function(x) dist(x, method="euclidean"),
          hclustfun= function(x) hclust(x, method="complete"),
          ColSideColors=dd)
par(lend = 1)
legend("topright", 
       legend = c("GNET", "DNET", "SINET", "PNET", "Unknown primary"),
       col = c("Pink", "Red", "Purple", "Green", "Yellow"),
       lwd = 10,
       lty = 1,
       cex=0.6)

legend("right", 
       legend = c("Sporadic", "MEN1", "VHL"),
       col = c("Black", "Grey", "Blue"),
       lwd = 10,
       lty = 1,
       cex=0.6)

dev.copy(png, paste0(outputPrefix, "-HEATMAP.png"))
dev.off()

###MEN1 only
condition <- targets$Site
cc_site <- gsub("GNET", "Pink", condition)
cc_site <- gsub("DNET", "Red", cc_site)
cc_site <- gsub("SINET", "Purple", cc_site)
cc_site <- gsub("Unknown primary", "Yellow", cc_site)
cc_site <- gsub("PNET", "Green", cc_site)

cc <- cc_site
heatmap.2(x, col=my_palette,
          scale="row", 
          key=T, 
          keysize=1, 
          symkey=T,
          density.info="none", 
          trace="none",
          #cexCol=0.5,  
          labRow=F, #rownames(x), #rownames(x), 
          margins = c(0,15),
          main="1500 Most Variably Methylated CpGs",
          distfun  = function(x) dist(x, method="manhattan"),
          hclustfun= function(x) hclust(x, method="complete"),
          ColSideColors=cc_site)
par(lend = 1)
legend("topright", 
       legend = c("GNET", "DNET", "SINET", "PNET", "Unknown primary"),
       col = c("Pink", "Red", "Purple", "Green", "Yellow"),
       lwd = 10,
       lty = 1,
       cex=0.6)

######################################
####   Added 7/26/2019  END
######################################

########
### DMRS
########

#Change pheno for each...
pheno_sporadic <- c(rep("Sporadic", 44), rep("Normal", 4))
pheno_men1 <- c(rep("MEN1", 42), rep("Normal", 4))
pheno_vhl <- c(rep("VHL", 10), rep("Normal", 4))

dmr_change <-  as.list(c("Sporadic_dmr_hyper", "Sporadic_dmr_hypo", "MEN1_dmr_hyper", "MEN1_dmr_hypo", "VHL_dmr_hyper", "VHL_dmr_hypo"))

bumphunter:::foreachCleanup()
registerDoParallel(cores = 2)

for(j in dmr_change){
  
  M <- get(j)
  
  if(grepl("Sporadic", j)){
    pheno <- pheno_sporadic
    x <- as.matrix(beta_samples[ , 6:18])}
  if(grepl("MEN1", j)){
    pheno <- pheno_men1
    x <- as.matrix(beta_samples[ , 6:19])}
  if(grepl("VHL", j)){
    pheno <- pheno_vhl
    x <- as.matrix(beta_samples[ , 6:19])}
  
  designMatrix <- model.matrix(~ pheno)  
  
  cl <- clusterMaker(M$chr, M$pos, maxGap = 300)
  chr <- M$chr
  pos <- M$pos
  
  #define DMRs
  dmrs <- bumphunter(x, design = designMatrix, pos=pos, chr=chr, cl, cutoff = 0.55,  B=500, type="M", nullMethod = "bootstrap", annotatedTranscripts = TRUE) 
  dmrs$table <- dmrs$table[order(dmrs$table$fwer),]
  saveRDS(dmrs$table, paste(j, "500_dmr.rds", sep="")) #### delete after dmrs are produced
  #plot(dmrs$table$p.value, -log10(dmrs$table$value), main="Volcano plot", xlab="Avg. methylation difference", ylab="-log10 p-value",xlim=c(-.5,.5))
}






#####################
### APC gene analysis
#####################

APC_cg <- c("cg17655851", "cg00577935", "cg08571859", "cg14511739", "cg22035501", "cg11613015", "cg14479889", "cg16970232", "cg03667968", "cg20311501", "cg23938220", "cg02511809", "cg12534150")

beta_APC <- beta_samples[rownames(beta_samples) %in% APC_cg,]

#APC methylated by syndrome
colnames(beta_APC) <- targets$Sample_Group

dd<- melt(beta_APC, value.name = "value")
ddfig <- ggplot(dd, aes(x=factor(Var2), y=value, color=Var2, palette = "jco")) + 
  geom_boxplot() + 
  ylab("Beta values") +
  xlab("Genetic predisposition") +
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
#  stat_compare_means(method = "anova") +
#  stat_compare_means(label = "p.signif", method = "t.test",
#                     ref.group = "Sporadic", hide.ns = TRUE) +
  theme_classic() +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 16),
      axis.title.x = element_text(size = 18), 
      axis.title.y = element_text(size = 18), 
      axis.text = element_text(size = 16),
      legend.position = c(.83, .8)
      ) 

ddfig

#APC methylated by NET site
colnames(beta_APC) <- targets$Site

dd<- melt(beta_APC, value.name = "value")
ddfig <- ggplot(dd, aes(x=factor(Var2), y=value, color = Var2)) + 
  geom_boxplot() + 
  ylab("Beta values") +
  xlab("NET site") +
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18), 
    axis.text = element_text(size = 16),
     ) 
ddfig

beta_APC_MEN1 <- beta_APC[, c(2,  3,  4,  6, 10, 11, 12, 13, 14, 19, 20, 21, 22, 26, 28, 29, 30, 34, 36, 37, 38, 42, 44, 45, 46, 50, 51, 53, 59, 61, 67, 68, 69, 74, 76, 77, 82, 85, 89, 90, 91, 93)]

dd<- melt(beta_APC_MEN1, value.name = "value")
ddfig <- ggplot(dd, aes(x=factor(Var2), y=value, color=Var2, palette = "jco")) + 
  geom_boxplot() + 
  ylab("Beta values") +
  xlab("Genetic predisposition") +
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
#  stat_compare_means(method = "anova") +
#  stat_compare_means(label = "p.signif", method = "t.test",
#                     ref.group = "GNET", hide.ns = TRUE) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18), 
    axis.text = element_text(size = 16),
    legend.position = "none"
  ) 

ddfig


#APC methylated by gender
gender <- c("Female", "Male", "Male", "Female", "Male", "Male", "Female", "Female", "Male", "Male", "Male", "Female", "Female", "Male","Female", "Male", "Female", "Female", "Male", "Male", "Female", "Male", "Female", "Male", "Male", "Male", "Male", "Male", "Female", "Female", "Male", "Male", "Male", "Male", "Male", "Male", "Male", "Female", "Male", "Male", "Male", "Male", "Male", "Male", "Male", "Female", "Female", "Male", "Female", "Male", "Male", "Male", "Male", "Female", "Female", "Female", "Female", "Male", "Male", "Male", "Male", "Female", "Female", "Male", "Male", "Male", "Male", "Female", "Male", "Female", "Female", "Female", "Male", "Male", "Female", "Female", "Female", "Male", "Female", "Male", "Female", "Male", "Female", "Male", "Female", "Male", "Female", "Male", "Male", "Male", "Female", "Male", "Female", "Male", "Male", "Female")

colnames(beta_APC) <- gender

dd<- melt(beta_APC, value.name = "value")
ddfig <- ggplot(dd, aes(x=factor(Var2), y=value, color=Var2)) + 
  geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_classic() +
  ylab("Beta values") +
  xlab("Gender") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18), 
    axis.text = element_text(size = 16),
  ) 
ddfig


beta_APC_MEN1 <- beta_APC[, c(2,  3,  4,  6, 10, 11, 12, 13, 14, 19, 20, 21, 22, 26, 28, 29, 30, 34, 36, 37, 38, 42, 44, 45, 46, 50, 51, 53, 59, 61, 67, 68, 69, 74, 76, 77, 82, 85, 89, 90, 91, 93)]

dd<- melt(beta_APC_MEN1, value.name = "value")
ddfig <- ggplot(dd, aes(x=factor(Var2), y=value, color=Var2, palette = "jco")) + 
  geom_boxplot() + 
  ylab("Beta values") +
  xlab("Genetic predisposition") +
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  stat_compare_means(method = "t.test") +
  #  stat_compare_means(label = "p.signif", method = "t.test",
  #                     ref.group = "Sporadic", hide.ns = TRUE) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18), 
    axis.text = element_text(size = 16),
    legend.position = c(.83, .8)
  ) 

ddfig


















  




##looping figures for all variales
data_var <- as.list(c("M_VHLdmp_hypo", "M_VHLdmp_hyper", "M_MEN1dmp_hypo", "M_MEN1dmp_hyper", "M_Sporadic_dmp_hypo", "M_Sporadic_dmp_hyper"))

#G_hypo", "G_hyper", "LN_hypo", "LN_hyper"))#, "vhl_hypo", "vhl_hyper", "sporadic_hypo", "sporadic_hyper", "men1_hypo", "men1_hyper" ))

for(i in data_var){
  data1 <- get(i)
    data1$Relation_to_Island <- replace(data1$Relation_to_Island, data1$Relation_to_Island == "OpenSea", "Open Sea") 
  data1$Relation_to_Island <- replace(data1$Relation_to_Island, data1$Relation_to_Island == "S_Shore", "Shore") 
  data1$Relation_to_Island <- replace(data1$Relation_to_Island, data1$Relation_to_Island == "N_Shore", "Shore") 
  data1$Relation_to_Island <- replace(data1$Relation_to_Island, data1$Relation_to_Island == "N_Shelf", "Shelf") 
  data1$Relation_to_Island <- replace(data1$Relation_to_Island, data1$Relation_to_Island == "Island", "Island") 
  data1$Relation_to_Island <- replace(data1$Relation_to_Island, data1$Relation_to_Island == "S_Shelf", "Shelf") 

#Summarizing n of methylation positions 
a <- dcast(data1, data1$qval ~ data1$Relation_to_Island )

sum_island_type <- c(paste(sum(a[,2])), paste(sum(a[,3])), paste(sum(a[,4])), paste(sum(a[,5])))

###Pie chart

#Font specification
t <- list(family = "sans serif", size = 20, color = 'black')
p <- plot_ly(a, labels = ~colnames(a[,2:5]), values = ~sum_island_type, type = 'pie',
             textposition = 'inside',
             textinfo = 'label+value',
             insidetextfont = list(color = '#FFFFFF'),
             marker = list(colors = colors,
                           line = list(color = '#FFFFFF', width = 1.5)),
             #The 'pull' attribute can also be used to create space between the sectors
             showlegend = FALSE) %>%
  layout(title = paste(i), font=t, xaxis = list(showgrid = FALSE, zeroline = FALSE,
                                                showticklabels = FALSE), 
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

print(p)
export(p = last_plot(), file = paste(i, "_pie.pdf"), selenium = NULL)
}

#delta M value charts
for(i in data_var){
  data2 <- get(i)

Relation_to_Island <- as.factor(data2$Relation_to_Island)

Island_levels <- levels(Relation_to_Island)
#labs <- c("Island", "Shelf", "Shore", "Open Sea", "Shelf", "Shore") 
#df1 <- melt(data.frame(Island_levels, labs))

data2$Relation_to_Island <- factor(data2$Relation_to_Island, levels=c("Island",  "N_Shelf", "N_Shore", "OpenSea", "S_Shelf", "S_Shore"), labels=c("Island", "Shelf", "Shore", "Open Sea", "Shelf ", "Shore "))

data2$Relation_to_Island_n <- factor(data2$Relation_to_Island, levels=c("Shelf", "Shore", "Island", "Open Sea", "Shore ", "Shelf "))


#Relation_to_Island <- droplevels(Relation_to_Island)

a <- ggplot(subset(data2, !(Relation_to_Island == "Open Sea")), 
            aes(x = pos, y = dM)) + 
  geom_line(aes(color =  Relation_to_Island_n)) + 
  facet_grid(. ~ Relation_to_Island_n) + 
  theme(legend.position = "none",
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")) + 
  stat_smooth(method = "loess") + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank())

pdf(paste(i, "dM.pdf", sep=""),5,1.33)
print(a)
dev.off()
}

#Figures sporadic vs. normal
##pie charts: https://plot.ly/r/pie-charts/
##plot using: xxxx + smoothing, from here: http://felixfan.github.io/stacking-plots-same-x/
##plot dm/island type

######
# DMRs
######
# For DMRs use bumphunter - First increase the threshold 
# (start with cutoff=0.2) to have low bumps count, and then 
# repeat with a high permutation count ( start with B=0, and 
# move to 1000)
# Finding differentially methylated regions (DMRs)




#make clusters for bumphunter
setwd("C:/Users/tirosha/Desktop/Work/Projects/Sanjit methylation/sanjit_methylation/image data")
M_VHLdmp_hypo <- readRDS("VHLhypo.rds")
M_VHLdmp_hyper <- readRDS("VHLhyper.rds")
M_MEN1dmp_hypo <- readRDS("MEN1hypo.rds")
M_MEN1dmp_hyper <- readRDS("MEN1hyper.rds")
M_Sporadic_dmp_hypo <- readRDS("Sporadic_hypo.rds")
M_Sporadic_dmp_hyper <- readRDS("Sporadic_hyper.rds")


#Change pheno for each...
#pheno <- c(rep("Sporadic", 9), rep("Normal", 4))
pheno_men1 <- c(rep("MEN1", 10), rep("Normal", 4))
pheno_sporadic <- c(rep("Sporadic", 9), rep("Normal", 4))
pheno_vhl <- c(rep("VHL", 10), rep("Normal", 4))

dmr_change <-  as.list(c("M_Sporadic_dmp_hyper", "M_Sporadic_dmp_hypo", "M_MEN1dmp_hyper", "M_MEN1dmp_hypo", "M_VHLdmp_hyper", "M_VHLdmp_hypo"))

#gset.funnorm <- readRDS("gsetFunnorm_MEN1.rds")
#gset.funnorm <- readRDS("gsetFunnorm_VHL.rds")

bumphunter:::foreachCleanup()
registerDoParallel(cores = 2)

for(j in dmr_change){
  
  M <- get(j)
  
  if(grepl("Sporadic", j)){
    pheno <- pheno_sporadic
    x <- as.matrix(M[ , 6:18])}
  if(grepl("MEN1", j)){
    pheno <- pheno_men1
    x <- as.matrix(M[ , 6:19])}
  if(grepl("VHL", j)){
    pheno <- pheno_vhl
    x <- as.matrix(M[ , 6:19])}
  
    designMatrix <- model.matrix(~ pheno)  

  cl <- clusterMaker(M$chr, M$pos, maxGap = 300)
  chr <- M$chr
  pos <- M$pos
  
  #define DMRs
  dmrs <- bumphunter(x, design = designMatrix, pos=pos, chr=chr, cl, cutoff = 0.55,  B=500, type="M", nullMethod = "bootstrap", annotatedTranscripts = TRUE) 
  dmrs$table <- dmrs$table[order(dmrs$table$fwer),]
  saveRDS(dmrs$table, paste(j, "500_dmr.rds", sep="")) #### delete after dmrs are produced
  #plot(dmrs$table$p.value, -log10(dmrs$table$value), main="Volcano plot", xlab="Avg. methylation difference", ylab="-log10 p-value",xlim=c(-.5,.5))
}


#Next - visualization 
# https://support.bioconductor.org/p/79684/
setwd("C:/Users/tirosha/Desktop/Work/Projects/Sanjit methylation/sanjit_methylation/image data")

Sporadic_hyper_dmrs <- readRDS("M_Sporadic_dmp_hyper500_dmr.rds")
Sporadic_hypo_dmrs <- readRDS("M_Sporadic_dmp_hypo500_dmr.rds")
MEN1_hyper_dmrs <- readRDS("M_MEN1dmp_hyper500_dmr.rds")
MEN1_hypo_dmrs <- readRDS("M_MEN1dmp_hypo500_dmr.rds")
VHL_hyper_dmrs <- readRDS("M_VHLdmp_hyper500_dmr.rds")
VHL_hypo_dmrs <- readRDS("M_VHLdmp_hypo500_dmr.rds")

dmrs <- Sporadic_hyper_dmrs
dmrs <- Sporadic_hypo_dmrs 
dmrs <- MEN1_hyper_dmrs

gset.funnorm <- readRDS("gsetFunnorm_Sporadic.rds")
groups <- c(rep("Sporadic", 9), rep("Normal", 4)) #c(targets$Sample_Group)

gset.funnorm <- readRDS("gsetFunnorm_VHL.rds")
groups <- c(rep("VHL", 10), rep("Normal", 4)) #c(targets$Sample_Group)

gset.funnorm <- readRDS("gsetFunnorm_MEN1.rds")
groups <- c(rep("MEN1", 10), rep("Normal", 4)) #c(targets$Sample_Group)
#gene/line/
#CDCA7L, MEN1, 12750, -500 and +300 in the final plot range
#RBM47, MEN1, 10551, -3000, +30
#SFRP5, sporadic, 363

    genome <- "hg19"
    dmr <- dmrs[rownames(dmrs)==363,] #was dmrs$table[i,]
    chrom <- dmr$chr
    start <- dmr$start 
    end <- dmr$end
    minbase <- start - 0.25 * (end - start)
    maxbase <- end + 0.25 * (end - start)
    pal <- c("#E41A1C", "#377EB8")
    
# Start building the tracks
    iTrack <- IdeogramTrack(genome = genome, 
                            chromosome = chrom, 
                            name = "")
    gTrack <- GenomeAxisTrack(col = "black", 
                              cex = 1, 
                              name = "", 
                              fontcolor = "black")
  
    #grTrack <- GeneRegionTrack(dmr, 
    #                           genome = genome, 
    #                       chromosome = chrom, 
    #                       rstart=start, rends=end,
    #                       name = "Gene Model",
    #                       transcriptAnnotation = "symbol")


    rTrack <- UcscTrack(genome = genome, 
                        chromosome = chrom, 
                        track = "NCBI RefSeq", 
                        from = minbase-1000,
                        to = maxbase+1000, 
                        trackType = "GeneRegionTrack",
                        rstarts = "exonStarts", 
                        rends = "exonEnds", 
                        gene = "name",
                        symbol = "name2", 
                        transcript = "name",
                        strand = "strand",
                        fill = "darkblue",
                        stacking = "squish", 
                        name = "RefSeq",
                        showId = TRUE, 
                        transcriptAnnotation = "symbol",
                        geneSymbol = TRUE)
    
# methylation data track
    gr <- granges(gset.funnorm)
    gr$beta <- getBeta(gset.funnorm)
    methTrack <- DataTrack(range = gr,
                           groups = groups,
                           genome = genome,
                           chromosome = chrom, 
                           ylim = c(-0.05, 1.05), 
                           col = pal,
                           type = c("a","p"), 
                           name = "DNA Meth.\n(beta value)",
                           background.panel = "white", 
                           legend = TRUE, 
                           cex.title = 0.8,
                           cex.axis = 0.8, 
                           cex.legend = 0.8)
    #DMR position data track
      dmrTrack <- AnnotationTrack(start = dmr$start, 
                                 end = dmr$end, 
                                genome = genome, 
                                name = "DMR",
                                chromosom = chrom)
    # Finally, plot the tracks
    tracks <- list(iTrack, gTrack, methTrack, dmrTrack, rTrack) 
    sizes <- c(2, 2, 5, 2, 3) #set up the relative sizes of the tracks
    
    pdf("SFRP5_CIMP_MEN1.pdf", 6, 10) #paste(1, "_", "DMRs_", i, ".pdf", sep=""), 5, 8)
    tiff("SFRP5_CIMP_MEN1.tiff", 1750, 3000, res = 500)
    plotTracks(tracks, 
               from = minbase - 1500, 
               to = maxbase + 30, 
               showTitle = TRUE, 
               add53 = TRUE,
               add35 = TRUE,
               grid = TRUE, 
               lty.grid = 3, 
               sizes = sizes,
               length(tracks))
    dev.off()




VHL_hypo_dmrs <- readRDS("M_VHLdmp_hypo500_dmr.rds")
VHL_hyper_dmrs <- readRDS("M_VHLdmp_hyper500_dmr.rds")
MEN1_hypo_dmrs <- readRDS("M_MEN1dmp_hypo500_dmr.rds")
MEN1_hyper_dmrs <- readRDS("M_MEN1dmp_hyper500_dmr.rds")
Sporadic_hypo_dmrs <- readRDS("M_Sporadic_dmp_hypo500_dmr.rds")
Sporadic_hyper_dmrs <- readRDS("M_Sporadic_dmp_hyper500_dmr.rds")

dmrs_datas <- as.list(c("VHL_hypo_dmrs", "VHL_hyper_dmrs", "MEN1_hypo_dmrs", "MEN1_hyper_dmrs", "Sporadic_hypo_dmrs", "Sporadic_hyper_dmrs"  ))

for(aa in dmrs_datas){
  dd <- get(aa)
  print(paste(aa, nrow(dd[dd$fwerArea < 0.05,])))
  write.xlsx(dd[dd$fwerArea < 0.05,], paste(aa, ".xlsx", sep=""))
}

#CIMP based on DMRs - used in final version
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
subject <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene, annotation="org.Hs.eg.db")

# annotatr package
for(p in dmrs_datas){
  genes_dmrs <- get(p)
  genes_dmrs <- genes_dmrs[(genes_dmrs$fwerArea < 0.05), ]
  gene_list <- paste(p, "_gene_list", sep="") 

    x <- matchGenes(genes_dmrs, subject, type = "any", promoterDist = 2500, skipExons = FALSE, verbose = TRUE)
  write.xlsx(x, paste(gene_list, ".xlsx", sep=""))
    saveRDS(x, paste(gene_list, ".rds", sep=""))

## add here annotation for ENCODE#######    
}




genes_dmrs_vhl_hypo <- readRDS("VHL_hypo_dmrs_gene_list.rds")
genes_dmrs_vhl_hyper <- readRDS("VHL_hyper_dmrs_gene_list.rds")
genes_dmrs_MEN1_hypo <- readRDS("MEN1_hypo_dmrs_gene_list.rds")
genes_dmrs_MEN1_hyper <- readRDS("MEN1_hyper_dmrs_gene_list.rds")
genes_dmrs_Sporadic_hypo <- readRDS("Sporadic_hyper_dmrs_gene_list.rds")
genes_dmrs_Sporadic_hyper <- readRDS("Sporadic_hypo_dmrs_gene_list.rds")

hyper <- list(VHL = genes_dmrs_vhl_hyper[,1], MEN1 = genes_dmrs_MEN1_hyper[,1], Sporadic = genes_dmrs_Sporadic_hyper[,1])

hypo <- list(VHL = genes_dmrs_vhl_hypo[,1], MEN1 = genes_dmrs_MEN1_hypo[,1], Sporadic = genes_dmrs_Sporadic_hypo[,1])

v1 <- Venn(hyper)
v1
v2 <- Venn(hypo)
v2

dev.off()
plot(v1, doWeights = FALSE, type = "circles", show = list(SetLabels = TRUE, Faces = FALSE))
dev.copy(png, "venn_diagram_hyper.png")
dev.off()
pdf("venn_diagram_hyper.pdf", 5, 5)
plot(v1, doWeights = FALSE, type = "circles", show = list(SetLabels = TRUE, Faces = FALSE))
dev.off()

pdf("venn_diagram_hypo.pdf", 5, 5)
plot(v2, doWeights = FALSE, type = "circles", show = list(SetLabels = TRUE, Faces = FALSE))
dev.off()
plot(v2, doWeights = FALSE, type = "circles", show = list(SetLabels = TRUE, Faces = FALSE))
dev.copy(png, "venn_diagram_hypo.png")
dev.off()

setwd("C:/Users/tirosha/Desktop/Work/Projects/Sanjit methylation/sanjit_methylation")
geneexp <- read.csv("genes_exp.csv")
ge_sporadic_hyper <- geneexp[geneexp$group == "sporadic" & geneexp$meth.status == "hyper",]
ge_sporadic_hypo <- geneexp[geneexp$group == "sporadic" & geneexp$meth.status == "hypo",]
ge_men1_hyper <- geneexp[geneexp$group == "men1" & geneexp$meth.status == "hyper",]
ge_men1_hypo <- geneexp[geneexp$group == "men1" & geneexp$meth.status == "hypo",]
ge_vhl_hyper <- geneexp[geneexp$group == "vhl" & geneexp$meth.status == "hyper",]
ge_vhl_hypo <- geneexp[geneexp$group == "vhl" & geneexp$meth.status == "hypo",]

dmrs_genes <- c("ge_sporadic_hyper", "ge_sporadic_hypo", "ge_men1_hyper",  "ge_men1_hypo",  "ge_vhl_hyper", "ge_vhl_hypo")

for(k in dmrs_genes){
  geneexp <- get(k)
  print(k)
  geneexp <- geneexp[order(log2(geneexp$Fold.Change)),]
  fig_ge <- ggplot(geneexp, aes(x = log2(geneexp$Fold.Change), 
                                y = geneexp$gene, 
                                size = geneexp$P.value, 
                                fill = geneexp$Regulation,
                                reorder(geneexp$gene, geneexp$Fold.Change))) +  
    geom_point(shape = 21) + 
    #ggtitle("Gene expression in gene with aberrated methylation") + 
    labs(x = "Log2 Fold Change", 
       y = "Gene Name",
       size = "P-value",
       fill = "Gene Regulation") + 
    scale_x_continuous(breaks=c(1,2,3,4,5,6), limits=c(0, 6)) +
  scale_size(range = c(2, 10)) +
    theme(
      #legend.position="bottom", #legend.direction="horizontal",
      ##legend.box = "horizontal",
      legend.key.size = unit(1, "cm"),
      legend.text = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.title = element_text(colour="black", size = 12),
      axis.text.x=element_text(colour="black", size = 12),
      axis.text.y=element_text(colour="black", size = 12),
      axis.line = element_line(size=1, colour = "black"),
      axis.ticks.length = unit(1, "mm"))
  tiff(paste(k, ".tiff", sep=""), 2200, 4000, res = 500)
  print(fig_ge)
  dev.off()
}




# annotatr (did not use this part eventually, everything is in the minfi island data)
#annots = c('hg19_cpgs', 'hg19_basicgenes', 'hg19_genes_intergenic',
#           'hg19_genes_intronexonboundaries')
#annotations = build_annotations(genome = 'hg19', annotations = annots)
#genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene,annotation="org.Hs.eg.db")
#dm_regions <- VHL_hyper_dmrs
#dm_annotated = annotate_regions(
#  regions = dm_regions,
#  annotations = annotations,
#  ignore.strand = TRUE,
#  quiet = FALSE)

##################################################
# Gene annotation from DMPs - for pathway analysis
##################################################

setwd("C:/Users/tirosha/Desktop/Work/Projects/Sanjit methylation/sanjit_methylation/image data")

M_VHLdmp_hypo <- readRDS("VHLhypo.rds")
M_VHLdmp_hyper <- readRDS("VHLhyper.rds")
M_MEN1dmp_hypo <- readRDS("MEN1hypo.rds")
M_MEN1dmp_hyper <- readRDS("MEN1hyper.rds")
M_Sporadic_dmp_hypo <- readRDS("Sporadic_hypo.rds")
M_Sporadic_dmp_hyper <- readRDS("Sporadic_hyper.rds")

DMP_datasets <- as.list(c("M_Sporadic_dmp_hyper", "M_Sporadic_dmp_hypo", "M_VHLdmp_hypo", "M_VHLdmp_hyper", "M_MEN1dmp_hyper",  "M_MEN1dmp_hypo"))

for(g in DMP_datasets){
  aa <- get(g)
  aa <- aa[aa$pval < 0.001,]
  #aa <- aa[aa$dM >3 | aa$dM < (-3), ] #Relation_to_Island != "OpenSea",]
    assign(g, aa)
}

#################
# Annotating DMPs
#################
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
subject <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene,annotation="org.Hs.eg.db", by = "gene", codingOnly = TRUE)

for(p in DMP_datasets){
  CIMP <- get(p)
  CIMP$start <- CIMP$pos
  CIMP$end <- CIMP$pos
  gene_list <- paste(p, "_gene_list", sep="") 
  x <- matchGenes(CIMP, subject, type = "any", promoterDist = 2500, skipExons = TRUE, verbose = FALSE)
  saveRDS(x, paste(gene_list, ".rds", sep=""))
  assign(gene_list, x)
}

#gene_list_VHLhyper <- readRDS("M_VHLdmp_hyper_gene_list.rds")
#gene_list_VHLhypo <- readRDS("M_VHLdmp_hypo_gene_list.rds")
#gene_list_MEN1hyper <- readRDS("M_MEN1dmp_hyper_gene_list.rds") 
#gene_list_MEN1hypo <- readRDS("M_MEN1dmp_hypo_gene_list.rds") 
gene_list_Sporadichyper <- readRDS("M_Sporadic_dmp_hyper_gene_list.rds")
gene_list_Sporadichypo <- readRDS("M_Sporadic_dmp_hypo_gene_list.rds")

gene_list_VHLhyper <- M_VHLdmp_hyper_gene_list
gene_list_VHLhypo <- M_VHLdmp_hypo_gene_list
gene_list_MEN1hyper <- M_MEN1dmp_hyper_gene_list
gene_list_MEN1hypo <- M_MEN1dmp_hypo_gene_list
gene_list_Sporadichypo <- M_Sporadic_dmp_hypo_gene_list
gene_list_Sporadichyper <- M_Sporadic_dmp_hyper_gene_list


gene_list_VHLhyper <- gene_list_VHLhyper[,1]
gene_list_VHLhypo <- gene_list_VHLhypo[,1]
gene_list_MEN1hyper <- gene_list_MEN1hyper[,1]
gene_list_MEN1hypo <- gene_list_MEN1hypo[,1]
gene_list_Sporadichypo <- gene_list_Sporadichypo[,1]
gene_list_Sporadichyper <- gene_list_Sporadichyper[,1]

#gene_list_VHLhyper <- as.data.frame(gene_list_VHLhyper)
#gene_list_MEN1hyper <- as.data.frame(gene_list_MEN1hyper)
gene_list_Sporadic_hyper <- as.data.frame(gene_list_Sporadichyper)
gene_list_Sporadic_hypo <- as.data.frame(gene_list_Sporadichypo)

#write.xlsx(gene_list_VHLhyper, "gene_list_VHLhyper.xlsx")
#write.xlsx(gene_list_MEN1hyper, "gene_list_MEN1hyper.xlsx")
#write.xlsx(gene_list_Sporadic_hyper, "gene_list_Sporadic_hyper.xlsx")

gene_list_VHL <- c(gene_list_VHLhyper, gene_list_VHLhypo)
gene_list_MEN1 <- as.character(c(gene_list_MEN1hyper, gene_list_MEN1hypo))
gene_list_Sporadic <- c(gene_list_Sporadichyper, gene_list_Sporadichypo)

#saveRDS(gene_list_VHL, "gene_list_VHL.rds")
#saveRDS(gene_list_MEN1, "gene_list_MEN1.rds")
#saveRDS(gene_list_Sporadic, "gene_list_Sporadic.rds")


#########
# Pathway
#########

setwd("C:/Users/tirosha/Desktop/Work/Projects/Sanjit methylation/sanjit_methylation/image data")

###Annotation and pathway start here:
#gene_list_VHL <- readRDS("gene_list_VHL.rds")
#gene_list_MEN1 <- readRDS("gene_list_MEN1.rds")
#gene_list_Sporadic <- readRDS("gene_list_Sporadic.rds")

GL <- as.list(c("gene_list_Sporadic","gene_list_VHL", "gene_list_MEN1"))
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
j <- "gene_list_Sporadic"
for(j in GL){
  geneList <- get(j)
  myHugoGeneNames <- as.vector(geneList)
  mapTab <- getBM(attributes = c("hgnc_symbol", "entrezgene"), filters = "hgnc_symbol", values = myHugoGeneNames, mart = ensembl, uniqueRows=TRUE)
  d <- as.data.frame(table(mapTab$entrezgene))
  geneList = d[,1]
  #names(geneList) = as.character(d[,1])
  
  geneList = sort(geneList, decreasing = TRUE)
  de <- geneList
x <- enrichPathway(gene=de, organism = "human", pAdjustMethod = "fdr", readable = FALSE, pvalueCutoff=0.05)
head(as.data.frame(x))
pdf(paste(j, "_pathway_analysis.pdf"),10 ,10)
#tiff(paste(j, "_pathway_analysis.tiff"), 3000, 6000, res = 500)

print(dotplot(x, showCategory = 15, font.size = 12))
dev.off()
}
dotplot(x, showCategory = 15, font.size = 12)

barplot(x, showCategory=8)
DOSE::dotplot(x, showCategory=10)
enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
cnetplot(x, categorySize="pvalue", foldChange=geneList)
#GSEA analysis
y <- gsePathway(geneList, nPerm=1000,
                minGSSize=1, pvalueCutoff=1,
                pAdjustMethod="BH", verbose=FALSE)
res <- as.data.frame(y)
head(res)
enrichMap(y) #not workig
gseaplot(y, geneSetID = "R-HSA-5357801")
viewPathway("E2F mediated regulation of DNA replication", readable=TRUE, foldChange=geneList)
#DOSE for visualization?
library(igraph)

#Cluster concensus analysis
results = ConsensusClusterPlus(MM,
                               maxK=6,
                               reps=50,
                               pItem=0.8,
                               pFeature=1,
                               title="Cluster consensus",
                               clusterAlg="hc",
                               distance="pearson",
                               #seed=1262118388.71279,
                               plot="png")


icl = calcICL(results,
              title="Cluster concencsus",
              plot="png")

cluster <- as.data.frame(results[[4]][["consensusClass"]])

clustering <- cbind(cluster, targets$Sample_Group , targets$Slide, targets$Array )



#Example to show that it now also works with just a single column or single row
mat <- matrix(1:100, byrow=T, nrow=10)
column_annotation <- sample(c("red", "blue", "green"), 10, replace=T)
column_annotation <- as.matrix(column_annotation)
colnames(column_annotation) <- c("Variable X")

row_annotation <- sample(c("red", "blue", "green"), 10, replace=T)
row_annotation <- as.matrix(t(row_annotation))
rownames(row_annotation) <- c("Variable Y")

heatmap.3(matix, ColSideColors=column_annotation)

# Revision 2020 06 07 ####

#saveRDS(targets, "targets.RDS")

targets <- targets[targets$Sample_Group == "MEN1", ] %>% droplevels()
matrix_men1 <- readRDS("beta_MEN1.rds")
beta <- matrix_men1
# MEN1 func PCA + heatmap only
targets <- targets[targets$Functional %in% c("gastrinoma", "NFPNET"),] %>% droplevels()
beta <- beta[, which(targets$Sample_Group == "MEN1" & targets$Functional %in% c("gastrinoma", "NFPNET"))]

condition <- targets$Functional
z <- rowFtests(beta, as.factor(condition), var.equal = TRUE)
select <- order(z$statistic, decreasing=T)[seq_len(min(2200,length(z$statistic)))]

pc <- prcomp(t(beta[select,]))

# set condition
scores <- data.frame(pc$x, condition)
levels(scores$condition) <- c("Gastrinoma", "Non functioning NET")
dev.off()

(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
  + geom_point(size = 4)
  + stat_ellipse(level = 0.90)
  + ggtitle("Principal Components Analysis by Syndrome - MEN1")
  + scale_colour_brewer(name = " ", palette = "Set1")
  #+ scale_color_manual(values = color_site)
  + theme(
    plot.title = element_text(face = 'bold'),
    legend.position = c(.1,.85),
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_blank(),
    axis.text.y = element_text(colour = "Black"),
    axis.text.x = element_text(colour = "Black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(color = 'black',fill = NA) 
  ))

cc_func <- scores$condition
levels(cc_func) <- c("Green", "Blue")  #gastrinoma greem, NFPNET blue
cc_func <- as.character(cc_func)
my_palette <- colorRampPalette(c("blue",'yellow'))(n=400)

x <- beta[select, ]
dev.off()
targets$Site
site <- c("DNET", "PNET")
func <- c("Gastrinoma", "NFPNET")

cc_func <- scores$condition
levels(cc_func) <- c("Green", "Blue")  #gastrinoma greem, NFPNET blue
cc_func <- as.character(cc_func)

cc_site <- targets$Site %>% factor()
levels(cc_site) <- c("Red", "Yellow")  #DNET red, PNET yellow
cc_site <- as.character(cc_site)

dd <- cbind(cc_func, cc_site)

dev.off()
pdf("hm_MEN1_func.pdf", 12,6,6)
  heatmap.plus(x, col=my_palette,
             scale="row", 
             key=T, 
             keysize=1, 
             symkey=T,
             density.info="none", 
             trace="none",
             #cexCol=0.5,  
             labRow=F, #rownames(x), #rownames(x), 
             margins = c(0,15),
             main="",
             distfun  = function(x) dist(x, method="manhattan"),
             hclustfun= function(x) hclust(x, method="complete"),
             ColSideColors=dd)
par(lend = 1)
legend("topright", 
       legend = c("DNET", "PNET"),
       col = c("Red", "Yellow"),
       lwd = 10,
       lty = 1,
       cex=0.7)

legend("right", 
       legend = c("Gastrinoma", "NFPNET"),
       col = c("Green", "Blue"),
       lwd = 10,
       lty = 1,
       cex=0.7)

dev.off()
dev.copy(png, paste0(outputPrefix, "-HEATMAP.png"))

## APC #######
setwd("C:/Users/tiros/Documents/Work/Projects/Methylation/raw data")
baseDir <- ("C:/Users/tiros/Documents/Work/Projects/Methylation/raw data")
targets <- read.metharray.sheet(baseDir, pattern = "Sample_Sheet.csv")
beta_samples <- readRDS("beta_96.rds")
APC_cg <- c("cg17655851", "cg00577935", "cg08571859", "cg14511739", "cg22035501", "cg11613015", "cg14479889", "cg16970232", "cg03667968", "cg20311501", "cg23938220", "cg02511809", "cg12534150")

beta_APC <- beta_samples[APC_cg,]
targets$Sample_Group[which(colMins(beta_APC) > 0.5)]

targets$Grade[which(colMins(beta_APC) > 0.5)] %>% summary.factor()
targets$AnyMet[which(colMins(beta_APC) > 0.5)] %>% summary.factor()

targets$AnyMet[which(colMins(beta_APC) <= 0.5)] %>% summary.factor()
targets$LN[which(colMins(beta_APC) > 0.5)] %>% summary.factor()
targets$Dis[which(colMins(beta_APC) > 0.5)] %>% summary.factor()
targets$LN[which(colMins(beta_APC) <= 0.5)] %>% summary.factor()
targets$Dis[which(colMins(beta_APC) <= 0.5)] %>% summary.factor()

targets$AnyMet[which(targets$AnyMet == "No" & targets$Sample_Group == "Sporadic")] %>% summary.factor()

# Comparison by staging ####
setwd("C:/Users/tiros/Documents/Work/Projects/Methylation/raw data")
baseDir <- ("C:/Users/tiros/Documents/Work/Projects/Methylation/raw data")
targets <- read.metharray.sheet(baseDir, pattern = "Sample_Sheet.csv")
beta <- readRDS("beta_96.rds")
APC_cg <- c("cg17655851", "cg00577935", "cg08571859", "cg14511739", "cg22035501", "cg11613015", "cg14479889", "cg16970232", "cg03667968", "cg20311501", "cg23938220", "cg02511809", "cg12534150")
beta_APC <- beta_samples[APC_cg,]

targets$APC <- rep(22, nrow(targets))
targets$APC[which(colMins(beta_APC) > 0.5)] <- "High"
targets$APC[which(colMins(beta_APC) <= 0.5)] <- "Low"
targets$APC_adv <- ifelse(targets$APC == "High" & targets$AnyMet == "Yes", "HighMet",
                          ifelse(targets$APC == "High" & targets$AnyMet == "No", "HighNoMet",
                                 ifelse(targets$APC == "Low" & targets$AnyMet == "Yes", "LowMet",
                                        ifelse(targets$APC == "Low" & targets$AnyMet == "No", "LowNoMet",NA)))) %>% factor()
            
condition <- targets$APC
z <- rowFtests(beta, as.factor(condition), var.equal = TRUE)
select <- order(z$statistic, decreasing=T)[seq_len(min(2200,length(z$statistic)))]

pc <- prcomp(t(beta[select,]))

# set condition
scores <- data.frame(pc$x, condition)
#levels(scores$condition) <- c("APC hypermethylation", "No APC hypermethylation")
dev.off()

(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
  + geom_point(size = 4)
  + stat_ellipse(level = 0.95)
  + ggtitle("Principal Components Analysis by APC methylation status")
  + scale_colour_brewer(name = " ", palette = "Set1")
  #+ scale_color_manual(values = color_site)
  + theme(
    plot.title = element_text(face = 'bold'),
    legend.position = c(.2,.9),
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_blank(),
    axis.text.y = element_text(colour = "Black"),
    axis.text.x = element_text(colour = "Black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(color = 'black',fill = NA) 
  ))

condition <- targets$AnyMet
z <- rowFtests(beta, as.factor(condition), var.equal = TRUE)
select <- order(z$statistic, decreasing=T)[seq_len(min(2200,length(z$statistic)))]
scores <- data.frame(pc$x, condition)

x <- beta[select, ]
dev.off()

cc_APC <- targets$APC %>% factor()
levels(cc_APC) <- c("Blue", "Red")  #APC high blue, APC low red
cc_APC <- as.character(cc_APC)

cc_adv <- targets$AnyMet %>% factor()
levels(cc_adv) <- c("Green", "Yellow")  #advanced yellow, localized green
cc_adv <- as.character(cc_adv)

dd <- cbind(cc_adv, cc_APC)
my_palette <- colorRampPalette(c("blue",'yellow'))(n=4000)

dev.off()
pdf("hm_advanced_APC.pdf", 10,6,8)
heatmap.plus(x, col=my_palette,
             scale="row", 
             key=T, 
             keysize=2, 
             symkey=T,
             density.info="none", 
             trace="none",
             cexCol=2,  
             labRow=F, #rownames(x), #rownames(x), 
             margins = c(0,15),
             main="",
             distfun  = function(x) dist(x, method="euclidean"),
             #hclustfun= function(x) hclust(x, method="complete"),
             ColSideColors=dd)
par(lend = 4)
legend("topright", 
       legend = c("Advanced", "Localized"),
       title = "Disease extent",
       col = c("Yellow", "Green"),
       lwd = 10,
       xjust = 0,
       lty = 1,
       cex=0.8)

legend("right", 
       legend = c("Yes", "No"),
       title = "APC hypermethylation",
       col = c("Blue", "Red"),
       lwd = 10,
       lty = 1,
       xjust = 0,
       cex=0.8)
dev.off()




# ARX PDX1 - probably not included in revision ######
Island_data <- readRDS("EPIC_cg_data.rds")
CG_ARX <- Island_data$Name[Island_data$UCSC_RefGene_Accession == "NM_139058" & grepl("Promoter", Island_data$Regulatory_Feature_Group) ] 

CG_PDX1 <- Island_data$Name[Island_data$UCSC_RefGene_Accession == "NM_000209" & grepl("Promoter", Island_data$Regulatory_Feature_Group)] 

beta_ARX <- beta_samples[rownames(beta_samples) %in% CG_ARX, ]
beta_PDX1 <- beta_samples[rownames(beta_samples) %in% CG_PDX1, ]
beta_ARX_PDX1 <- rbind(beta_ARX, beta_PDX1)

beta_ARX_PDX1_long <- melt(beta_ARX_PDX1, id.vars=c("genes"))


means_ARX <- colMedians(as.matrix(beta_ARX))
means_PDX1 <- colMedians(as.matrix(beta_PDX1))
ProgReg <- targets$AnyMet
type <- targets$Site
prog_PDX1_ARX <- as.data.frame(cbind(ProgReg, means_ARX, means_PDX1, means_PDX1/means_ARX,  type, targets$Sample_Group))

prog_PDX1_ARX <- prog_PDX1_ARX[prog_PDX1_ARX$V6 == "MEN1", ]

fill <- "#4271AE"
lines <- "#1F3552"

h <- ggplot(prog_PDX1_ARX, aes(x=factor(prog_PDX1_ARX$ProgReg), y=as.numeric(as.character(prog_PDX1_ARX$V4))))
h + geom_violin(colour = lines, fill = fill, size =1, scale = "width") +
  facet_wrap(~factor(prog_PDX1_ARX$type)) +
  ylab("Beta values") +
  xlab("Gene symbol") +
  #  stat_summary(aes(group = "genes"), fun.y=mean, geom="line", colour="green") 
  geom_jitter(shape=16, position=position_jitter(0.25)) + 
  stat_compare_means(label = "p.signif", method = "wilcox.test") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14), 
    axis.text = element_text(size = 12),
    legend.position = "NONE") 

beta_samples <- readRDS("beta_96.rds") %>% as.data.frame()
