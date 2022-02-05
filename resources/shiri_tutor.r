resources_dir <- as.name("~/Documents/meth/resources")
setwd(paste(resources_dir))

BiocManager::install(c("DESeq2")) # install.packages("DESeq2") # doesn't work
install.packages("magick") # takes more than 60 seconds to download
install.packages("gplots")
install.packages("ggforce")
BiocManager::install(c("Vennerable")) # install.packages("Vennerable") # is too  old

source("deseq2_packages.R") # can't use 
source("basics.R") # Error in library("Vennerable") : there is no package called ‘Vennerable’
source("pathway_and_annotation.R")
source("heatmaps_packages.R")
source("methylation.R")
source("statistics_packages.r")
library("methylclock")

levine <- read.csv("D:/Users/user/Downloads/levine_beta.csv")
horvath <- read.csv("D:/Users/user/Downloads/horvath_beta.csv")
mat_96 <- readRDS("D:/Users/user/Desktop/beta_96.rds")
metadata <- read.csv("D:/Users/user/Desktop/metadata.csv")
condition <- metadata$Sample_Group

#cg_horvath <- coefHorvath$CpGmarker
#saveRDS(cg_horvath, "cg_horvath.rds")
#cg_levine <- coefLevine$CpGmarker
#saveRDS(cg_levine, "cg_levine.rds")

#beta_horvath <- mat_96[which(rownames(mat_96) %in% cg_horvath), ]
#beta_levine <- mat_96[which(rownames(mat_96) %in% cg_levine), ]
#age <- metadata$age
#clocks <- DNAmAge(mat_96, age=age)

clocks$Horvath
clocks$Hannum
clocks$Levine

matrix <- mat_96
matrix <- beta_horvath
matrix <- beta_levine

z <- rowFtests(matrix, as.factor(condition), var.equal = FALSE)
z <- z[!is.na(z$p.value), ] %>% droplevels()
matrix <- matrix[!is.na(z$p.value), ] 
select <- order(z$p.value, decreasing=F)[1:length(z$p.value[z$p.value<0.05 & !is.na(z$p.value)])]

matrix1 <- matrix[select , ] %>% t()
pc <- prcomp(matrix1)
scores <- data.frame(pc$x, condition)
#library(ggallin)
#dev.off()
pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition)))) + 
  geom_point(size = 4) + stat_ellipse(level = 0.95) + ggtitle("Principal Component Analysis") + 
  scale_colour_brewer(name = " ", palette = "Set1") + 
  #  scale_y_continuous(trans=pseudolog10_trans) + scale_x_continuous(trans=pseudolog10_trans) +
  theme(
    plot.title = element_text(face = 'bold'),
    legend.position = c(.87,.9),
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
    panel.background = element_rect(color = 'black',fill = NA)) 
pcaplot