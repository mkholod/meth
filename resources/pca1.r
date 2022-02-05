# have beta

library(methylclockData)
# library(methylclock)
library(stringr)
library(genefilter)

load_DNAm_Clocks_data()

#### HORVATH
horvath_cgs_sorted <- str_sort(coefHorvath[-1,1])
row_names_beta_sorted <- str_sort(row.names(beta))
horvath_beta <- beta[row_names_beta_sorted %in% horvath_cgs_sorted,]


#### LEVINE
levine_cgs <- coefLevine[-1,1]
levine_cgs_sorted <- str_sort(levine_cgs)
levine_beta_row_names <- row_names_beta_sorted %in% levine_cgs_sorted
levine_beta <- beta[row_names_beta_sorted %in% levine_cgs_sorted,]

all_data = myDNAmAge_with_acceleration_age_with_metadata
condition <- all_data$Sample_Group

########
matrix <- horvath_beta
matrix <- levine_beta
matrix <- beta

z <- rowFtests(matrix, as.factor(condition), var.equal = FALSE)
z <- z[!is.na(z$p.value), ] %>% droplevels()
matrix <- matrix[!is.na(z$p.value), ] 
select <- order(z$p.value, decreasing=F)[1:length(z$p.value[z$p.value<0.05 & !is.na(z$p.value)])]

matrix1 <- matrix[select , ] %>% t()
pc <- prcomp(matrix1)
scores <- data.frame(pc$x, condition)

ggbiplot(pc, ellipse=TRUE, groups=condition, var.axes=FALSE)

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
