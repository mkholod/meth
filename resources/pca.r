# mtcars[,c(1:7,10,11)]
# my_data.pca <- prcomp(data_with_epi_toc[,c(2:5, 10:13, 31, 48:55)], center = TRUE,scale. = TRUE)
my_data.pca <- prcomp(data_with_epi_toc[,c(2:5, 10:13, 31, 48:51, 54:55)], center = TRUE,scale. = TRUE)
summary(my_data.pca)

# Importance of components:
#   PC1    PC2    PC3    PC4    PC5     PC6     PC7     PC8     PC9    PC10     PC11
# Standard deviation     2.3159 1.9055 1.6944 1.4754 0.8423 0.36579 0.26138 0.18771 0.09602 0.03957 0.001565
# Proportion of Variance 0.3576 0.2421 0.1914 0.1451 0.0473 0.00892 0.00455 0.00235 0.00061 0.00010 0.000000
# Cumulative Proportion  0.3576 0.5996 0.7910 0.9362 0.9835 0.99238 0.99693 0.99928 0.99990 1.00000 1.000000
# PC12      PC13      PC14      PC15
# Standard deviation     1.192e-15 3.317e-16 2.371e-16 8.581e-17
# Proportion of Variance 0.000e+00 0.000e+00 0.000e+00 0.000e+00
# Cumulative Proportion  1.000e+00 1.000e+00 1.000e+00 1.000e+00

library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)

ggbiplot(my_data.pca)
# ggbiplot(my_data.pca, labels=rownames(data_with_epi_toc)) # very noisy
ggbiplot(my_data.pca, ellipse=TRUE, groups=data_with_epi_toc$Sample_Group)
ggbiplot(my_data.pca, ellipse=TRUE, groups=data_with_epi_toc$Site)
ggbiplot(my_data.pca, ellipse=TRUE, groups=data_with_epi_toc$AnyMet)
ggbiplot(my_data.pca, ellipse=TRUE, groups=data_with_epi_toc$LN)
ggbiplot(my_data.pca, ellipse=TRUE, groups=data_with_epi_toc$Dis)
ggbiplot(my_data.pca, ellipse=TRUE, groups=data_with_epi_toc$Grade)
ggbiplot(my_data.pca, ellipse=TRUE, groups=data_with_epi_toc$Functional)
ggbiplot(my_data.pca, ellipse=TRUE, groups=data_with_epi_toc$Sample_Well)

horvath.pca <- prcomp(data_with_epi_toc[,c(2:5, 10:13, 31, 48:51, 54:55)], center = TRUE,scale. = TRUE)

