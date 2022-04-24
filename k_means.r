# from https://www.youtube.com/watch?v=NKQpVU1LTm8

install.packages("factoextra")
library(factoextra)

horvath_levine_csv_dir <- file.path(getwd(), "csv")
horvath_beta_csv_path <- file.path(horvath_levine_csv_dir, "horvath_beta.csv") 
horvath_beta <- read.csv(file=horvath_beta_csv_path, row.names = 1) 
horvath_beta_t <- data.frame(t(horvath_beta))
horvath_beta_t_scale <- scale(horvath_beta_t)
horvath_beta_t_scale_dist <- dist(horvath_beta_t_scale)

fviz_nbclust(horvath_beta_t_scale, kmeans, method = "wss") +
  labs(subtitle = "Elbow method")

km.out <- kmeans(horvath_beta_t_scale, centers = 6, nstart = 100)
#print(km.out)

km.clusters <- km.out$cluster
rownames(horvath_beta_t_scale) <- paste(metadata_file$Sample_Group, 1:dim(horvath_beta_t)[1], sep = "_")
#fviz_cluster(list(data=horvath_beta_t_scale, cluster = km.clusters))
table(km.clusters, metadata_file$Sample_Group)

#######

# km.clusters MEN1 Sporadic VHL
          # 1    2        1   7
          # 2    3       16   0
          # 3   26        3   0
          # 4    3        3   3
          # 5    2        5   0
          # 6    6       16   0

# km.clusters MEN1 Sporadic VHL
          # 1    2        5   0
          # 2    3        3   3
          # 3   15       17   0
          # 4   20       18   0
          # 5    2        1   7

# km.clusters MEN1 Sporadic VHL
          # 1   31       34   0
          # 2    3        3   3
          # 3    6        1   7
          # 4    2        6   0

# km.clusters MEN1 Sporadic VHL
          # 1    2        6   0
          # 2   37       34   7
          # 3    3        4   3