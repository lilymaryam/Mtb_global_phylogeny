library(stringr)
rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)

# Reference clusters
file1 <- read.csv(args[1], sep = "\t", header = T)

# Clusters from placement phylogeny

file2 <- read.csv(args[2], sep = "\t", header = T)

matching <- c()

# We rename the matching clusters the same way to eventually get the Cohen's Kappa 

df_rename <- as.data.frame(matrix(ncol = 2, nrow = 1))
colnames(df_rename) <- c("Renamed_cluster", "Samples")

count <- 0
vec_cluster_file1 <- c()
vec_cluster_file2 <- c()

for(i in 1:nrow(file1)){
  row1 <- file1[i,]
  cluster1 <- row1$Cluster
  samples1 <- unlist(str_split(row1$Samples, ","))
  vec_cluster_file1 <- c(vec_cluster_file1, samples1)
  for(j in 1:nrow(file2)){
    row2 <- file2[j,]
    cluster2 <- row2$Cluster
    samples2 <- unlist(str_split(row2$Samples, ","))
    vec_cluster_file2 <- c(vec_cluster_file2, samples2)
 
 if(all(sort(samples1) == sort(samples2))){
      print(paste("Reference cluster ", cluster1, "matches ", cluster2, sep = ""))
      count <- count+1
      matching <- c(matching, cluster1)
      new_name <- paste("Renamed_cluster_", count, sep ="")
      df_rename[count,1] <- new_name
      df_rename[count,2] <- paste(unlist(str_split(samples1, ",")), collapse = ",")
    }
  }
}

non_matching <- file1$Cluster[!(file1$Cluster %in% matching)]


vec_cluster_file1 <- sort(unique(vec_cluster_file1))
vec_cluster_file2 <- sort(unique(vec_cluster_file2))
vec_combined <- sort(unique(c(vec_cluster_file1, vec_cluster_file2))) # Samples in cluster


# Non-shared clusters

clusters_file1 <- c()
clusters_file2 <- c()
sample_vec <- c()
count <- 1
count_1 <- 0
count_2 <- 0
count_3 <- 0 

for(sample in vec_combined){
  sample_vec <- c(sample_vec, sample)

  if(sum(grepl(sample, vec_cluster_file1)) == 0){
    clusters_file1 <- c(clusters_file1, "Not_in_cluster_in_regular_phylo")
    count <- count+1
    count_1 <- count_1+1
    next
  }
  cluster <- df_rename[grepl(paste("\\b",sample,"\\b", sep = ""), df_rename$Samples),]$Renamed_cluster
  if(length(cluster) == 1){

    count <- count + 1
    clusters_file1 <- c(clusters_file1,cluster)
    count_2 <- count_2 + 1
  } else{
    cluster <- file1[grepl(paste("\\b",sample,"\\b", sep = ""), file1$Samples),]$Cluster

    clusters_file1 <- c(clusters_file1, cluster)
    count_3 <- count_3 + 1
    count <- count + 1
  }
}

for(sample in vec_combined){
  print(sample)
  if(sum(grepl(sample, vec_cluster_file2)) == 0){
    clusters_file2 <- c(clusters_file2, "Not_in_cluster_in_placed_phylo")
    count <- count+1
    next
  }
  cluster <- df_rename[grepl(paste("\\b",sample,"\\b", sep = ""), df_rename$Samples),]$Renamed_cluster
  if(length(cluster) == 1){
    clusters_file2 <- c(clusters_file2,cluster)
  } else{
    cluster <- file2[grepl(paste("\\b",sample,"\\b", sep = ""), file2$Samples),]$Cluster
    if(length(cluster) == 0){
      clusters_file2 <- c(clusters_file2, paste(sample, "_", count, sep = ""))
    }
    clusters_file2 <- c(clusters_file2, cluster)
    count <- count + 1
  }
}

library(irr)

cluster_df <- data.frame(Samples = sample_vec, cluster_1 = clusters_file1, cluster_2 = clusters_file2)

# We save a table containing all sample in transmission, and the cluster they're assigned too. If the cluster matches between both approaches,
# the cluster name will be the same (to ease cohen's kappa calculation). For non-matching clusters, it's needed to revise them and check if they
# are completely different or if they highly overlap, in which case the cluster ID will modified to be the same for that sample, since
# the sample is correctly assigned to the cluster in both cases despite being the cluster a bit different

write.csv(x = cluster_df, paste("temp_table.csv", file1, sep = ""), row.names = F)

# Once we've revised the table, we import it to get the Cohen's Kappa value

cluster_df <- read.csv("temp_table.csv", sep = ",")

kappa_result <- kappa2(cluster_df[,c(2,3)], "unweighted")
