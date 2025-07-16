library(stringr)

rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)

file1 <- read.csv(args[1], sep = "\t", header = T) 
file2 <- read.csv(args[2], sep = "\t", header = T)

matching_file1 <- c()
matching_file2 <- c()

# We generate a new file where we rename the matching nodes. First, we save the matching clusters and then the most similar ones.

count <- 0
vec_cluster_file1 <- c()
vec_cluster_file2 <- c()

for(i in 1:nrow(file1)){
  row1 <- file1[i,]
  cluster1 <- row1$Cluster
  samples1 <- unlist(str_split(row1$Samples, ","))

  for(j in 1:nrow(file2)){
    row2 <- file2[j,]
    cluster2 <- row2$Cluster
    samples2 <- unlist(str_split(row2$Samples, ","))
    
    if(all(sort(samples1) == sort(samples2))){
      print(paste("Cluster ",cluster1," matches cluster ",cluster2, sep =""))
      matching_file1 <- c(matching_file1, cluster1)
      matching_file2 <- c(matching_file2, cluster2)
    } 
  }
}

print(paste("Standard pipeline Clusters: ",nrow(file1), sep =""))
print(paste("UShER clusters: ",nrow(file2), sep =""))
print(paste("Equal clusters between both approaches: ", length(matching_file1), sep =""))

# Non-matching clusters

nonmatching_file1 <- setdiff(file1$Cluster, nonmatching_file1)
nonmatching_file2 <- setdiff(file2$Cluster, nonmatching_file2)

print("Non-matching clusters:")
print(paste("Standard pipeline clusters: ", paste(nonmatching_file1, collapse = ", "), sep =""))
print(paste("UShER clusters: ", paste(nonmatching_file2,collapse = ", "), sep =""))
