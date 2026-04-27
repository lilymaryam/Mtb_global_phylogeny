# This script calculates the Adjusted Rand Index (ARI) and Precision/Recall metrics
# for comparing two clustering results: a gold standard (GS) and a method (M2)

### ADJUSTED RAND-INDEX

library(tidyverse)
library(mclust)
 
# Let's load the data
gs_df <- read.csv("clusters_from_multifasta_GS.tsv", sep = "\t") # GS clusters
m2_df <- read.csv("clusters_from_ushersampled_phylogeny_M2.tsv", sep = "\t") # M2 clusters

# We convert the cluster assignments into a long format where each sample has its own row
gs_long <- gs_df %>% separate_rows(Samples, sep = ",")
m2_long <- m2_df %>% separate_rows(Samples, sep = ",")

# Now we create a combined dataframe to compare the cluster assignments for each sample
all_samples <- unique(c(gs_long$Samples, m2_long$Samples))
comparison_df <- data.frame(Sample = all_samples)

# We join the cluster assignments from both GS and M2 to this combined dataframe
comparison_df <- comparison_df %>%
  left_join(gs_long, by = c("Sample" = "Samples")) %>%
  rename(GS_Cluster = Cluster) %>%
  left_join(m2_long, by = c("Sample" = "Samples")) %>%
  rename(M2_Cluster = Cluster)

# Some samples might be missing from either GS or M2, we can assign them to a unique cluster 
# (e.g., "S_GS_X" or "S_M2_X") to ensure they are treated as singletons in the adjusted Rand-index calculation
comparison_df <- comparison_df %>%
  mutate(
    GS_Cluster = ifelse(is.na(GS_Cluster), paste0("S_GS_", row_number()), GS_Cluster),
    M2_Cluster = ifelse(is.na(M2_Cluster), paste0("S_M2_", row_number()), M2_Cluster)
  )

# Now we can calculate the Adjusted Rand Index
ari_score <- adjustedRandIndex(comparison_df$GS_Cluster, comparison_df$M2_Cluster)

cat(paste("Adjusted Rand Index:", round(ari_score, 4)))



### PRECISION/RECALL

library(tidyverse)

# Let's load the data again (or you can reuse the already loaded dataframes)
gs_raw <- read.csv("clusters_from_multifasta_GS.tsv", sep = "\t")
m2_raw <- read.csv("clusters_from_ushersampled_phylogeny_M2.tsv", sep = "\t")

# Function to calculate the number of pairs in a cluster of size n
count_pairs <- function(n) { (n * (n - 1)) / 2 }

# Total pairs in GS
gs_total_pairs <- gs_raw %>%
  mutate(n = str_count(Samples, ",") + 1) %>%
  summarise(total = sum(count_pairs(n))) %>%
  pull(total)

# Total pairs in M2
m2_total_pairs <- m2_raw %>%
  mutate(n = str_count(Samples, ",") + 1) %>%
  summarise(total = sum(count_pairs(n))) %>%
  pull(total)

# To calculate true positives, we need to find how many pairs of samples are clustered together in both GS and M2.

gs_long <- gs_raw %>% separate_rows(Samples, sep = ",")
m2_long <- m2_raw %>% separate_rows(Samples, sep = ",")

tp <- inner_join(gs_long, m2_long, by = "Samples") %>%
  group_by(Cluster.x, Cluster.y) %>%
  tally() %>%
  filter(n > 1) %>%
  mutate(pairs = count_pairs(n)) %>%
  summarise(sum_tp = sum(pairs)) %>%
  summarise(grand_tp = sum(sum_tp)) %>%
  pull(grand_tp)

# Now we can calculate precision and recall
precision <- tp / m2_total_pairs
recall <- tp / gs_total_pairs

cat(sprintf("Precision: %.4f\nRecall:    %.4f", precision, recall))



