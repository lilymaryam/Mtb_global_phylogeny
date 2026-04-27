### RAND-INDEX

# Install mclust if you haven't: install.packages("mclust")
library(tidyverse)
library(mclust)
 
# 1. Load data
gs_df <- read.csv("/Users/fran/Downloads/zvalencia_multifasta_h37rv.clusters_10.tsv", sep = "\t")
m2_df <- read.csv("/Users/fran/Downloads/valencia_dataset.masked10.clusters_10.tsv", sep = "\t")

# 2. Reshape to long format (Sample | ClusterID)
gs_long <- gs_df %>% separate_rows(Samples, sep = ",")
m2_long <- m2_df %>% separate_rows(Samples, sep = ",")

# 3. Create a master list of all samples
all_samples <- unique(c(gs_long$Samples, m2_long$Samples))
comparison_df <- data.frame(Sample = all_samples)

# 4. Join the cluster IDs
comparison_df <- comparison_df %>%
  left_join(gs_long, by = c("Sample" = "Samples")) %>%
  rename(GS_Cluster = Cluster) %>%
  left_join(m2_long, by = c("Sample" = "Samples")) %>%
  rename(M2_Cluster = Cluster)

# 5. Handle singletons (NAs)
# Samples not in a cluster are unique entities. 
# We assign them a unique dummy ID so they aren't grouped together.
comparison_df <- comparison_df %>%
  mutate(
    GS_Cluster = ifelse(is.na(GS_Cluster), paste0("S_GS_", row_number()), GS_Cluster),
    M2_Cluster = ifelse(is.na(M2_Cluster), paste0("S_M2_", row_number()), M2_Cluster)
  )

# 6. Calculate ARI
ari_score <- adjustedRandIndex(comparison_df$GS_Cluster, comparison_df$M2_Cluster)

cat(paste("Adjusted Rand Index:", round(ari_score, 4)))



### PRECISION/RECALL

library(tidyverse)

# 1. Load Data
gs_raw <- read.csv("/Users/fran/Downloads/zvalencia_multifasta_h37rv.clusters_10.tsv", sep = "\t")
m2_raw <- read.csv("/Users/fran/Downloads/valencia_dataset.masked10.clusters_10.tsv", sep = "\t")

# Helper: calculate n*(n-1)/2
count_pairs <- function(n) { (n * (n - 1)) / 2 }

# 2. Total pairs that SHOULD exist (Gold Standard)
gs_total_pairs <- gs_raw %>%
  mutate(n = str_count(Samples, ",") + 1) %>%
  summarise(total = sum(count_pairs(n))) %>%
  pull(total)

# 3. Total pairs Method 2 CLAIMED exist
m2_total_pairs <- m2_raw %>%
  mutate(n = str_count(Samples, ",") + 1) %>%
  summarise(total = sum(count_pairs(n))) %>%
  pull(total)

# 4. Calculate True Positives (TP)
# We find the intersection: how many samples from GS_Cluster_A are also in M2_Cluster_B?
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

# 5. Metrics
precision <- tp / m2_total_pairs
recall <- tp / gs_total_pairs
f1 <- 2 * (precision * recall) / (precision + recall)

cat(sprintf("Precision: %.4f\nRecall:    %.4f\nF1-Score:  %.4f\n", precision, recall, f1))



