library(ggplot2)

df <- read.csv("dataframe_barplots_placement_stats.tsv", sep = "\t")

order_vec <- c("LP_consistency","LP_2_nodes","GP_40_consistency",
               "GP_40_2_nodes","GP_100_consistency","GP_100_2_nodes")

df$analysis <- factor(df$analysis, levels = order_vec)

order_vec_country <- c("Spain", "Thailand", "Moldova")
df$country <- factor(df$country, levels = order_vec_country)

ggplot(df, aes(fill = analysis, y = value, x = country)) + 
  geom_bar(position=position_dodge(width = 0.6), stat="identity", width = 0.5) +
  theme_bw() + 
  scale_fill_manual(values = c("#cad2c5","#a7bea9","#84a98c",
                               "#52796f","#354f52","#2f3e46"))

### Barplot agrupado por categorías de análisis y separado entre consistency
### y 2 nodes (el primer plot serían 3 grupos con 3 barras cada uno (1 por país))
### y el segundo sería un grupo de 3 barras

df <- read.csv("dataframe_barplots_placement_stats.tsv", sep = "\t")

df_consistency <- df[df$analysis %in% c("LP_40_consistency", "GP_40_consistency", "GP_all_consistency"),]
order_vec_analysis <- c("LP_40_consistency","GP_40_consistency","GP_all_consistency")
order_vec_country <- c("Spain","Thailand","Moldova")
df_consistency$analysis <- factor(df_consistency$analysis, levels = order_vec_analysis)
df_consistency$country <- factor(df_consistency$country, levels = order_vec_country)

df_2_nodes <- df[df$analysis %in% c("LP_40_2_nodes", "GP_40_2_nodes", "GP_all_2_nodes"),]
order_vec_analysis <- c("LP_40_2_nodes","GP_40_2_nodes","GP_all_2_nodes")
df_2_nodes$analysis <- factor(df_2_nodes$analysis, levels = order_vec_analysis)
df_2_nodes$country <- factor(df_2_nodes$country, levels = order_vec_country)



ggplot(df_consistency, aes(fill = country, y = value, x = analysis)) + 
  geom_bar(position=position_dodge(width = 0.6), stat="identity", width = 0.5) +
  theme_bw() + scale_fill_manual(values = c("#264653","#2a9d8f","#e9c46a")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(df_2_nodes, aes(fill = country, y = value, x = analysis)) + 
  geom_bar(position=position_dodge(width = 0.6), stat="identity", width = 0.5) +
  theme_bw() + scale_fill_manual(values = c("#264653","#2a9d8f","#e9c46a")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
