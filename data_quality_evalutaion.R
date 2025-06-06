library(tidyverse)
library(FactoMineR)  # For PCA
library(factoextra)  # For PCA visualization
library(readxl)

{
  # Raw area and IS and QC correction
  # row: features, column: samples
  pos_raw <- read.csv('../DATA/peakmat_features/batch2_reorganize_20241107/batch2_pos_mix_correction_include_qc.csv',
                      check.names = FALSE)  # data columns: 2-1369
  neg_raw <- read.csv('../DATA/peakmat_features/batch2_reorganize_20241107/batch2_neg_mix_correction_include_qc.csv',
                      check.names = FALSE)  # data columns: 2-1374
  pos_sequence <- read.csv('../DATA/对外测试的积分原始数据/测样序列/20220919 xinhua sequence pos 2051(20220407).csv',
                           skip = 1) %>%
    mutate(order = 1:nrow(.))
  neg_sequence <- read.csv('../DATA/对外测试的积分原始数据/测样序列/20220919 xinhua sequence neg 2051(20220506).csv',
                           skip=1) %>%
    mutate(order = 1:nrow(.))
  
  pos_labels <- colnames(pos_raw)[2:1369] %>%
    tibble(column_name = .) %>%
    mutate(
      group = case_when(
        str_starts(column_name, "QC") ~ "QC",
        str_starts(column_name, "S") ~ "Sample"
      ))
  
  neg_labels <- colnames(neg_raw)[2:1374] %>%
    tibble(column_name = .) %>%
    mutate(
      group = case_when(
        str_starts(column_name, "QC") ~ "QC",
        str_starts(column_name, "S") ~ "Sample"
      ))
  
  # Compute QC rsd and fiter the feature greater than 0.3
  pos_qc_data <- pos_raw[, pos_labels$column_name[pos_labels$group == 'QC']]
  pos_rsd <- apply(pos_qc_data, 1, function(x) sd(x) / mean(x))
  pos_rsd_cpd <- data.frame(unique_features = pos_raw$id,
                            qc_rsd = pos_rsd) %>%
    mutate(unique_features = paste0('P', unique_features)) %>%
    inner_join(read_excel('../DATA/statistics/identified_compounds_information_batch2.xlsx', sheet=1) %>%
                 select(unique_features, Name),
               by='unique_features') %>%
    filter(qc_rsd <= 0.3)
  
  neg_qc_data <- neg_raw[, neg_labels$column_name[neg_labels$group == 'QC']]
  neg_rsd <- apply(neg_qc_data, 1, function(x) sd(x) / mean(x))
  
  neg_rsd_cpd <- data.frame(unique_features = neg_raw$id,
                            qc_rsd = neg_rsd) %>%
    mutate(unique_features = paste0('N', unique_features)) %>%
    inner_join(read_excel('../DATA/statistics/identified_compounds_information_batch2.xlsx', sheet=1) %>%
                 select(unique_features, Name),
               by='unique_features') %>%
    filter(qc_rsd <= 0.3)
  
  
  # 1. normalization: sum the rows, each value should divide by the sum
  pos_normalized <- pos_raw %>%
    mutate(id = paste0('P', id)) %>%
    filter(id %in% pos_rsd_cpd$unique_features)
    
  pos_normalized[2:1369] <- sweep(pos_normalized[2:1369], 1, rowSums(pos_normalized[2:1369]),
                                  FUN = "/")
  
  neg_normalized <- neg_raw %>%
    mutate(id = paste0('N', id)) %>%
    filter(id %in% neg_rsd_cpd$unique_features)
  neg_normalized[2:1374] <- sweep(neg_normalized[2:1374], 1, rowSums(neg_normalized[2:1374]),
                                  FUN = "/")
  
  # 2. Perform PCA and make the score plot, label the samples with group of QC and samples
  
  
  pos_data_pca <- pos_normalized[, 2:1369] %>%
    t() %>%
    as.data.frame()
  pos_pca_result <- PCA(pos_data_pca, graph=FALSE)
  
  variance_explained <- pos_pca_result$eig[, 2]  # Percentage of variance explained for each PC
  pc1_var <- round(variance_explained[1], 2)    # Variance explained by PC1
  pc2_var <- round(variance_explained[2], 2)    # Variance explained by PC2
  
  pos_pca_scores <- as.data.frame(pos_pca_result$ind$coord)
  pos_pca_scores$Group <- pos_labels$group
  
  ggplot(pos_pca_scores, aes(x = Dim.1, y = Dim.2, color = Group)) +
    scale_color_manual(values = c('QC'='blue',
                                  'Sample' = 'gray28')) +
    geom_point(size=0.5) +
    stat_ellipse(aes(group = Group), linetype = "dashed", size = 0.7, level = 0.95) +
    theme_minimal() +
    labs(
      x = paste0("PC1 (", pc1_var, "%)"),
      y = paste0("PC2 (", pc2_var, "%)")
    )
  
  ggsave("../article/pictures/20241107_reorganize/data_quality/is_qc_exclude_some_pca_pos.png",
         width = 5,
         height =5,
         units = c("in"),
         dpi = 300,
         bg='white')
  
  neg_data_pca <- neg_normalized[, 2:1374] %>%
    t() %>%
    as.data.frame()
  neg_pca_result <- PCA(neg_data_pca, graph=FALSE)
  
  variance_explained <- neg_pca_result$eig[, 2]  # Percentage of variance explained for each PC
  pc1_var <- round(variance_explained[1], 2)    # Variance explained by PC1
  pc2_var <- round(variance_explained[2], 2)    # Variance explained by PC2
  
  neg_pca_scores <- as.data.frame(neg_pca_result$ind$coord)
  neg_pca_scores$Group <- neg_labels$group
  
  ggplot(neg_pca_scores, aes(x = Dim.1, y = Dim.2, color = Group)) +
    scale_color_manual(values = c('QC'='blue',
                                  'Sample' = 'gray28')) +
    geom_point(size=0.5) +
    stat_ellipse(aes(group = Group), linetype = "dashed", size = 0.7, level = 0.95) +
    theme_minimal() +
    labs(
      x = paste0("PC1 (", pc1_var, "%)"),
      y = paste0("PC2 (", pc2_var, "%)")
    )
  
  ggsave("../article/pictures/20241107_reorganize/data_quality/is_qc_exclude_some_pca_neg.png",
         width = 5,
         height =5,
         units = c("in"),
         dpi = 300,
         bg='white')
  
  # 3. get the first PC of PCA, plot it with the order of the sample analysis in 'File.Name' column of sequence table
  pos_pc1_scores <- as.data.frame(pos_pca_result$ind$coord[, 1]) %>%
    rownames_to_column(var = 'File.Name') %>%
    inner_join(pos_sequence, by='File.Name') %>%
    rename(PC1 = `pos_pca_result$ind$coord[, 1]`) %>%
    mutate(Sample.Type = if_else(Sample.Type == 'Unknown', 'Sample', Sample.Type)) %>%
    filter(Sample.Type == 'QC') %>%
    arrange(order) %>%
    mutate(File.Name = factor(File.Name, levels = unique(File.Name)))
  
  mean_pc1 <- mean(pos_pca_result$ind$coord[, 1])
  sd_pc1 <- sd(pos_pca_result$ind$coord[, 1])
  
  ggplot(pos_pc1_scores, aes(x = File.Name, y = PC1, color=Sample.Type)) +
    geom_point() +
    scale_color_manual(values = c('QC'='blue',
                                  'Sample' = 'gray28')) +
    geom_hline(yintercept = mean_pc1, linetype = "solid", color = "black", size = 0.6, alpha = 0.8) +  # Mean line
    geom_hline(yintercept = mean_pc1 + sd_pc1, linetype = "dashed", color = "blue", size = 0.6, alpha = 0.8) +  # +1 SD
    geom_hline(yintercept = mean_pc1 - sd_pc1, linetype = "dashed", color = "blue", size = 0.6, alpha = 0.8) +  # -1 SD
    geom_hline(yintercept = mean_pc1 + 2 * sd_pc1, linetype = "dotted", color = "red", size = 0.6, alpha = 0.8) +  # +2 SD
    geom_hline(yintercept = mean_pc1 - 2 * sd_pc1, linetype = "dotted", color = "red", size = 0.6, alpha = 0.8) +  # -2 SD
    theme_minimal() +
    labs(
      x = "Analysis Order",
      y = "PC1"
    ) +
    theme(axis.text.x = element_blank(),
          legend.title = element_blank(),
          legend.position = 'None')
  ggsave("../article/pictures/20241107_reorganize/data_quality/is_qc_pc1_exclude_some_pos_qc.png",
         width = 5.3,
         height =5,
         units = c("in"),
         dpi = 300,
         bg='white')
  
  
  pos_pc2_scores <- as.data.frame(pos_pca_result$ind$coord[, 2]) %>%
    rownames_to_column(var = 'File.Name') %>%
    inner_join(pos_sequence, by='File.Name') %>%
    rename(PC2 = `pos_pca_result$ind$coord[, 2]`) %>%
    mutate(Sample.Type = if_else(Sample.Type == 'Unknown', 'Sample', Sample.Type)) %>%
    filter(Sample.Type == 'QC') %>%
    arrange(order) %>%
    mutate(File.Name = factor(File.Name, levels = unique(File.Name)))
  
  mean_pc2 <- mean(pos_pca_result$ind$coord[, 2])
  sd_pc2 <- sd(pos_pca_result$ind$coord[, 2])
  
  ggplot(pos_pc2_scores, aes(x = File.Name, y = PC2, color=Sample.Type)) +
    geom_point() +
    scale_color_manual(values = c('QC'='blue',
                                  'Sample' = 'gray28')) +
    geom_hline(yintercept = mean_pc2, linetype = "solid", color = "black", size = 0.6, alpha = 0.8) +  # Mean line
    geom_hline(yintercept = mean_pc2 + sd_pc2, linetype = "dashed", color = "blue", size = 0.6, alpha = 0.8) +  # +1 SD
    geom_hline(yintercept = mean_pc2 - sd_pc2, linetype = "dashed", color = "blue", size = 0.6, alpha = 0.8) +  # -1 SD
    geom_hline(yintercept = mean_pc2 + 2 * sd_pc2, linetype = "dotted", color = "red", size = 0.6, alpha = 0.8) +  # +2 SD
    geom_hline(yintercept = mean_pc2 - 2 * sd_pc2, linetype = "dotted", color = "red", size = 0.6, alpha = 0.8) +  # -2 SD
    theme_minimal() +
    labs(
      x = "Analysis Order",
      y = "PC2"
    ) +
    theme(axis.text.x = element_blank(),
          legend.title = element_blank(),
          legend.position = 'None')
  ggsave("../article/pictures/20241107_reorganize/data_quality/is_qc_pc2_exclude_some_pos_qc.png",
         width = 5.3,
         height =5,
         units = c("in"),
         dpi = 300,
         bg='white')
  
  
  neg_pc1_scores <- as.data.frame(neg_pca_result$ind$coord[, 1]) %>%
    rownames_to_column(var = 'File.Name') %>%
    inner_join(neg_sequence, by='File.Name') %>%
    rename(PC1 = `neg_pca_result$ind$coord[, 1]`) %>%
    mutate(Sample.Type = if_else(Sample.Type == 'Unknown', 'Sample', Sample.Type)) %>%
    filter(Sample.Type == 'QC') %>%
    arrange(order) %>%
    mutate(File.Name = factor(File.Name, levels = unique(File.Name)))
  
  mean_pc1 <- mean(neg_pca_result$ind$coord[, 1])
  sd_pc1 <- sd(neg_pca_result$ind$coord[, 1])
  
  ggplot(neg_pc1_scores, aes(x = File.Name, y = PC1, color=Sample.Type)) +
    geom_point() +
    scale_color_manual(values = c('QC'='blue',
                                  'Sample' = 'gray28')) +
    geom_hline(yintercept = mean_pc1, linetype = "solid", color = "black", size = 0.6, alpha = 0.8) +  # Mean line
    geom_hline(yintercept = mean_pc1 + sd_pc1, linetype = "dashed", color = "blue", size = 0.6, alpha = 0.8) +  # +1 SD
    geom_hline(yintercept = mean_pc1 - sd_pc1, linetype = "dashed", color = "blue", size = 0.6, alpha = 0.8) +  # -1 SD
    geom_hline(yintercept = mean_pc1 + 2 * sd_pc1, linetype = "dotted", color = "red", size = 0.6, alpha = 0.8) +  # +2 SD
    geom_hline(yintercept = mean_pc1 - 2 * sd_pc1, linetype = "dotted", color = "red", size = 0.6, alpha = 0.8) +  # -2 SD
    theme_minimal() +
    labs(
      x = "Analysis Order",
      y = "PC1"
    ) +
    theme(axis.text.x = element_blank(),
          legend.title = element_blank(),
          legend.position = 'None')
  
  ggsave("../article/pictures/20241107_reorganize/data_quality/is_qc_pc1_exclude_some_neg_qc.png",
         width = 5.3,
         height =5,
         units = c("in"),
         dpi = 300,
         bg='white')
  
  neg_pc2_scores <- as.data.frame(neg_pca_result$ind$coord[, 2]) %>%
    rownames_to_column(var = 'File.Name') %>%
    inner_join(neg_sequence, by='File.Name') %>%
    rename(PC2 = `neg_pca_result$ind$coord[, 2]`) %>%
    mutate(Sample.Type = if_else(Sample.Type == 'Unknown', 'Sample', Sample.Type)) %>%
    filter(Sample.Type == 'QC') %>%
    arrange(order) %>%
    mutate(File.Name = factor(File.Name, levels = unique(File.Name)))
  
  mean_pc2 <- mean(neg_pca_result$ind$coord[, 2])
  sd_pc2 <- sd(neg_pca_result$ind$coord[, 2])
  
  ggplot(neg_pc2_scores, aes(x = File.Name, y = PC2, color=Sample.Type)) +
    geom_point() +
    scale_color_manual(values = c('QC'='blue',
                                  'Sample' = 'gray28')) +
    geom_hline(yintercept = mean_pc2, linetype = "solid", color = "black", size = 0.6, alpha = 0.8) +  # Mean line
    geom_hline(yintercept = mean_pc2 + sd_pc2, linetype = "dashed", color = "blue", size = 0.6, alpha = 0.8) +  # +1 SD
    geom_hline(yintercept = mean_pc2 - sd_pc2, linetype = "dashed", color = "blue", size = 0.6, alpha = 0.8) +  # -1 SD
    geom_hline(yintercept = mean_pc2 + 2 * sd_pc2, linetype = "dotted", color = "red", size = 0.6, alpha = 0.8) +  # +2 SD
    geom_hline(yintercept = mean_pc2 - 2 * sd_pc2, linetype = "dotted", color = "red", size = 0.6, alpha = 0.8) +  # -2 SD
    theme_minimal() +
    labs(
      x = "Analysis Order",
      y = "PC2"
    ) +
    theme(axis.text.x = element_blank(),
          legend.title = element_blank(),
          legend.position = 'None')
  
  ggsave("../article/pictures/20241107_reorganize/data_quality/is_qc_pc2_exclude_some_neg_qc.png",
         width = 5.3,
         height =5,
         units = c("in"),
         dpi = 300,
         bg='white')
  
  
  # 4. Compute the RSD of the QC samples, compute the percent of the features less than 0.1, 0.2, 0.3 and so on, plot it
  pos_qc_data <- pos_normalized[, pos_labels$column_name[pos_labels$group == 'QC']]
  pos_rsd <- apply(pos_qc_data, 1, function(x) sd(x) / mean(x))
  pos_rsd_cpd <- data.frame(unique_features = pos_normalized$id,
                            qc_rsd = pos_rsd) %>%
    mutate(unique_features = paste0('P', unique_features)) %>%
    inner_join(read_excel('../DATA/statistics/identified_compounds_information_batch2.xlsx', sheet=1) %>%
                 select(unique_features, Name),
               by='unique_features')
  
  rsd_bins <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, Inf)
  rsd_labels <- c('0-10%', '10-20%', '20-30%', '30-40%', '40%-50%', '>50%')
  
  df_rsd <- data.frame(RSD = pos_rsd) %>%
    mutate(RSD_Range = cut(RSD, breaks = rsd_bins, labels = rsd_labels, right = FALSE)) %>%
    group_by(RSD_Range) %>%
    summarise(Count = n()) %>%
    mutate(
      Cumulative_Count = cumsum(Count),
      Cumulative_Percentage = (Cumulative_Count / sum(Count)) * 100
    )
  
  ggplot(df_rsd, aes(x = RSD_Range, y = Count)) +
    geom_bar(stat = "identity", fill = '#1f78b4', color = "black") +
    geom_line(aes(y = Cumulative_Percentage, group = 1), color = 'orange', size = 1) +
    geom_point(aes(y = Cumulative_Percentage), color = 'orange', size = 2) +
    scale_y_continuous(
      name = 'Count',
      sec.axis = sec_axis(~ ., name = 'Cumulative Percentage (%)')
    ) +
    labs(
      x = "RSD in QC Samples",
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
      legend.position = "none"  # Remove legend
    )
  
  ggsave("../article/pictures/20241107_reorganize/data_quality/is_qc_rsd_exclude_some_pos.png",
         width = 5,
         height =5,
         units = c("in"),
         dpi = 300,
         bg='white')
  
  
  neg_qc_data <- neg_normalized[, neg_labels$column_name[neg_labels$group == 'QC']]
  neg_rsd <- apply(neg_qc_data, 1, function(x) sd(x) / mean(x))
  
  neg_rsd_cpd <- data.frame(unique_features = neg_normalized$id,
                            qc_rsd = neg_rsd) %>%
    mutate(unique_features = paste0('N', unique_features)) 
  # %>%
  #   inner_join(read_excel('../DATA/statistics/identified_compounds_information_batch2.xlsx', sheet=1) %>%
  #                select(unique_features, Name),
  #              by='unique_features')
  
  rsd_bins <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, Inf)
  rsd_labels <- c('0-10%', '10-20%', '20-30%', '30-40%', '40%-50%', '>50%')
  
  df_rsd <- data.frame(RSD = neg_rsd) %>%
    mutate(RSD_Range = cut(RSD, breaks = rsd_bins, labels = rsd_labels, right = FALSE)) %>%
    group_by(RSD_Range) %>%
    summarise(Count = n()) %>%
    mutate(
      Cumulative_Count = cumsum(Count),
      Cumulative_Percentage = (Cumulative_Count / sum(Count)) * 100
    )
  
  ggplot(df_rsd, aes(x = RSD_Range, y = Count)) +
    geom_bar(stat = "identity", fill = '#1f78b4', color = "black") +
    geom_line(aes(y = Cumulative_Percentage, group = 1), color = 'orange', size = 1) +
    geom_point(aes(y = Cumulative_Percentage), color = 'orange', size = 2) +
    scale_y_continuous(
      name = 'Count',
      sec.axis = sec_axis(~ ., name = 'Cumulative Percentage (%)')
    ) +
    labs(
      x = "RSD in QC Samples",
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
      legend.position = "none"  # Remove legend
    )
  ggsave("../article/pictures/20241107_reorganize/data_quality/is_qc_rsd_exclude_some_neg.png",
         width = 5,
         height =5,
         units = c("in"),
         dpi = 300,
         bg='white')
  
}