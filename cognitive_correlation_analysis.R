library(tidyverse)
library(readxl)

# Correlation analysis between intelligence and compounds
# Batch 2
cpd_data <- read.csv('../DATA/peakmat_features/batch2_reorganize_20241107/batch2_compound_data.csv')
sample_id <- read_excel('../DATA/sample_id.xlsx')
wisc <- read_excel('../DATA/m48WISC01_hxq.xlsx')
cpd_info <- read_excel('../DATA/statistics/identified_compounds_information_batch2.xlsx', sheet=1)

merge_df <- cpd_data %>%
  inner_join(sample_id %>% select(ID, sample_name), by='sample_name') %>%
  inner_join(wisc %>% select(ID, VCI, VSI, FRI, WMI, PSIRevised, FSIQ), by='ID')
# change blank to mean in VSI, WMI, PSIRevised
merge_df$VSI[is.na(merge_df$VSI)] <- mean(merge_df$VSI, na.rm=TRUE)
merge_df$WMI[is.na(merge_df$WMI)] <- mean(merge_df$WMI, na.rm=TRUE)
merge_df$PSIRevised[is.na(merge_df$PSIRevised)] <- mean(merge_df$PSIRevised, na.rm=TRUE)

# ggqqplot(merge_df$PFT0011, ylab='metabolite')
cor_cpd_iq <- cor(merge_df[, 2:500], merge_df[, 502:507], method='spearman')

# Get the correlation coefficients and p-values
cor_coefficients <- matrix(NA, ncol = 6, nrow = 499)
p_values <- matrix(NA, ncol = 6, nrow = 499)

for (i in 2:500) {
  for (j in 502:507) {
    test_result <- cor.test(merge_df[[i]], merge_df[[j]], method = 'spearman')
    cor_coefficients[i-1, j-501] <- test_result$estimate
    p_values[i-1, j-501] <- test_result$p.value
  }
}
cor_coefficients_df <- data.frame(cor_coefficients)
p_values_df <- data.frame(p_values)

colnames(cor_coefficients_df) <- paste0(c('VCI', 'VSI', 'FRI', 'WMI', 'PSI', 'FSIQ'), "_cor")
rownames(cor_coefficients_df) <- colnames(merge_df)[2:500]

colnames(p_values_df) <- paste0(c('VCI', 'VSI', 'FRI', 'WMI', 'PSI', 'FSIQ'), "_p")
rownames(p_values_df) <- colnames(merge_df)[2:500]

combined_df <- cbind(cor_coefficients_df, p_values_df)

combined_df$VCI_ap <- p.adjust(combined_df$VCI_p, method = 'BH')
combined_df$VSI_ap <- p.adjust(combined_df$VSI_p, method = 'BH')
combined_df$FRI_ap <- p.adjust(combined_df$FRI_p, method = 'BH')
combined_df$WMI_ap <- p.adjust(combined_df$WMI_p, method = 'BH')
combined_df$PSI_ap <- p.adjust(combined_df$PSI_p, method = 'BH')
combined_df$FSIQ_ap <- p.adjust(combined_df$FSIQ_p, method = 'BH')

combined_df <- cbind(unique_features = rownames(combined_df), combined_df)
combined_df <- combined_df %>%
  inner_join(cpd_info %>% select(Name, unique_features), by='unique_features')

write.csv(combined_df, file='./20241107_reanalysis/correlation_metab_iq_batch2.csv', row.names = FALSE)

# Correlation between four intelligence scores
library(corrplot)
corr_M <- cor(merge_df[, 507:512], method='spearman')
colnames(corr_M)[colnames(corr_M) == "PSIRevised"] <- "PSI"
rownames(corr_M)[rownames(corr_M) == "PSIRevised"] <- "PSI"
tiff(filename = "../article/pictures/20241107_reorganize/correlation/iq_correlation_batch2.tiff",
     width = 3.5, height = 3.5, units = "in", res = 300)
corrplot(corr_M, method = 'number', type='upper')
dev.off()

# Significant compounds in Batch 2
combined_df <- read.csv('./20241107_reanalysis/correlation_metab_iq_batch2.csv')
sig_cpd <- combined_df %>%
  filter(VCI_ap < 0.05 | VSI_ap < 0.05 | PSI_ap < 0.05 |FRI_ap < 0.05|WMI_ap < 0.05|FSIQ_ap < 0.05) 

# Heatmap
library(ComplexHeatmap)

sig_cpd <- combined_df %>%
  filter(VCI_ap < 0.05 | VSI_ap < 0.05 | PSI_ap < 0.05 |FRI_ap < 0.05|WMI_ap < 0.05|FSIQ_ap < 0.05) %>%
  inner_join(cpd_info %>% select(Name, Class), by='Name') %>%
  mutate(Name = ifelse(Name == 'Perfluorononanoic acid', 'PFNA', Name)) %>%
  mutate(Name = ifelse(Name == 'Perfluorodecanoic acid', 'PFDA', Name)) %>%
  mutate(Name = ifelse(Name == 'Perfluoroundecanoic acid', 'PFUnDA', Name)) %>%   # Perfluorooctanesulfonic acid, PFOS
  mutate(Name = ifelse(Name == 'Linoleic acid', 'FA 18:2', Name)) %>%
  mutate(Class = ifelse(Class == 'Glycerophosphocholines', 'PC', Class)) %>%
  mutate(FSIQ_cor_sign = ifelse(FSIQ_cor >= 0, "Positive", "Negative")) %>%
  mutate(Class = ifelse(Name == 'Ethyl 2Z,4E-decadienoic acid', 'Fatty acyls', Class),
         Class = ifelse(str_ends(Name, 'ethanolamide'), 'NAE', Class)) %>%
  # mutate(Class = factor(Class)) %>%
  arrange(desc(FSIQ_cor_sign), Class, Name)
  
sig_cpd_cor <- t(sig_cpd[,2:7])
row.names(sig_cpd_cor) <- c('VCI', 'VSI', 'FRI', 'WMI', 'PSI', 'FSIQ')
colnames(sig_cpd_cor) <- sig_cpd$Name

sig_cpd_p <- t(sig_cpd[, 8:13])

class_cols = list(Class = c(`Amino acids` = '#a6cee3',
                            `Fatty acids` = '#33a02c',
                            `Organic acids` = '#fdbf6f',
                            Others = '#ff7f00',
                            PC = '#003c30',
                            Xenobiotics  = '#b15928',
                            `Fatty acyls` = '#fb9a99',
                            NAE = '#cab2d6'
                            ))
column_ha = HeatmapAnnotation(Class=sig_cpd$Class,
                              annotation_legend_param = list(
                                Class = list(
                                  title='Class',
                                  at = unique(sig_cpd$Class),
                                  labels=unique(sig_cpd$Class),
                                  nrow=1,
                                  direction='horizontal'
                                )
                              ),
                              show_annotation_name = FALSE,
                              col = class_cols,
                              text = anno_text(colnames(sig_cpd_cor), rot=-45, just='left'),
                              show_legend = FALSE
)
ht = Heatmap(sig_cpd_cor, cluster_rows = FALSE,show_column_names = FALSE,cluster_columns = FALSE,
             heatmap_legend_param =list(title="Spearman Correlation Coefficient", direction='horizontal'),
             cell_fun = function(j, i, x, y, width, height, fill) {
               if (sig_cpd_p[i, j] >= 0.01 & sig_cpd_p[i, j] < 0.05) {
                 grid.text('*', x, y, gp = gpar(fontsize=10))
               }
               if (sig_cpd_p[i, j] >= 0.001 & sig_cpd_p[i, j] < 0.01) {
                 grid.text('**', x, y, gp = gpar(fontsize=10))
               }
               if (sig_cpd_p[i, j] < 0.001) {
                 grid.text('***', x, y, gp = gpar(fontsize=10))
               }
             },
             bottom_annotation = column_ha,
             row_names_side = 'right',
             show_heatmap_legend = FALSE)

default_color_fun <- ht@matrix_color_mapping@col_fun
ht_legend = Legend(col_fun = default_color_fun, title = "Spearman Correlation Coefficient", 
                   direction = "horizontal")

class_legend = Legend(labels = unique(sig_cpd$Class), title = "Class", 
                      legend_gp = gpar(fill = unlist(class_cols$Class)),
                      direction = "horizontal", nrow = 1)

# Combine the legends using packLegend
combined_legend = packLegend(ht_legend, class_legend, direction = "horizontal", gap = unit(5, "mm"))

png("../article/pictures/20241107_reorganize/correlation/corr_heatmap_20241118.png",
    width=15, height=5, units="in", res=300)
draw(ht,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     annotation_legend_list = combined_legend,  # Correctly set the combined legend here
     merge_legends = TRUE
)
dev.off()
