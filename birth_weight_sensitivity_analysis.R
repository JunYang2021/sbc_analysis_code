{
  # Batch2 data
  
  birth_data <- read.csv('../DATA/statistics/clinical_birth_data_batch2.csv')
  cpd_data <- read.csv('../DATA/peakmat_features/batch2_reorganize_20241107/batch2_compound_data.csv')
  cpd_info <- read_excel('../DATA/statistics/identified_compounds_information_batch2.xlsx', sheet=1)
  cl_data <- read.csv('../DATA/statistics/clinical_data_batch2.csv')
  
  merge_data <- cpd_data %>%
    inner_join(cl_data, by='sample_name') %>%
    inner_join(birth_data, by='sample_name') %>%
    filter(!is.na(brwe)) %>%
    mutate(growth = factor(if_else(brwe > 4000, 'Macrosomia', 'Non macrosomia'),
                           levels=c('Non macrosomia', 'Macrosomia'))) %>%
    select(-c('VCI', 'VSI', 'FRI', 'WMI', 'PSIRevised', 'FSIQ'))
  cols_to_factor <- c('chsex_new', 'yqyzdc_educton', 'yqyzdc_income', 'yqyzdc_folv4_1', 'yzhdc_FOLIC1',
                      'yzhdc_SMKN3', 'fmbs_CBH46')
  merge_data[cols_to_factor] <- lapply(merge_data[cols_to_factor], as.factor)
  cl_data[cols_to_factor] <- lapply(cl_data[cols_to_factor], as.factor)
  
  merge_data$yqyzdc_educton <- factor(merge_data$yqyzdc_educton,
                                      levels = c(2, 3, 4, 5, 6),
                                      labels = c('Low', 'Low','Low','Medium', 'High'))
  merge_data$yqyzdc_income <- factor(merge_data$yqyzdc_income,
                                     levels = c(1,2,3,4,5,6,7,8),
                                     labels =c('Low','Low','Low','Low','Low',
                                               'Medium','High', 'High'))
}

{
  # Regression after removing some samples
  lower_bound <- quantile(merge_data$brwe, 0.05, na.rm = TRUE)
  upper_bound <- quantile(merge_data$brwe, 0.95, na.rm = TRUE)
  
  rm_merge_data <- merge_data[merge_data$brwe >= lower_bound & merge_data$brwe <= upper_bound, ]
  imp_ini <- mice(rm_merge_data, m = 1, print=F, maxit = 1, method = "pmm")
  pred <- imp_ini$pred
  pred
  pred[, 1:502] <- 0
  meth <- imp_ini$meth
  meth
  meth[c('chsex_new', 'yzhdc_FOLIC1', 'yzhdc_SMKN3', 'fmbs_CBH46')] <- 'logreg'
  meth[c('yqyzdc_educton', 'yqyzdc_income', 'yqyzdc_folv4_1')] <- 'polr'
  
  imp <- mice(rm_merge_data, meth = meth, pred=pred, print=F, m=10, seed=2024)
  
  compound_names <- names(rm_merge_data)[2:502]
  confounder_names <- names(rm_merge_data)[c(503, 508,512)]
  
  pol_sheet <- data.frame()
  for (compound_var in compound_names) {
    formula <- paste("brwe ~", paste(c(compound_var, confounder_names), collapse = " + "))
    fit <- with(imp, lm(as.formula(formula)))
    pooled_result <- pool(fit)
    pol <- summary(pooled_result, conf.int=TRUE) %>%
      filter(term==compound_var)
    pol_sheet <- rbind(pol_sheet, pol)
  }
  
  pol_sheet <- pol_sheet %>%
    mutate(adjusted_p = p.adjust(p.value, method='BH')) %>%
    inner_join(cpd_info %>% select(unique_features, Name), by = c('term' = 'unique_features'))
  write.csv(pol_sheet, './20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe_remove_extreme.csv', row.names = FALSE) # Linear correlation compounds
}

{
  # Plot of significant compounds in different models (Sensitivity analysis)
  model_1_data <- read.csv('./20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe.csv', check.names = FALSE) %>%
    mutate(Model = 'Model_1')
  model_2_data <- read.csv('./20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe_no_adjust_bmi.csv', check.names = FALSE) %>%
    mutate(Model = 'Model_2')
  model_3_data <- read.csv('./20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe_no_adjust_yunzhou.csv', check.names = FALSE) %>%
    mutate(Model = 'Model_3')
  model_4_data <- read.csv('./20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe_no_adjustment.csv', check.names = FALSE) %>%
    mutate(Model = 'Model_4')
  model_5_data <- read.csv('./20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe_remove_extreme.csv', check.names = FALSE) %>%
    mutate(Model = 'Model_5')
  
  sig_cpds <- model_1_data[model_1_data[['p.value']] < 0.05, ]$term
  
  plot_data <- rbind(
    model_1_data[model_1_data$term %in% sig_cpds, ],
    model_2_data[model_2_data$term %in% sig_cpds, ],
    model_3_data[model_3_data$term %in% sig_cpds, ],
    model_4_data[model_4_data$term %in% sig_cpds, ],
    model_5_data[model_5_data$term %in% sig_cpds, ]
  ) %>%
    inner_join(cpd_info %>% select(Name, Class), by=c("Name" = "Name")) %>%
    group_by(Name) %>%  # Group by the compound
    mutate(y = ifelse(row_number() == 1, 1, 0)) %>% # Add 'y', 1 for the first row, 0 for others
    ungroup() %>%
    mutate(Significance = ifelse(p.value < 0.05, "Significant", "Not Significant")) %>%
    arrange(estimate, Class) %>%
    mutate(Name = factor(Name, levels = unique(Name)))
  
  ggplot(plot_data, aes(x = Name, y = estimate, shape = Model, color = Significance)) +
    geom_pointrange(aes(ymin =`2.5 %`, ymax = `97.5 %`), 
                    position = position_dodge(width = 0.7)) +  # Separate dependent variables to avoid overlap
    geom_point(size = 0.01, position = position_dodge(width = 0.7)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("Not Significant" = "black", "Significant" = "#7570b3")) +
    scale_shape_manual(values = c(16, 17, 18, 15, 14)) + 
    # coord_flip() +
    theme_minimal() +
    labs(title = "",
         y = "Beta Estimate (95% CI)", 
         x = "") +
    theme(panel.grid.major.x = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.line.x = element_blank(),
          axis.text.x = element_text(color = "black", size=10, angle=45, hjust=1, vjust=1),
          axis.title.y = element_text(color = "black", size=9),
          axis.title.x = element_blank(),
          strip.text = element_blank(),
          legend.position = 'bottom')
  ggsave("../article/pictures/20241107_reorganize/brwe/sentivity_brwe.png",
         width = 15, height=8, units = 'in', dpi = 300, bg='white')
    
}