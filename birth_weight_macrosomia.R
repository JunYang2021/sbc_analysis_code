library(tidyverse)
library(readxl)
library(patchwork)


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
  # Construct lm model to get variance explained (for covarites)
  library(caret)
  influence_factors <- c('chsex_new', 'yqyzdc_educton', 'yqyzdc_income',
                         'yqyzdc_folv4_1', 'maternal_age', 'before_bmi', 'yzhdc_FOLIC1',
                         'yzhdc_SMKN3','fmbs_CBH46', 'yunzhou')
  
  train_control <- trainControl(method = "repeatedcv",      # Cross-validation method
                                number = 10,                # Number of folds
                                repeats = 100,              # Number of repeats
                                savePredictions = TRUE)
  
  results <- data.frame(Factor = character(),
                        DependentVar = character(),
                        TrainR2 = numeric(),
                        ValidR2 = numeric(),
                        AdjR2 = numeric(),
                        pvalue=numeric())
  
  for (cur_fac in influence_factors) {
    print(cur_fac)
    # Perform which model (only one IV)
    merge_data_com <- merge_data %>%
      filter(!is.na(get(cur_fac)))
    formula <- as.formula(paste("brwe ~", cur_fac))
    cv_model <- train(formula, data = merge_data_com, 
                      method = "lm", trControl = train_control)
    valid_r2 <- mean(cv_model$resample$Rsquared)
    
    # Fit the model on the full training data (non-CV)
    full_model <- lm(formula, data = merge_data_com)
    print(summary(full_model))
    
    # Calculate R2 on training set
    train_predictions <- predict(full_model, merge_data_com)
    train_actuals <- merge_data_com$brwe
    ss_total <- sum((train_actuals - mean(train_actuals))^2)
    ss_residual <- sum((train_actuals - train_predictions)^2)
    train_r2 <- 1 - (ss_residual / ss_total)
    
    p_value <- coef(summary(full_model))[2, "Pr(>|t|)"]
    
    n <- nrow(merge_data_com)
    k <- 1
    adj_r2 <- 1 - (1 - valid_r2) * (n - 1) / (n - k - 1)
    
    # Save results
    results <- rbind(results, data.frame(Factor = cur_fac, 
                                         DependentVar = "brwe", 
                                         TrainR2 = train_r2, 
                                         ValidR2 = valid_r2, 
                                         AdjR2 = adj_r2,
                                         pvalue=p_value))
  }
  # Significant covarites: "chsex_new", "maternal_age", "before_bmi", "yunzhou", "fmbs_CBH46"
  write.csv(results, "./20241107_reanalysis/brwe_modeling/adjusted_r2_results.csv",
            row.names = FALSE)
}


{
  # Construct lm model to get correlations (for covarites)
  influence_factors <- c('chsex_new', 'yqyzdc_educton', 'yqyzdc_income',
                         'yqyzdc_folv4_1', 'maternal_age', 'before_bmi', 'yzhdc_FOLIC1',
                         'yzhdc_SMKN3','fmbs_CBH46', 'yunzhou')
  results_list <- list()
  for (cur_fac in influence_factors) {
    print(cur_fac)
    # Perform which model (only one IV)
    merge_data_com <- merge_data %>%
      filter(!is.na(get(cur_fac)))
    formula <- as.formula(paste("brwe ~", cur_fac))
    full_model <- lm(formula, data = merge_data_com)
    print(summary(full_model))
    model_summary <- summary(full_model)
    coefficients <- model_summary$coefficients
    conf_int <- confint(full_model)
    
    # Remove intercept row
    coefficients <- coefficients[-1, , drop = FALSE]
    conf_int <- conf_int[-1, , drop = FALSE]
    
    # Store results
    results_list[[cur_fac]] <- data.frame(
      Variable = cur_fac,
      Term = rownames(coefficients),
      Estimate = coefficients[, 1],
      Conf.Low = conf_int[, 1],
      Conf.High = conf_int[, 2],
      P.Value = coefficients[, 4],
      N = nrow(merge_data_com)
    )
    
  }
  results <- do.call(rbind, results_list) %>%
    mutate(combined = paste0(sprintf("%.3f", Estimate), " (", sprintf("%.3f", Conf.Low), ", ", sprintf("%.3f", Conf.High), ")"))
  # Significant covarites: "chsex_new", "maternal_age", "before_bmi", "yunzhou", "Maternal education"
  write.csv(results, "./20241107_reanalysis/brwe_modeling/covariates_brwe_lm.csv",
            row.names = FALSE)
}

{
  # Perform linear regression on birth weight in Batch 2
  library(mice)
  imp_ini <- mice(merge_data, m = 1, print=F, maxit = 1, method = "pmm")
  pred <- imp_ini$pred
  pred
  pred[, 1:500] <- 0
  meth <- imp_ini$meth
  meth
  meth[c('chsex_new', 'yzhdc_FOLIC1', 'yzhdc_SMKN3', 'fmbs_CBH46')] <- 'logreg'
  meth[c('yqyzdc_educton', 'yqyzdc_income', 'yqyzdc_folv4_1')] <- 'polr'
  
  imp <- mice(merge_data, meth = meth, pred=pred, print=F, m=10, seed=2024)
  
  compound_names <- names(merge_data)[2:500]
  confounder_names <- names(merge_data)[c(501, 506,510)]
  
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
  write.csv(pol_sheet, './20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe.csv', row.names = FALSE) # Linear correlation compounds
  
}

{
  #multiple linear regression on covariates (adjusted)
  confounder_names <- names(merge_data)[c(503:512)]

  formula <- paste("brwe ~", paste(confounder_names, collapse = " + "))
  fit <- with(imp, lm(as.formula(formula)))
  pooled_result <- pool(fit)
  pol <- summary(pooled_result, conf.int=TRUE)[-1, , drop = FALSE] %>%
    select(term, estimate, `2.5 %`, `97.5 %`, p.value)

  write.csv(pol, "./20241107_reanalysis/brwe_modeling/covariates_brwe_lm_adjusted.csv",
            row.names = FALSE)
}

{
  # Perform linear regression on birth weight in Batch 2 （No adjustment of BMI）
  compound_names <- names(merge_data)[2:500]
  confounder_names <- names(merge_data)[c(501, 510)]
  
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
  write.csv(pol_sheet, './20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe_no_adjust_bmi.csv', row.names = FALSE) # Linear correlation compounds
  
  # No adjustement of gestational age
  compound_names <- names(merge_data)[2:500]
  confounder_names <- names(merge_data)[c(501, 506)]
  
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
  write.csv(pol_sheet, './20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe_no_adjust_yunzhou.csv', row.names = FALSE)
  
  # No adjustment
  compound_names <- names(merge_data)[2:500]
  
  pol_sheet <- data.frame()
  for (compound_var in compound_names) {
    formula <- paste("brwe ~", compound_var)
    fit <- with(imp, lm(as.formula(formula)))
    pooled_result <- pool(fit)
    pol <- summary(pooled_result, conf.int=TRUE) %>%
      filter(term==compound_var)
    pol_sheet <- rbind(pol_sheet, pol)
  }
  
  pol_sheet <- pol_sheet %>%
    mutate(adjusted_p = p.adjust(p.value, method='BH')) %>%
    inner_join(cpd_info %>% select(unique_features, Name), by = c('term' = 'unique_features'))
  write.csv(pol_sheet, './20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe_no_adjustment.csv', row.names = FALSE) # Linear correlation compounds
  
  # UpSet Plot  用313新电脑（R: 4.3.2; ComplexUpset: 1.3.3; ggplot2: 3.5.1）
  cpd_info <- read_excel('../DATA/statistics/identified_compounds_information_batch2.xlsx', sheet=1)
  model.1 <- read.csv('./20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe.csv') %>%
    filter(p.value < 0.05)
  model.2 <- read.csv('./20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe_no_adjust_bmi.csv') %>%
    filter(p.value < 0.05)
  model.3 <- read.csv('./20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe_no_adjust_yunzhou.csv') %>%
    filter(p.value < 0.05)
  model.4 <- read.csv('./20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe_no_adjustment.csv') %>%
    filter(p.value < 0.05)
  
  upset_cpd <- cpd_info %>%
    select(Name, Class, Class_pie) %>%
    mutate(Model_1 = if_else(Name %in% model.1$Name, TRUE, FALSE),
           Model_2 = if_else(Name %in% model.2$Name, TRUE, FALSE),
           Model_3 = if_else(Name %in% model.3$Name, TRUE, FALSE),
           Model_4 = if_else(Name %in% model.4$Name, TRUE, FALSE)) %>%
    filter(Model_1 == TRUE | Model_2 == TRUE |Model_3 == TRUE | Model_4 == TRUE) %>%
    mutate(Class_pie= factor(Class_pie,
                             levels = c('Fatty acids',
                                        'Amino acids',
                                        'Bile acids',
                                        'Carnitines',
                                        
                                        'Fatty acyls',
                                        'Xenobiotics',
                                        'Lipids',
                                        
                                        'Steroids',
                                        
                                        'Organic acids',
                                        'Nucleotides',
                                        'Others'
                             )))
  
  cols <-  colnames(upset_cpd)[4:7]
  
  upset(upset_cpd, 
        cols,
        mode='inclusive_intersection',
        base_annotations=list(
          'Intersection size'=intersection_size(
            mode='inclusive_intersection',
            counts=TRUE,
            mapping=aes(fill=Class_pie)
          ) + 
            scale_fill_brewer(palette='Set3') +
            labs(fill = 'Class') 
        ),
        set_sizes=FALSE,
        sort_intersections=FALSE,
        intersections=list(
          'Model_1',
          'Model_2',
          'Model_3',
          'Model_4',
          c('Model_1', 'Model_2'),
          c('Model_1', 'Model_3'),
          c('Model_1', 'Model_4'),
          c('Model_2', 'Model_3'),
          c('Model_2', 'Model_4'),
          c('Model_3', 'Model_4'),
          c('Model_1', 'Model_2', 'Model_3','Model_4')),
        sort_sets=FALSE,
        name='Significant compounds (p < 0.05) in different models'
  )
  ggsave('../article/pictures/20241107_reorganize/brwe/upset_brwe_models_comparison.png', width = 6, height = 5, units = 'in',
         dpi=300)
}

{
  # Dot-whisper plot of different models
  cpd_info <- read_excel('../DATA/statistics/identified_compounds_information_batch2.xlsx', sheet=1)
  
  model_3_data <- read.csv('./20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe.csv', check.names = FALSE) %>%
    # filter(p.value < 0.05) %>%
    mutate(Model = 'Model_3')
  model_2_data <- read.csv('./20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe_no_adjust_bmi.csv', check.names = FALSE) %>%
    # filter(term %in% model_3_data$term) %>%
    mutate(Model = 'Model_2')
  model_1_data <- read.csv('./20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe_no_adjustment.csv', check.names = FALSE) %>%
    # filter(term %in% model_3_data$term) %>%
    mutate(Model = 'Model_1')
  
  class_palette <- c("Amino acids" = '#a6cee3',
                     "Bile acids" = '#1f78b4',
                     "Carnitines" = '#b2df8a',
                     "Fatty acids" = '#33a02c',
                     "Fatty acyls" = '#fb9a99',
                     "PA" = '#d9d9d9',
                     "PC" = '#003c30',
                     "PE" = '#ffff99',
                     "PG"= '#8e0152',
                     "PI" = '#313695',
                     "Nucleotides" = '#e31a1c',
                     "Organic acids" = '#fdbf6f',
                     "Others" = '#ff7f00',
                     "Sphingolipids" = '#cab2d6',
                     "Steroids" = '#6a3d9a',
                     "Xenobiotics"  = '#b15928')
  
  plot_data <- rbind(
    model_1_data, model_2_data, model_3_data
  ) %>%
    inner_join(cpd_info %>% select(Name, Class), by=c("Name" = "Name")) %>%
    group_by(Name) %>%  # Group by the compound
    mutate(y = ifelse(row_number() == 1, 1, 0)) %>% # Add 'y', 1 for the first row, 0 for others
    ungroup() %>%
    mutate(Significance = ifelse(p.value < 0.05, "Significant", "Not Significant")) %>%
    arrange(Class, estimate) %>%
    mutate(Name = factor(Name, levels = unique(Name))) %>%
    mutate(Class = if_else(Class == 'Glycerophosphates', 'PA', Class),
           Class = if_else(Class == 'Glycerophosphocholines', 'PC', Class),
           Class = if_else(Class == 'Glycerophosphoethanolamines', 'PE', Class),
           Class = if_else(Class == 'Glycerophosphoglycerols', 'PG', Class),
           Class = if_else(Class == 'Glycerophosphoinositols', 'PI', Class))
  
  plot_function <- function(class_list) {
    plot_data_c <- plot_data %>%
      filter(Class %in% class_list)
    dtplot <- ggplot(plot_data_c, aes(x = Name, y = estimate, shape = Model, color = Significance)) +
    geom_pointrange(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                    position = position_dodge(width = 0.3)) +  # Separate dependent variables to avoid overlap
    geom_point(size = 0.01, position = position_dodge(width = 0.7)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("Not Significant" = "black", "Significant" = "#7570b3")) +
    scale_shape_manual(values = c(16, 17, 18)) + 
    # coord_flip() +
    theme_minimal() +
    labs(title = "",
         y = "", 
         x = "Beta Estimate (95% CI)") +
    theme(panel.grid.major.x = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.line.x = element_blank(),
          axis.text.x = element_blank(), 
          axis.title.x = element_blank(),
          strip.text = element_blank())
  
  p_y <-
    ggplot(plot_data_c) +
    geom_col(aes(x = Name, y = y, fill = Class)) +
    scale_fill_manual(values = class_palette) +
    coord_cartesian(expand = F) +
    # coord_flip() +
    theme(axis.text.y = element_blank(), 
          axis.text.x = element_text(color = "black", size=10, angle=45, hjust=1, vjust=1),
          axis.title.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.background = element_blank(), 
          plot.background = element_blank())
  
  library(patchwork)
  dtplot + p_y +
    plot_layout(heights = c(1, .07), guides='collect') &
    theme(legend.position = "bottom")
  }
  
  plot_function(c('Fatty acids'))
  plot_function(c("Carnitines",'Fatty acyls'))
  ggsave("../article/pictures/20241107_reorganize/brwe/dot_whisker_plot_three_models_car.png", 
         width = 15, height =6.6, units = c("in"), dpi = 300, bg='white')
  plot_function(c('PC'))
  ggsave("../article/pictures/20241107_reorganize/brwe/dot_whisker_plot_three_models_pc.png", 
         width = 15, height =6.6, units = c("in"), dpi = 300, bg='white')
  plot_function(c('PE', 'PI', 'PG', 'PA'))
  ggsave("../article/pictures/20241107_reorganize/brwe/dot_whisker_plot_three_models_pe.png", 
         width = 15, height =6.6, units = c("in"), dpi = 300, bg='white')
  plot_function(c("Amino acids","Organic acids"))
  ggsave("../article/pictures/20241107_reorganize/brwe/dot_whisker_plot_three_models_aa.png", 
         width = 15, height =6.6, units = c("in"), dpi = 300, bg='white')
  plot_function(c("Bile acids", "Steroids"))
  ggsave("../article/pictures/20241107_reorganize/brwe/dot_whisker_plot_three_models_ba.png", 
         width = 15, height =6.6, units = c("in"), dpi = 300, bg='white')
  plot_function(c("Nucleotides","Sphingolipids", "Xenobiotics"))
  ggsave("../article/pictures/20241107_reorganize/brwe/dot_whisker_plot_three_models_nucle.png", 
         width = 15, height =6.6, units = c("in"), dpi = 300, bg='white')
  plot_function(c("Others"))
  ggsave("../article/pictures/20241107_reorganize/brwe/dot_whisker_plot_three_models_other.png", 
         width = 15, height =6.6, units = c("in"), dpi = 300, bg='white')
  
  
  # Plot for FA displaying trend of C atoms
  plot_data_c <- plot_data %>%
    filter(Class == 'Fatty acids') %>%
    inner_join(read_excel('../DATA/statistics/FA_structures.xlsx'), by = c('Name' = 'FA_name')) %>%
    filter(sum_comp != '_' & keep==1) %>%
    mutate(sum_comp = factor(sum_comp, levels = unique(sum_comp)))
  dtplot <- ggplot(plot_data_c, aes(x = sum_comp, y = estimate, shape = Model, color = Significance)) +
    geom_pointrange(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                    position = position_dodge(width = 0.3)) +  # Separate dependent variables to avoid overlap
    geom_point(size = 0.01, position = position_dodge(width = 0.7)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("Not Significant" = "black", "Significant" = "#7570b3")) +
    scale_shape_manual(values = c(16, 17, 18)) + 
    theme_minimal() +
    labs(title = "",
         y = "", 
         x = "Beta Estimate (95% CI)") +
    theme(panel.grid.major.x = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.line.x = element_blank(),
          # axis.text.x = element_blank(),
          axis.text.x = element_text(color = "black", size=8, angle=90,hjust=1.2, vjust=0.8),
          axis.title.x = element_blank(),
          strip.text = element_blank())
  
  p_y <-
    ggplot(plot_data_c) +
    geom_point(aes(x = sum_comp, y = DB_number, color = "#Double bonds")) +
    geom_line(aes(x = sum_comp, y = DB_number, group = 1, color = "#Double bonds"), linewidth=0.6) +
    geom_point(aes(x =sum_comp, y = C_number, color = "#Carbon atoms")) +
    geom_line(aes(x = sum_comp, y = C_number, group = 1, color = "#Carbon atoms"), linewidth=0.6) +
    scale_color_manual(values = c("#Double bonds" = "#d7191c", 
                                  "#Carbon atoms" = "#2c7bb6")) +
    coord_cartesian(expand = F) +
    theme(axis.text.y = element_blank(), 
          axis.text.x = element_text(color = "black", size=10, angle=45, hjust=1, vjust=1),
          axis.title.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.background = element_blank(), 
          plot.background = element_blank()) +
    labs(color = "Structure")

  dtplot + p_y +
    plot_layout(heights = c(1, .2), guides='collect') &
    theme(legend.position = "bottom")
  
  ggsave("../article/pictures/20241107_reorganize/brwe/dot_whisker_plot_three_models_FA.png", 
         width = 16, height =6.6, units = c("in"), dpi = 300, bg='white')
  
  dtplot
  ggsave("../article/pictures/20241107_reorganize/brwe/dot_whisker_plot_three_models_no_number.png", 
         width = 10, height =4, units = c("in"), dpi = 300, bg='white')
  
}

{
  # Plot of the significant compounds with FDR filtering
  batch2_sig <- read.csv('./20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe.csv', check.names = FALSE) 
  filtered_data <- batch2_sig %>%
    filter(adjusted_p < 0.1) %>%
    mutate(effect = 1)   # Negative estimated
  class_palette <- c("Amino acids" = '#a6cee3',
                     "Bile acids" = '#1f78b4',
                     "Carnitines" = '#b2df8a',
                     "Fatty acids" = '#33a02c',
                     "Fatty acyls" = '#fb9a99',
                     "PA" = '#d9d9d9',
                     "PC" = '#003c30',
                     "PE" = '#ffff99',
                     "PG"= '#8e0152',
                     "PI" = '#313695',
                     "Nucleotides" = '#e31a1c',
                     "Organic acids" = '#fdbf6f',
                     "Others" = '#ff7f00',
                     "Sphingolipids" = '#cab2d6',
                     "Steroids" = '#6a3d9a',
                     "Xenobiotics"  = '#b15928')
  filtered_data[filtered_data$Name %in% c('LPE 20:4',
                                          'LPG O-13:1',
                                          'Carnitine',
                                          'LPE 18:2 (isomer)',
                                          'Leucyl-Proline',
                                          'Deoxycholic acid',
                                          'Acetylcholine',
                                          'Leucylphenylalanine'), 'effect'] <- -1
  filtered_data[6, 'Name'] <- 'Glycerophosphocholine' # "Glycerophosphocholine (isomer)"
  filtered_data[7, 'Name'] <- 'LPE 18:2' # "LPE 18:2 (isomer)"
  
  plot_data <- filtered_data %>%
    inner_join(cpd_info %>% select(Name, Class), by=c("Name" = "Name")) %>%
    mutate(Name = if_else(Name== "Perfluoroundecanoic acid", 'PFUnDA', Name)) %>%
    mutate(Class = case_when(
      Class == "Glycerophosphates" ~ "PA",
      Class == "Glycerophosphocholines" ~ "PC",
      Class == "Glycerophosphoethanolamines" ~ "PE",
      Class == "Glycerophosphoglycerols" ~ "PG",
      Class == "Glycerophosphoinositols" ~ "PI",
      TRUE ~ Class
    )) %>%
    arrange(effect, Class) %>%# First sort by 'effect' and then by 'Class'
    mutate(Name = factor(Name, levels = unique(Name)))
  
  ggplot(plot_data, aes(x = Name, y = estimate)) +
    geom_pointrange(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                    position = position_dodge(width = 0.6)) +  # Separate dependent variables to avoid overlap
    geom_point(size = 0.01, position = position_dodge(width = 0.6)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    theme_classic() +
    labs(title = "",
         x = "", 
         y = "Beta Estimate (95% CI)") +
    theme(
          axis.text.x = element_text(color = "black", size=10, angle=45, hjust=1, vjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color = "black", size=9))
  
  ggsave("../article/pictures/20241107_reorganize/brwe/dot_whisker_plot_batch2_FDR_0_1_sig_cpds.png",
         width = 8, height =4.5, units = c("in"), dpi = 300, bg='white')
}

{
  # Output for pathway analysis
  cpd_info <- read_excel('../DATA/statistics/identified_compounds_information_batch2.xlsx', sheet=1)
  sig_cpd <- read.csv('./20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe.csv', check.names = FALSE) %>%
    filter(p.value < 0.05) %>%
    inner_join(cpd_info %>% select(Name, HMDB), by='Name') %>%
    mutate(Effect = if_else(estimate > 0, 'Up', 'Down')) %>%
    filter(Effect=='Up')
  
  for (id in sig_cpd$HMDB) {
    cat(id, "\n")
  }  # Metabanalyst: metabolites, results: ../article/pictures/20241107_reorganize/pathway_analysis/brwe_up_metabolites
  
  sig_cpd <- read.csv('./20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe.csv', check.names = FALSE) %>%
    filter(p.value < 0.05) %>%
    inner_join(cpd_info %>% select(Name, HMDB), by='Name') %>%
    mutate(Effect = if_else(estimate > 0, 'Up', 'Down')) %>%
    filter(Effect=='Down')
  
  for (id in sig_cpd$HMDB) {
    cat(id, "\n")
  } # 5个通路中有三个是因为valine，没有参考价值
  
  # Remake the dot plot
  enrich_res <- read.csv('../article/pictures/20241107_reorganize/pathway_analysis/brwe_up_metabolites/msea_ora_result.csv') %>%
    mutate(enrichment_ratio = hits / expected) %>%
    mutate(logP = -log10(Raw.p)) %>%
    mutate(X = reorder(X, logP)) 
  
  ggplot(enrich_res, aes(x = logP, y = X)) +
    geom_point(aes(size = enrichment_ratio, color = Raw.p)) + 
    scale_color_gradient(low = "red", high = "yellow") +  # Gradient color for p-values
    scale_size_continuous(range = c(2, 10)) +            # Adjust point sizes
    theme_minimal() +                                    # Clean theme
    labs(
      x = expression(-log[10](p-value)), 
      y = NULL,
      color = "P-value",
      size = "Enrichment Ratio"
    ) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.y = element_text(size = 9, hjust = 1, color='black'),
      legend.position = c(.75, .37)
    )
  
  ggsave('../article/pictures/20241107_reorganize/brwe/batch2_sig_up_compounds_brwe_pathway_kegg.png', 
         width =5.55, height =4.95, units = c("in"), dpi = 300, bg='white')
  
  
}


{
  # 画出68个显著代谢物的类别分布和趋势分布
  batch2_sig <- read.csv('./20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe.csv', check.names = FALSE) %>%
    filter(p.value < 0.05) %>%
    inner_join(cpd_info %>% select(Name, Class), by=c("Name" = "Name")) %>%
    mutate(Effect = if_else(estimate > 0, 'Positive', 'Negative')) %>%
    mutate(Effect=factor(Effect, levels=c('Positive', 'Negative')))
    
  
  summary_counts <- batch2_sig %>%
    mutate(Class = case_when(
      Class == "Glycerophosphates" ~ "PA",
      Class == "Glycerophosphocholines" ~ "PC",
      Class == "Glycerophosphoethanolamines" ~ "PE",
      Class == "Glycerophosphoglycerols" ~ "PG",
      Class == "Glycerophosphoinositols" ~ "PI",
      TRUE ~ Class
    )) %>%
    mutate(Class=factor(Class, levels=c("Fatty acids",
                                        "Fatty acyls",
                                        "PC",
                                        "PE",
                                        "PA",
                                        "PI",
                                        "PG",
                                        "Amino acids",
                                        "Steroids",
                                        "Carnitines",
                                        "Bile acids",
                                        "Organic acids",
                                        "Xenobiotics",
                                        "Others"))) %>%
    group_by(Class, Effect) %>%
    summarise(count = n()) %>%
    arrange(Class, desc(Effect)) %>%
    mutate(lab_ypos = cumsum(count) - 0.5 * count)
  
  ggplot(summary_counts, aes(x = Class, y = count, fill = Effect)) +
    geom_bar(stat = "identity") +
    labs(
      x = "",
      y = "Count"
    ) +
    geom_text(
      aes(y = lab_ypos, label = count, group = Effect),
      color = "white",
      size=2.5
    ) + 
    scale_fill_manual(values = c("#e41a1c", "#377eb8")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = 'black'),
          legend.position = c(.5, .85))
  ggsave("../article/pictures/20241107_reorganize/brwe/batch2_sig_compounds_brwe.png", 
         width = 3.7, height =3.3, units = c("in"), dpi = 300, bg='white')
  

}

{
  # 代谢物对体重的variance explained (elastic net)
  {
    # Construct machine learning model to get variance explained
    library(caret)
    merge_data <- cpd_data %>%
      inner_join(birth_data, by='sample_name') %>%
      filter(!is.na(brwe))
    batch2_sig <- read.csv('./20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe.csv', check.names = FALSE) %>%
      filter(p.value < 0.05)
    
    train_control <- trainControl(method = "repeatedcv",      # Cross-validation method
                                  number = 10,                # Number of folds
                                  repeats = 100,              # Number of repeats
                                  savePredictions = TRUE)
    
    results <- read.csv("./20241107_reanalysis/brwe_modeling/adjusted_r2_results.csv")
    
    # Perform Lasso regression
    formula <- as.formula(paste("brwe ~", paste(batch2_sig$term, collapse = " + ")))
    cv_model <- train(formula,
                      data = merge_data, 
                      method = "glmnet", 
                      trControl = train_control)
    valid_r2 <- mean(cv_model$resample$Rsquared)  # 4.2 %
    
    n <- nrow(merge_data)
    k <- 76
    adj_r2 <- 1 - (1 - valid_r2) * (n - 1) / (n - k - 1)
    
    train_predictions <- predict(cv_model, newdata = merge_data)
    
    # Actual values
    y_actual <- merge_data$brwe
    
    # Compute R-squared for the training set
    sst <- sum((y_actual - mean(y_actual))^2)  # Total sum of squares
    sse <- sum((y_actual - train_predictions)^2)  # Residual sum of squares
    train_r2 <- 1 - (sse / sst) # 7.6 % 
  }
}

{
  # 性别的分层分析
  # Run Line 5->23 to get merge_data
  boy_merge_data <- merge_data %>%
    filter(chsex_new == 1)   # Number: 569
  girl_merge_data <- merge_data %>%
    filter(chsex_new == 2)  # Number: 556
  
  imp_ini <- mice(boy_merge_data, m = 1, print=F, maxit = 1, method = "pmm")
  pred <- imp_ini$pred
  pred
  pred[, 1:500] <- 0
  meth <- imp_ini$meth
  meth
  meth[c('yzhdc_FOLIC1', 'yzhdc_SMKN3', 'fmbs_CBH46')] <- 'logreg'
  meth[c('yqyzdc_educton', 'yqyzdc_income', 'yqyzdc_folv4_1')] <- 'polr'
  
  imp_boy <- mice(boy_merge_data, meth = meth, pred=pred, print=F, m=10, seed=2026)
  plot(imp_boy)
  stripplot(imp_boy, yqyzdc_educton~.imp, pch=5, cex=1)
  
  imp_girl <- mice(girl_merge_data, meth = meth, pred=pred, print=F, m=10, seed=2026)
  plot(imp_girl)
  
  compound_names <- names(merge_data)[2:500]
  confounder_names <- names(merge_data)[c(506, 510)]
  
  boy_pol_sheet <- data.frame()
  for (compound_var in compound_names) {
    formula <- paste("brwe ~", paste(c(compound_var, confounder_names), collapse = " + "))
    fit <- with(imp_boy, lm(as.formula(formula)))
    pooled_result <- pool(fit)
    pol <- summary(pooled_result, conf.int=TRUE) %>%
      filter(term==compound_var)
    boy_pol_sheet <- rbind(boy_pol_sheet, pol)
  }
  
  boy_pol_sheet <- boy_pol_sheet %>%
    mutate(adjusted_p = p.adjust(p.value, method='BH')) %>%
    inner_join(cpd_info %>% select(unique_features, Name), by = c('term' = 'unique_features'))
  
  girl_pol_sheet <- data.frame()
  for (compound_var in compound_names) {
    formula <- paste("brwe ~", paste(c(compound_var, confounder_names), collapse = " + "))
    fit <- with(imp_girl, lm(as.formula(formula)))
    pooled_result <- pool(fit)
    pol <- summary(pooled_result, conf.int=TRUE) %>%
      filter(term==compound_var)
    girl_pol_sheet <- rbind(girl_pol_sheet, pol)
  }
  
  girl_pol_sheet <- girl_pol_sheet %>%
    mutate(adjusted_p = p.adjust(p.value, method='BH')) %>%
    inner_join(cpd_info %>% select(unique_features, Name), by = c('term' = 'unique_features'))
  
  boy_girl_str_sheet <- bind_rows(
    boy_pol_sheet %>% mutate(`Fetal sex` = 'Boy'),
    girl_pol_sheet %>% mutate(`Fetal sex` = 'Girl')
  )
  
  write.csv(boy_girl_str_sheet, 
            file="./20241107_reanalysis/brwe_modeling/multiple_regression_batch2_gender_stratification_brwe.csv", row.names = FALSE)
  
  
  # Plot scatter plots of stratified analysis
  library(ggrepel)
  boy_girl_str_sheet <- read.csv("./20241107_reanalysis/brwe_modeling/multiple_regression_batch2_gender_stratification_brwe.csv",
                                 check.names = FALSE)
  data_stratified <- boy_girl_str_sheet %>%
    mutate(z = estimate/std.error) %>%
    select(term, Name, p.value, z, `Fetal sex`) %>%
    pivot_wider(
      names_from = `Fetal sex`,
      values_from = c(p.value, z),
      names_glue = "{`Fetal sex`} {.value}"
    ) %>%
    mutate(
      Shape = case_when(
        `Boy p.value` < 0.05 & `Girl p.value` < 0.05 ~ "Both significant",
        `Boy p.value` < 0.05 ~ "Significant in males",
        `Girl p.value` < 0.05 ~ "Significant in females",
        TRUE ~ "Not significant"
      )
    ) %>%
    mutate(
      Quadrant = factor(case_when(
        `Boy z` > 0 & `Girl z` > 0 ~ 'I',
        `Boy z` < 0 & `Girl z` > 0 ~ 'II',
        `Boy z` < 0 & `Girl z` < 0 ~ 'III',
        `Boy z` > 0 & `Girl z` < 0 ~ 'IV',
      ))
    ) %>%
    filter(Shape != 'Not significant') %>%
    inner_join(cpd_info %>% select(Name, Class_pie), by='Name') %>%
    mutate(Class_pie = factor(Class_pie, levels = c('Amino acids',
                                                    'Bile acids',
                                                    'Carnitines',
                                                    'Fatty acids',
                                                    'Fatty acyls',
                                                    'Lipids',
                                                    'Steroids',
                                                    'Organic acids',
                                                    'Nucleotides',
                                                    'Xenobiotics',
                                                    'Others')))
  
  class_palette <- c("Amino acids" = '#a6cee3',
                     "Bile acids" = '#1f78b4',
                     "Carnitines" = '#b2df8a',
                     "Fatty acids" = '#33a02c',
                     "Fatty acyls" = '#fb9a99',
                     "Lipids" = '#cab2d6',
                     "Nucleotides" = '#e31a1c',
                     "Organic acids" = '#fdbf6f',
                     "Others" = '#ff7f00',
                     "Steroids" = '#6a3d9a',
                     "Xenobiotics"  = '#b15928')
  
  ggplot(data_stratified, aes(x = `Boy z`, y = `Girl z`, shape = Shape, color=Class_pie)) +
    geom_point(size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "gray") +
    geom_hline(yintercept = 0, color = "black") +
    geom_vline(xintercept = 0, color = "black") +
    geom_text_repel(aes(label = Name),size = 2, max.overlaps = 9, color='black',
                    min.segment.length = Inf) + # Add labels
    # geom_text(check_overlap = TRUE, color='black', size=2)+
    scale_shape_manual(values = c('Both significant' = 16,
                                  'Significant in males'=17,
                                  'Significant in females'=15,
                                  'Not significant'=18)) +
    scale_color_manual(values = class_palette) + 
    labs(
      x = "Z-scores of birth weight (male offspring)",
      y = "Z-scores of birth weight (female offspring)",
      shape = "Significance",
      color='Class'
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      text = element_text(size = 12),
      panel.border = element_blank()
    )
  print(summary(data_stratified$Quadrant))
  # I  II III  IV 
  # 43   5  17  17 
  ggsave("../article/pictures/20241107_reorganize/brwe/stratified_gender.png", 
         width = 6,
         height =5, 
         units = c("in"), 
         dpi = 300,
         bg='white')
  
  # class distribution
  fa_data_stratified <- data_stratified %>%
    filter(Class_pie == 'Fatty acids')
  
  ggplot(fa_data_stratified, aes(x = `Boy z`, y = `Girl z`, shape = Shape, color=Class_pie)) +
    geom_point(size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") + 
    geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 0, color = "black") +
    geom_vline(xintercept = 0, color = "black") +
    geom_text_repel(aes(label = Name), size = 3, max.overlaps = 10, color='black') + # Add labels
    scale_shape_manual(values = c(16, 17, 15, 18)) +
    scale_color_manual(
      values = c(
        'Amino acids'='#a6cee3',
        'Bile acids'='#1f78b4',
        'Carnitines'='#b2df8a',
        'Fatty acids'='#33a02c',
        'Fatty acyls'='#fb9a99',
        'Lipids'='#e31a1c',
        'Steroids'='#fdbf6f',
        'Organic acids'='#ff7f00',
        'Nucleotides'='#cab2d6',
        'Xenobiotics'='#6a3d9a',
        'Others'='#b15928'
      )) + 
    labs(
      x = "Z-scores of birth weight (boys offspring)",
      y = "Z-scores of birth weight (girls offspring)",
      shape = "Significance",
      color='Class'
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      text = element_text(size = 12),
      panel.border = element_blank()
    )
  ggsave("../article/pictures/20241107_reorganize/brwe/stratified_gender_fa.png", 
         width = 6,
         height =5, 
         units = c("in"), 
         dpi = 300,
         bg='white')
  
  
  # Gender interaction term
  gender_interaction_term <- data.frame(term=sapply(cpd_info$unique_features, function(x) paste0(x, ":chsex_new2")),
                                        term_name=sapply(cpd_info$Name, function(x) paste0(x, ":gender")))
  
  library(mice)
  imp_ini <- mice(merge_data, m = 1, print=F, maxit = 1, method = "pmm")
  pred <- imp_ini$pred
  pred
  pred[, 1:500] <- 0
  meth <- imp_ini$meth
  meth
  meth[c('chsex_new', 'yzhdc_FOLIC1', 'yzhdc_SMKN3', 'fmbs_CBH46')] <- 'logreg'
  meth[c('yqyzdc_educton', 'yqyzdc_income', 'yqyzdc_folv4_1')] <- 'polr'
  
  imp <- mice(merge_data, meth = meth, pred=pred, print=F, m=10, seed=2024)
  
  compound_names <- names(merge_data)[2:500]
  confounder_names <- names(merge_data)[c(501, 506, 510)]
  
  pol_sheet <- data.frame()
  for (compound_var in compound_names) {
    interaction_term <- paste(compound_var, '*', 'chsex_new')
    formula <- paste("brwe ~", paste(c(compound_var,interaction_term, confounder_names), collapse = " + "))
    fit <- with(imp, lm(as.formula(formula)))
    pooled_result <- pool(fit)
    pol <- summary(pooled_result, conf.int=TRUE) %>%
      filter(term==paste(compound_var, ":chsex_new2", sep = ''))
    pol_sheet <- rbind(pol_sheet, pol)
  }
  
  pol_sheet <- pol_sheet %>%
    mutate(adjusted_p = p.adjust(p.value, method='BH')) %>%
    inner_join(gender_interaction_term, by='term')
  write.csv(pol_sheet, 
            './20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe_gender_stratification.csv',
            row.names = FALSE) # Linear correlation compounds
  
  
  # Plot the dot-whisper plot
  sig_interacted <- read.csv('./20241107_reanalysis/brwe_modeling/mlr_compounds_batch2_brwe_gender_stratification.csv',
                             check.names = FALSE) %>%
    filter(p.value < 0.05) %>%
    mutate(feature = gsub(':chsex_new2', "", term))   # 不删除本来不显著的
  
  boy_girl_str_sheet <- read.csv("./20241107_reanalysis/brwe_modeling/multiple_regression_batch2_gender_stratification_brwe.csv",
                                 check.names = FALSE) %>%
    mutate(`Fetal sex` = ifelse(`Fetal sex` == "Boy", "Male",
                              ifelse(`Fetal sex` == "Girl", "Female", `Fetal sex`)))
  
  # Dot-whisper plot
  library(patchwork)
  class_palette <- c("Amino acids" = '#a6cee3',
                     "Bile acids" = '#1f78b4',
                     "Carnitines" = '#b2df8a',
                     "Fatty acids" = '#33a02c',
                     "Fatty acyls" = '#fb9a99',
                     "Glycerophosphates" = '#d9d9d9',
                     "Glycerophosphocholines" = '#003c30',
                     "PE" = '#ffff99',
                     "Glycerophosphoglycerols"= '#8e0152',
                     "Glycerophosphoinositols" = '#313695',
                     "Nucleotides" = '#e31a1c',
                     "Organic acids" = '#fdbf6f',
                     "Others" = '#ff7f00',
                     "Sphingolipids" = '#cab2d6',
                     "Steroids" = '#6a3d9a',
                     "Xenobiotics"  = '#b15928')
  
  plot_data <- boy_girl_str_sheet %>%
    filter(term %in% sig_interacted$feature) %>%
    mutate(Significance = ifelse(p.value < 0.05, "Significant", "Not Significant")) %>%
    inner_join(cpd_info %>% select(unique_features, Class), by=c("term" = "unique_features")) %>%
    group_by(Name) %>%  # Group by the compound
    mutate(effect = if_else(estimate[`Fetal sex` == "Male"] < estimate[`Fetal sex` == "Female"], -1, 1)) %>%
    mutate(y = ifelse(row_number() == 1, 1, 0)) %>% # Add 'y', 1 for the first row, 0 for others
    ungroup() %>%
    arrange(effect, Class, Name) %>%
    mutate(Name = if_else(Name == '(3S,5R,6R,7E)-3,5,6-Trihydroxy-7-megastigmen-9-one',
                          "CHEBI:168487", Name),
           Name = if_else(Name == '3beta,7alpha-Dihydroxy-5-cholestenoate',
                          "3,7-Dihydroxy-5-cholestenoic acid", Name),
           Name = if_else(Name == 'Lithocholic acid 3-O-sulfate',
                          "Sulfolithocholate", Name)) %>%
    mutate(Name = factor(Name, levels = unique(Name))) %>%
    mutate(Class = if_else(Class == 'Glycerophosphoethanolamines', 'PE', Class))
  
  dtplot <- ggplot(plot_data, aes(x = Name, y = estimate, shape = `Fetal sex`, color = Significance)) +
    geom_pointrange(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                    position = position_dodge(width = 0.3)) +  # Separate dependent variables to avoid overlap
    geom_point(size = 0.01, position = position_dodge(width = 0.3)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("Not Significant" = "black", "Significant" = "#7570b3")) +
    scale_shape_manual(values = c("Male"=16, "Female"=17)) + 
    theme_minimal() +
    labs(title = "",
         x = "", 
         y = "Beta Estimate (95% CI)") +
    theme(panel.grid.major.x = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color = "black", size=9),
          legend.position = 'bottom',
          legend.direction = "horizontal"
    )
  
  p_y <-
    ggplot(plot_data) +
    geom_col(aes(x = Name, y = y, fill = Class)) +
    scale_fill_manual(values = class_palette) +
    coord_cartesian(expand = F) +
    theme(axis.text.y = element_blank(), 
          axis.text.x = element_text(color = "black", size=10, angle=45, hjust=1, vjust=1),
          axis.title.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.background = element_blank(), 
          plot.background = element_blank())
  
  dtplot + p_y +
    plot_layout(heights = c(1, .15), guides='collect') &
    theme(legend.position = 'none')
  ggsave("../article/pictures/20241107_reorganize/brwe/dot_whisker_plot_batch2_gender_interacted_cpd_stra.png",
         width =4.63, height =5, units = c("in"), dpi = 300, bg='white')
  
  library(cowplot)
  legend <- cowplot::get_legend(dtplot)
  cowplot::plot_grid(legend, ncol = 1)
  ggsave("../article/pictures/20241107_reorganize/stratified/sex_stratification_legend.png", width = 6, height = 0.25)
  
}

{
  # 通过巨大儿不良结局研究风险 
  # Run Line 5->23 to get merge_data
  macro_data <- merge_data %>%
    filter(growth == 'Macrosomia')   # 由于与体重相关的变量中，只有maternal_age和孕前BMI存在少量缺失值，所以可以用均值代替一下
  merge_data$maternal_age[is.na(merge_data$maternal_age)] <- mean(merge_data$maternal_age, na.rm=TRUE)
  merge_data$before_bmi[is.na(merge_data$before_bmi)] <- mean(merge_data$before_bmi, na.rm=TRUE)
  merge_data <- merge_data %>%
    filter(!is.na(yunzhou))
  
  # Samples matching
  library(MatchIt)
  
  m.out0 <- matchit(growth ~ chsex_new + yunzhou + maternal_age + before_bmi, 
                    data = merge_data,
                    method=NULL, distance = 'glm') # 其他指标缺失值太多，不适合再考虑
  summary(m.out0)
  m.out1 <- matchit(growth ~ chsex_new + yunzhou + maternal_age + before_bmi, 
                    data = merge_data,
                    method='nearest', distance = 'glm', caliper = 0.1)
  m.out1
  summary(m.out1, un = FALSE)
  tiff(filename = "../article/pictures/20241107_reorganize/brwe/match_sample_brwe.tiff",
       width = 5, height = 5, units = "in", res = 300)
  
  plot(m.out1, type = "jitter", interactive = FALSE)
  dev.off()
  
  plot(m.out1, type = "density", interactive = FALSE,
       which.xs = ~chsex_new + yunzhou + maternal_age + before_bmi)
  plot(summary(m.out1))
  
  m.data <- match.data(m.out1)
  
  # Binary logistic regression
  metabolite_names <- names(m.data)[2:505]
  
  results <- data.frame(term = character(), 
                        Odds_Ratio = numeric(), 
                        ORCI_Lower = numeric(),
                        ORCI_Upper = numeric(),
                        p_Value = numeric(), 
                        stringsAsFactors = FALSE)
  for (metab_cpd in metabolite_names) {
    formula <- paste0('growth ~', metab_cpd, "+ chsex_new + yunzhou + maternal_age + before_bmi")
    model <- glm(as.formula(formula), data = m.data, family = binomial)
    
    coef <- coef(summary(model))[metab_cpd, "Estimate"]
    se <- coef(summary(model))[metab_cpd, "Std. Error"]
    z_value <- qnorm(0.975) # For 95% confidence interval
    ci_lower <- coef - z_value * se
    ci_upper <- coef + z_value * se
    
    odds_ratio <- exp(coef)
    orci_lower <- exp(ci_lower)
    orci_upper <- exp(ci_upper)
    
    # Extract p-value
    p_value <- coef(summary(model))[metab_cpd, "Pr(>|z|)"]
    
    results <- rbind(results, data.frame(term = metab_cpd, 
                                         Odds_Ratio = odds_ratio, 
                                         ORCI_Lower = orci_lower, 
                                         ORCI_Upper = orci_upper,
                                         p_Value = p_value))
  }
  results <- results %>%
    mutate(adjust_p = p.adjust(p_Value, method='BH')) %>%
    inner_join(cpd_info %>% select(unique_features, Name), by = c('term' = 'unique_features'))
  write.csv(results, './20241107_reanalysis/brwe_modeling/macrosomia_logistic.csv', row.names = FALSE)
  
  # Mann-Whitney U test
  results <- data.frame(term = character(), 
                        p_Value = numeric(), 
                        stringsAsFactors = FALSE)
  for (metab_cpd in metabolite_names) {
    test <- wilcox.test(m.data[[metab_cpd]] ~ m.data$growth)
    
    results <- rbind(results, data.frame(term = metab_cpd, 
                                         p_Value = test$p.value))
  }
  results <- results %>%
    mutate(adjust_p = p.adjust(p_Value, method='BH')) %>%
    inner_join(cpd_info %>% select(unique_features, Name), by = c('term' = 'unique_features'))
  write.csv(results, './20241107_reanalysis/brwe_modeling/macrosomia_u-test.csv', row.names = FALSE)
  
  log_result <- read.csv('./20241107_reanalysis/brwe_modeling/macrosomia_logistic.csv') %>%
    filter(p_Value < 0.05) %>% # Number: 73
    filter(adjust_p < 0.2) %>%  # Number: 23
    inner_join(cpd_info %>% select(Name, Class), by=c("Name" = "Name"))
  
  utest_result <- read.csv('./20241107_reanalysis/brwe_modeling/macrosomia_u-test.csv') %>%
    filter(p_Value < 0.05)  # Number: 65
  
  # Plot heatmap for log_result
  library(ComplexHeatmap)
  m_data <- m.data %>%
    mutate(growth = as.character(growth),               # Convert growth to character
           growth = if_else(growth == "Non macrosomia", "Control", growth)) %>%
    mutate(growth = factor(growth, levels = c("Control", "Macrosomia")))
  heatmap_data <- m.data %>%
    mutate(growth = as.character(growth),               # Convert growth to character
           growth = if_else(growth == "Non macrosomia", "Control", growth)) %>%
    mutate(growth = factor(growth, levels = c("Control", "Macrosomia"))) %>%
    select(log_result$term) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "rowname") %>% # Add row names as a column
    inner_join(cpd_info %>% select(unique_features, Name), by = c("rowname" = "unique_features")) %>%
    mutate(Name = if_else(Name == 'Glycerophosphocholine (isomer)', 'Glycerophosphocholine', Name)) %>%
    mutate(Name = if_else(Name == 'Octadecanedioate', 'FA 18:1;O2', Name)) %>%
    column_to_rownames(var='Name') %>%
    select(-rowname) %>%
    as.matrix()
  
  growth_group <- m_data$growth
  
  col_annotation <- HeatmapAnnotation(
    Growth = growth_group,
    col = list(Growth = c("Control" = "#377eb8", "Macrosomia" = "#e41a1c")), # Custom colors for binary groups
    show_legend = FALSE
  )
  
  # Generate the heatmap
  png("../article/pictures/20241107_reorganize/brwe/macrosomia_sample_matching_logistic.png",
      width=5, height=5.7, units="in", res=300)
  Heatmap(
    heatmap_data,
    top_annotation = col_annotation,               # Add the column annotation
    name='Abundance',
    show_row_names = TRUE,
    show_column_names = FALSE,                     # Hide column names (optional)
    cluster_rows = TRUE,                           # Enable row clustering
    cluster_columns = FALSE,    
    column_split = growth_group,
    # show_heatmap_legend = FALSE
    heatmap_legend_param = list(
      title_position='lefttop-rot'
    )
    )
  dev.off()
  
  # Dot whisper plot
  library(patchwork)
  
  class_palette <- c("Amino acids" = '#a6cee3',
                     "Bile acids" = '#1f78b4',
                     "Carnitines" = '#b2df8a',
                     "Fatty acids" = '#33a02c',
                     "Fatty acyls" = '#fb9a99',
                     "PA" = '#d9d9d9',
                     "PC" = '#003c30',
                     "PE" = '#ffff99',
                     "Glycerophosphoglycerols"= '#8e0152',
                     "Glycerophosphoinositols" = '#313695',
                     "Nucleotides" = '#e31a1c',
                     "Organic acids" = '#fdbf6f',
                     "Others" = '#ff7f00',
                     "Sphingolipids" = '#cab2d6',
                     "Steroids" = '#6a3d9a',
                     "Xenobiotics"  = '#b15928')
  
  plot_data <- log_result %>%
    mutate(Name = if_else(Name == 'Glycerophosphocholine (isomer)', 'Glycerophosphocholine', Name)) %>%
    mutate(Name = if_else(Name == 'Octadecanedioate', 'FA 18:1;O2', Name)) %>%
    mutate(effect = if_else(Name %in% c('Catechol Sulfate', 'LPA 22:6', 'LPA 18:2'), -1, 1)) %>%
    mutate(Class = case_when(
      Class == "Glycerophosphates" ~ "PA",
      Class == "Glycerophosphocholines" ~ "PC",
      Class == "Glycerophosphoethanolamines" ~ "PE",
      Class == "Glycerophosphoglycerols" ~ "PG",
      Class == "Glycerophosphoinositols" ~ "PI",
      TRUE ~ Class
    )) %>%
    mutate(Class=factor(Class, levels=c("Fatty acids",
                                        "Fatty acyls",
                                        "PC",
                                        "PE",
                                        "PA",
                                        "Organic acids",
                                        "Others"))) %>%
    arrange(effect, desc(Class)) %>%
    mutate(Name = factor(Name, levels = unique(Name)))
    
  dtplot <- ggplot(plot_data, aes(x = Name, y = Odds_Ratio)) +
    geom_pointrange(aes(ymin = ORCI_Lower, ymax = ORCI_Upper), size=0.4, shape=18) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    coord_flip() +
    theme_minimal() +
    labs(title = "",
         x = "", 
         y = "OR (95% CI)") +
    theme(panel.grid.major.y = element_blank(),
          axis.ticks.y = element_blank(), 
          axis.line.y = element_blank(),
          axis.text.y = element_blank(), 
          axis.title.y = element_blank(),
          strip.text = element_blank())
  
  p_y <-
    ggplot(plot_data) +
    geom_col(aes(x = Name, y = 1, fill = Class)) +
    scale_fill_manual(values = class_palette) +
    coord_cartesian(expand = F) +
    coord_flip() +
    theme(axis.text.x = element_blank(), 
          axis.text.y = element_text(color = "black", size=10), 
          axis.title.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.background = element_blank(), 
          plot.background = element_blank(),
          legend.position = "none")
  
  p_y + dtplot +
    plot_layout(widths = c(.04, 1), guides='collect') &
    theme()
  ggsave("../article/pictures/20241107_reorganize/brwe/macrosomia_sample_matching_odds_ratio.tiff",
         width = 5, height =5.7, units = c("in"), dpi = 300, bg='white')
}