library(tidyverse)
library(readxl)
library(mice)
library(patchwork)

# Perform linear regression in Batch 2
cpd_data <- read.csv('../DATA/peakmat_features/batch2_reorganize_20241107/batch2_compound_data.csv')
cpd_info <- read_excel('../DATA/statistics/identified_compounds_information_batch2.xlsx', sheet=1)
cl_data <- read.csv('../DATA/statistics/clinical_data_batch2.csv')

cl_data$VSI[is.na(cl_data$VSI)] <- mean(cl_data$VSI, na.rm=TRUE)
cl_data$WMI[is.na(cl_data$WMI)] <- mean(cl_data$WMI, na.rm=TRUE)
cl_data$PSIRevised[is.na(cl_data$PSIRevised)] <- mean(cl_data$PSIRevised, na.rm=TRUE)

birth_data <- read.csv('../DATA/statistics/clinical_birth_data_batch2.csv')

merge_data <- cpd_data %>%
  inner_join(cl_data, by='sample_name') %>%
  inner_join(birth_data, by='sample_name') %>%
  select(-fmbs_CBH46)   # 妊娠期糖尿病不是智力相关因素

cols_to_factor <- c('chsex_new', 'yqyzdc_educton', 'yqyzdc_income', 'yqyzdc_folv4_1', 'yzhdc_FOLIC1',
                    'yzhdc_SMKN3')
merge_data[cols_to_factor] <- lapply(merge_data[cols_to_factor], as.factor)

merge_data$yqyzdc_educton <- factor(merge_data$yqyzdc_educton,
                                    levels = c(2, 3, 4, 5, 6),
                                    labels = c('Low', 'Low','Low','Medium', 'High'))
merge_data$yqyzdc_income <- factor(merge_data$yqyzdc_income,
                                   levels = c(1,2,3,4,5,6,7,8),
                                   labels =c('Low','Low','Low','Low','Low',
                                             'Medium','High', 'High'))
imp_ini <- mice(merge_data, m = 1, print=F, maxit = 1, method = "pmm")
pred <- imp_ini$pred
pred
pred[, 1:500] <- 0
meth <- imp_ini$meth
meth
meth[c('chsex_new', 'yzhdc_FOLIC1', 'yzhdc_SMKN3')] <- 'logreg'
meth[c('yqyzdc_educton', 'yqyzdc_income', 'yqyzdc_folv4_1')] <- 'polr'

imp <- mice(merge_data, meth = meth, pred=pred, print=F, m=10, seed=2024)

compound_names <- names(merge_data)[2:500]
confounder_names <- names(merge_data)[c(501:504,506,507, 516)]

dependent_vars <- c("VCI", "VSI", "FRI", "WMI", "PSI", "FSIQ")
pol_sheets <- list()

for (var in dependent_vars) {
  pol_sheet <- data.frame()
  for (compound_var in compound_names) {
    if (var == "PSI") {
    formula <- paste("PSIRevised", "~", paste(c(compound_var, confounder_names), collapse = " + "))
  }
    else {
    formula <- paste(var, "~", paste(c(compound_var, confounder_names), collapse = " + "))
    }
    fit <- with(imp, lm(as.formula(formula)))
    pooled_result <- pool(fit)
    pol <- summary(pooled_result, conf.int=TRUE) %>%
      filter(term==compound_var)
    pol_sheet <- rbind(pol_sheet, pol)
  }
  pol_sheets[[var]] <- pol_sheet
}

merge_pol_sheet <- pol_sheets[['VCI']] %>%
  select(term, estimate, p.value, `2.5 %`, `97.5 %`) %>%
  rename('VCI estimate' = estimate, 'VCI p' = p.value, 'VCI lower' = `2.5 %`, 'VCI upper' = `97.5 %`) %>%
  inner_join(pol_sheets[['VSI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
             %>% rename('VSI estimate' = estimate, 'VSI p' = p.value, 'VSI lower' = `2.5 %`, 'VSI upper' = `97.5 %`), by='term') %>%
  inner_join(pol_sheets[['FRI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
             %>% rename('FRI estimate' = estimate, 'FRI p' = p.value, 'FRI lower' = `2.5 %`, 'FRI upper' = `97.5 %`), by='term') %>%
  inner_join(pol_sheets[['WMI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
             %>% rename('WMI estimate' = estimate, 'WMI p' = p.value, 'WMI lower' = `2.5 %`, 'WMI upper' = `97.5 %`), by='term') %>%
  inner_join(pol_sheets[['PSI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
             %>% rename('PSI estimate' = estimate, 'PSI p' = p.value, 'PSI lower' = `2.5 %`, 'PSI upper' = `97.5 %`), by='term') %>%
  inner_join(pol_sheets[['FSIQ']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
             %>% rename('FSIQ estimate' = estimate, 'FSIQ p' = p.value, 'FSIQ lower' = `2.5 %`, 'FSIQ upper' = `97.5 %`), by='term') %>%
  inner_join(cpd_info %>% select(unique_features, Name), by = c('term' = 'unique_features')) %>%
  mutate(`VCI FDR` = p.adjust(`VCI p`, method='fdr'),
         `VSI FDR` = p.adjust(`VSI p`, method='fdr'),
         `FRI FDR` = p.adjust(`FRI p`, method='fdr'),
         `WMI FDR` = p.adjust(`WMI p`, method='fdr'),
         `PSI FDR` = p.adjust(`PSI p`, method='fdr'),
         `FSIQ FDR` = p.adjust(`FSIQ p`, method='fdr'))
write.csv(merge_pol_sheet, './20241107_reanalysis/mlr_compounds_batch2_20250307.csv', row.names = FALSE) # Linear correlation compounds
{
  #multiple linear regression on covariates
  confounder_names <- names(merge_data)[c(503:510, 517:518)]
  dependent_vars <- c("VCI", "VSI", "FRI", "WMI", "PSI", "FSIQ")
  pol_sheets <- list()
  
  for (var in dependent_vars) {
      if (var == "PSI") {
        formula <- paste("PSIRevised", "~", paste(confounder_names, collapse = " + "))
        print(formula)
      }
      else {
        formula <- paste(var, "~", paste(confounder_names, collapse = " + "))
        print(formula)
      }
      fit <- with(imp, lm(as.formula(formula)))
      pooled_result <- pool(fit)
      pol <- summary(pooled_result, conf.int=TRUE)[-1, , drop = FALSE] %>%
        mutate(DV = var) %>%
        select(DV, term, estimate, `2.5 %`, `97.5 %`, p.value)
    pol_sheets[[var]] <- pol
  }
  
  results <- do.call(rbind, pol_sheets) %>%
    pivot_wider(
      names_from = DV,
      values_from = c(estimate, `2.5 %`, `97.5 %`, p.value)
    ) %>%
    mutate(VCI_combined = paste0(sprintf("%.3f", estimate_VCI), " (", sprintf("%.3f", `2.5 %_VCI`), ", ", sprintf("%.3f", `97.5 %_VCI`), ")")) %>%
    mutate(VSI_combined = paste0(sprintf("%.3f", estimate_VSI), " (", sprintf("%.3f", `2.5 %_VSI`), ", ", sprintf("%.3f", `97.5 %_VSI`), ")")) %>%
    mutate(FRI_combined = paste0(sprintf("%.3f", estimate_FRI), " (", sprintf("%.3f", `2.5 %_FRI`), ", ", sprintf("%.3f", `97.5 %_FRI`), ")")) %>%
    mutate(WMI_combined = paste0(sprintf("%.3f", estimate_WMI), " (", sprintf("%.3f", `2.5 %_WMI`), ", ", sprintf("%.3f", `97.5 %_WMI`), ")")) %>%
    mutate(PSI_combined = paste0(sprintf("%.3f", estimate_PSI), " (", sprintf("%.3f", `2.5 %_PSI`), ", ", sprintf("%.3f", `97.5 %_PSI`), ")")) %>%
    mutate(FSIQ_combined = paste0(sprintf("%.3f", estimate_FSIQ), " (", sprintf("%.3f", `2.5 %_FSIQ`), ", ", sprintf("%.3f", `97.5 %_FSIQ`), ")"))
  write.csv(results, "./20241107_reanalysis/iq_modeling/covariates_iq_lm_adjusted.csv",
            row.names = FALSE)
}

{
  # Different models
  # Model_1
  compound_names <- names(merge_data)[2:500]
  
  dependent_vars <- c("VCI", "VSI", "FRI", "WMI", "PSI", "FSIQ")
  pol_sheets <- list()
  
  for (var in dependent_vars) {
    pol_sheet <- data.frame()
    for (compound_var in compound_names) {
      if (var == "PSI") {
        formula <- paste("PSIRevised", "~", compound_var)
      }
      else {
        formula <- paste(var, "~", compound_var)
      }
      fit <- with(imp, lm(as.formula(formula)))
      pooled_result <- pool(fit)
      pol <- summary(pooled_result, conf.int=TRUE) %>%
        filter(term==compound_var)
      pol_sheet <- rbind(pol_sheet, pol)
    }
    pol_sheets[[var]] <- pol_sheet
  }
  
  merge_pol_sheet <- pol_sheets[['VCI']] %>%
    select(term, estimate, p.value, `2.5 %`, `97.5 %`) %>%
    rename('VCI estimate' = estimate, 'VCI p' = p.value, 'VCI lower' = `2.5 %`, 'VCI upper' = `97.5 %`) %>%
    inner_join(pol_sheets[['VSI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
               %>% rename('VSI estimate' = estimate, 'VSI p' = p.value, 'VSI lower' = `2.5 %`, 'VSI upper' = `97.5 %`), by='term') %>%
    inner_join(pol_sheets[['FRI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
               %>% rename('FRI estimate' = estimate, 'FRI p' = p.value, 'FRI lower' = `2.5 %`, 'FRI upper' = `97.5 %`), by='term') %>%
    inner_join(pol_sheets[['WMI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
               %>% rename('WMI estimate' = estimate, 'WMI p' = p.value, 'WMI lower' = `2.5 %`, 'WMI upper' = `97.5 %`), by='term') %>%
    inner_join(pol_sheets[['PSI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
               %>% rename('PSI estimate' = estimate, 'PSI p' = p.value, 'PSI lower' = `2.5 %`, 'PSI upper' = `97.5 %`), by='term') %>%
    inner_join(pol_sheets[['FSIQ']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
               %>% rename('FSIQ estimate' = estimate, 'FSIQ p' = p.value, 'FSIQ lower' = `2.5 %`, 'FSIQ upper' = `97.5 %`), by='term') %>%
    inner_join(cpd_info %>% select(unique_features, Name), by = c('term' = 'unique_features')) %>%
    mutate(`VCI FDR` = p.adjust(`VCI p`, method='fdr'),
           `VSI FDR` = p.adjust(`VSI p`, method='fdr'),
           `FRI FDR` = p.adjust(`FRI p`, method='fdr'),
           `WMI FDR` = p.adjust(`WMI p`, method='fdr'),
           `PSI FDR` = p.adjust(`PSI p`, method='fdr'),
           `FSIQ FDR` = p.adjust(`FSIQ p`, method='fdr'))
  write.csv(merge_pol_sheet, './20241107_reanalysis/mlr_compounds_batch2_no_adjustment.csv', row.names = FALSE)
  
  # Model_2 (no_adjust_BMI)
  compound_names <- names(merge_data)[2:500]
  confounder_names <- names(merge_data)[c(501:504,507, 516)]
  
  dependent_vars <- c("VCI", "VSI", "FRI", "WMI", "PSI", "FSIQ")
  pol_sheets <- list()
  
  for (var in dependent_vars) {
    pol_sheet <- data.frame()
    for (compound_var in compound_names) {
      if (var == "PSI") {
        formula <- paste("PSIRevised", "~", paste(c(compound_var, confounder_names), collapse = " + "))
      }
      else {
        formula <- paste(var, "~", paste(c(compound_var, confounder_names), collapse = " + "))
      }
      fit <- with(imp, lm(as.formula(formula)))
      pooled_result <- pool(fit)
      pol <- summary(pooled_result, conf.int=TRUE) %>%
        filter(term==compound_var)
      pol_sheet <- rbind(pol_sheet, pol)
    }
    pol_sheets[[var]] <- pol_sheet
  }
  
  merge_pol_sheet <- pol_sheets[['VCI']] %>%
    select(term, estimate, p.value, `2.5 %`, `97.5 %`) %>%
    rename('VCI estimate' = estimate, 'VCI p' = p.value, 'VCI lower' = `2.5 %`, 'VCI upper' = `97.5 %`) %>%
    inner_join(pol_sheets[['VSI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
               %>% rename('VSI estimate' = estimate, 'VSI p' = p.value, 'VSI lower' = `2.5 %`, 'VSI upper' = `97.5 %`), by='term') %>%
    inner_join(pol_sheets[['FRI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
               %>% rename('FRI estimate' = estimate, 'FRI p' = p.value, 'FRI lower' = `2.5 %`, 'FRI upper' = `97.5 %`), by='term') %>%
    inner_join(pol_sheets[['WMI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
               %>% rename('WMI estimate' = estimate, 'WMI p' = p.value, 'WMI lower' = `2.5 %`, 'WMI upper' = `97.5 %`), by='term') %>%
    inner_join(pol_sheets[['PSI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
               %>% rename('PSI estimate' = estimate, 'PSI p' = p.value, 'PSI lower' = `2.5 %`, 'PSI upper' = `97.5 %`), by='term') %>%
    inner_join(pol_sheets[['FSIQ']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
               %>% rename('FSIQ estimate' = estimate, 'FSIQ p' = p.value, 'FSIQ lower' = `2.5 %`, 'FSIQ upper' = `97.5 %`), by='term') %>%
    inner_join(cpd_info %>% select(unique_features, Name), by = c('term' = 'unique_features')) %>%
    mutate(`VCI FDR` = p.adjust(`VCI p`, method='fdr'),
           `VSI FDR` = p.adjust(`VSI p`, method='fdr'),
           `FRI FDR` = p.adjust(`FRI p`, method='fdr'),
           `WMI FDR` = p.adjust(`WMI p`, method='fdr'),
           `PSI FDR` = p.adjust(`PSI p`, method='fdr'),
           `FSIQ FDR` = p.adjust(`FSIQ p`, method='fdr'))
  write.csv(merge_pol_sheet, './20241107_reanalysis/mlr_compounds_batch2_no_adjust_bmi.csv', row.names = FALSE)
  # Model_3 (no_adjust_income)
  # compound_names <- names(merge_data)[2:500]
  # confounder_names <- names(merge_data)[c(501:502,504, 506,507, 516)]
  # 
  # dependent_vars <- c("VCI", "VSI", "FRI", "WMI", "PSI", "FSIQ")
  # pol_sheets <- list()
  # 
  # for (var in dependent_vars) {
  #   pol_sheet <- data.frame()
  #   for (compound_var in compound_names) {
  #     if (var == "PSI") {
  #       formula <- paste("PSIRevised", "~", paste(c(compound_var, confounder_names), collapse = " + "))
  #     }
  #     else {
  #       formula <- paste(var, "~", paste(c(compound_var, confounder_names), collapse = " + "))
  #     }
  #     fit <- with(imp, lm(as.formula(formula)))
  #     pooled_result <- pool(fit)
  #     pol <- summary(pooled_result, conf.int=TRUE) %>%
  #       filter(term==compound_var)
  #     pol_sheet <- rbind(pol_sheet, pol)
  #   }
  #   pol_sheets[[var]] <- pol_sheet
  # }
  # 
  # merge_pol_sheet <- pol_sheets[['VCI']] %>%
  #   select(term, estimate, p.value, `2.5 %`, `97.5 %`) %>%
  #   rename('VCI estimate' = estimate, 'VCI p' = p.value, 'VCI lower' = `2.5 %`, 'VCI upper' = `97.5 %`) %>%
  #   inner_join(pol_sheets[['VSI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
  #              %>% rename('VSI estimate' = estimate, 'VSI p' = p.value, 'VSI lower' = `2.5 %`, 'VSI upper' = `97.5 %`), by='term') %>%
  #   inner_join(pol_sheets[['FRI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
  #              %>% rename('FRI estimate' = estimate, 'FRI p' = p.value, 'FRI lower' = `2.5 %`, 'FRI upper' = `97.5 %`), by='term') %>%
  #   inner_join(pol_sheets[['WMI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
  #              %>% rename('WMI estimate' = estimate, 'WMI p' = p.value, 'WMI lower' = `2.5 %`, 'WMI upper' = `97.5 %`), by='term') %>%
  #   inner_join(pol_sheets[['PSI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
  #              %>% rename('PSI estimate' = estimate, 'PSI p' = p.value, 'PSI lower' = `2.5 %`, 'PSI upper' = `97.5 %`), by='term') %>%
  #   inner_join(pol_sheets[['FSIQ']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
  #              %>% rename('FSIQ estimate' = estimate, 'FSIQ p' = p.value, 'FSIQ lower' = `2.5 %`, 'FSIQ upper' = `97.5 %`), by='term') %>%
  #   inner_join(cpd_info %>% select(unique_features, Name), by = c('term' = 'unique_features')) %>%
  #   mutate(`VCI FDR` = p.adjust(`VCI p`, method='fdr'),
  #          `VSI FDR` = p.adjust(`VSI p`, method='fdr'),
  #          `FRI FDR` = p.adjust(`FRI p`, method='fdr'),
  #          `WMI FDR` = p.adjust(`WMI p`, method='fdr'),
  #          `PSI FDR` = p.adjust(`PSI p`, method='fdr'),
  #          `FSIQ FDR` = p.adjust(`FSIQ p`, method='fdr'))
  # write.csv(merge_pol_sheet, './20241107_reanalysis/mlr_compounds_batch2_no_adjust_income.csv', row.names = FALSE)
  
  # Model_4 (no_adjust_income&education)
  compound_names <- names(merge_data)[2:500]
  confounder_names <- names(merge_data)[c(501,504, 506,507, 516)]
  
  dependent_vars <- c("VCI", "VSI", "FRI", "WMI", "PSI", "FSIQ")
  pol_sheets <- list()
  
  for (var in dependent_vars) {
    pol_sheet <- data.frame()
    for (compound_var in compound_names) {
      if (var == "PSI") {
        formula <- paste("PSIRevised", "~", paste(c(compound_var, confounder_names), collapse = " + "))
      }
      else {
        formula <- paste(var, "~", paste(c(compound_var, confounder_names), collapse = " + "))
      }
      fit <- with(imp, lm(as.formula(formula)))
      pooled_result <- pool(fit)
      pol <- summary(pooled_result, conf.int=TRUE) %>%
        filter(term==compound_var)
      pol_sheet <- rbind(pol_sheet, pol)
    }
    pol_sheets[[var]] <- pol_sheet
  }
  
  merge_pol_sheet <- pol_sheets[['VCI']] %>%
    select(term, estimate, p.value, `2.5 %`, `97.5 %`) %>%
    rename('VCI estimate' = estimate, 'VCI p' = p.value, 'VCI lower' = `2.5 %`, 'VCI upper' = `97.5 %`) %>%
    inner_join(pol_sheets[['VSI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
               %>% rename('VSI estimate' = estimate, 'VSI p' = p.value, 'VSI lower' = `2.5 %`, 'VSI upper' = `97.5 %`), by='term') %>%
    inner_join(pol_sheets[['FRI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
               %>% rename('FRI estimate' = estimate, 'FRI p' = p.value, 'FRI lower' = `2.5 %`, 'FRI upper' = `97.5 %`), by='term') %>%
    inner_join(pol_sheets[['WMI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
               %>% rename('WMI estimate' = estimate, 'WMI p' = p.value, 'WMI lower' = `2.5 %`, 'WMI upper' = `97.5 %`), by='term') %>%
    inner_join(pol_sheets[['PSI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
               %>% rename('PSI estimate' = estimate, 'PSI p' = p.value, 'PSI lower' = `2.5 %`, 'PSI upper' = `97.5 %`), by='term') %>%
    inner_join(pol_sheets[['FSIQ']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
               %>% rename('FSIQ estimate' = estimate, 'FSIQ p' = p.value, 'FSIQ lower' = `2.5 %`, 'FSIQ upper' = `97.5 %`), by='term') %>%
    inner_join(cpd_info %>% select(unique_features, Name), by = c('term' = 'unique_features')) %>%
    mutate(`VCI FDR` = p.adjust(`VCI p`, method='fdr'),
           `VSI FDR` = p.adjust(`VSI p`, method='fdr'),
           `FRI FDR` = p.adjust(`FRI p`, method='fdr'),
           `WMI FDR` = p.adjust(`WMI p`, method='fdr'),
           `PSI FDR` = p.adjust(`PSI p`, method='fdr'),
           `FSIQ FDR` = p.adjust(`FSIQ p`, method='fdr'))
  write.csv(merge_pol_sheet, './20241107_reanalysis/mlr_compounds_batch2_no_adjust_income_education.csv', row.names = FALSE)
  
  
  # Plot of different models (VCI)
  cpd_info <- read_excel('../DATA/statistics/identified_compounds_information_batch2.xlsx', sheet=1)
  
  model_3_data <- read.csv('./20241107_reanalysis/mlr_compounds_batch2_20241118.csv', check.names = FALSE) %>%
    # filter(p.value < 0.05) %>%
    mutate(Model = 'Model_3')
  model_2_data <- read.csv('./20241107_reanalysis/mlr_compounds_batch2_no_adjust_bmi.csv', check.names = FALSE) %>%
    # filter(term %in% model_3_data$term) %>%
    mutate(Model = 'Model_2')
  model_1_data <- read.csv('./20241107_reanalysis/mlr_compounds_batch2_no_adjustment.csv', check.names = FALSE) %>%
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
    mutate(Significance = ifelse(`VCI p` < 0.05, "Significant", "Not Significant")) %>%
    arrange(Class, `VCI estimate`) %>%
    mutate(Name = factor(Name, levels = unique(Name))) %>%
    mutate(Class = if_else(Class == 'Glycerophosphates', 'PA', Class),
           Class = if_else(Class == 'Glycerophosphocholines', 'PC', Class),
           Class = if_else(Class == 'Glycerophosphoethanolamines', 'PE', Class),
           Class = if_else(Class == 'Glycerophosphoglycerols', 'PG', Class),
           Class = if_else(Class == 'Glycerophosphoinositols', 'PI', Class))
  
  plot_function <- function(class_list) {
    plot_data_c <- plot_data %>%
      filter(Class %in% class_list)
    dtplot <- ggplot(plot_data_c, aes(x = Name, y = `VCI estimate`, shape = Model, color = Significance)) +
      geom_pointrange(aes(ymin = `VCI lower`, ymax = `VCI upper`), 
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
  
  plot_function(c("Carnitines",'Fatty acyls'))
  ggsave("../article/pictures/20241107_reorganize/mlr/dot_whisker_plot_three_models_car.png", 
         width = 15, height =6.6, units = c("in"), dpi = 300, bg='white')
  plot_function(c('PC'))
  ggsave("../article/pictures/20241107_reorganize/mlr/dot_whisker_plot_three_models_pc.png", 
         width = 15, height =6.6, units = c("in"), dpi = 300, bg='white')
  plot_function(c('PE', 'PI', 'PG', 'PA'))
  ggsave("../article/pictures/20241107_reorganize/mlr/dot_whisker_plot_three_models_pe.png", 
         width = 15, height =6.6, units = c("in"), dpi = 300, bg='white')
  plot_function(c("Amino acids","Organic acids"))
  ggsave("../article/pictures/20241107_reorganize/mlr/dot_whisker_plot_three_models_aa.png", 
         width = 15, height =6.6, units = c("in"), dpi = 300, bg='white')
  plot_function(c("Bile acids", "Steroids"))
  ggsave("../article/pictures/20241107_reorganize/mlr/dot_whisker_plot_three_models_ba.png", 
         width = 15, height =6.6, units = c("in"), dpi = 300, bg='white')
  plot_function(c("Nucleotides","Sphingolipids", "Xenobiotics"))
  ggsave("../article/pictures/20241107_reorganize/mlr/dot_whisker_plot_three_models_nucle.png", 
         width = 15, height =6.6, units = c("in"), dpi = 300, bg='white')
  plot_function(c("Others"))
  ggsave("../article/pictures/20241107_reorganize/mlr/dot_whisker_plot_three_models_other.png", 
         width = 15, height =6.6, units = c("in"), dpi = 300, bg='white')
  
  # Plot for FA displaying trend of C atoms
  plot_data_c <- plot_data %>%
    filter(Class == 'Fatty acids') %>%
    inner_join(read_excel('../DATA/statistics/FA_structures.xlsx'), by = c('Name' = 'FA_name')) %>%
    filter(sum_comp != '_' & keep==1) %>%
    mutate(sum_comp = factor(sum_comp, levels = unique(sum_comp)))
  dtplot <- ggplot(plot_data_c, aes(x = sum_comp, y = `VCI estimate`, shape = Model, color = Significance)) +
    geom_pointrange(aes(ymin = `VCI lower`, ymax = `VCI upper`), 
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
          axis.text.x = element_blank(),
          # axis.text.x = element_text(color = "black", size=8, angle=90,hjust=1.2, vjust=0.8),
          legend.position = c(0.5, 0.9),
          legend.direction = 'horizontal',
          axis.title.x = element_blank(),
          strip.text = element_blank())
  
  p_y <-
    ggplot(plot_data_c) +
    geom_point(aes(x = sum_comp, y = DB_number, color = "#Double bonds")) +
    geom_line(aes(x = sum_comp, y = DB_number, group = 1, color = "#Double bonds"), linewidth=0.6) +
    geom_point(aes(x = sum_comp, y = C_number, color = "#Carbon atoms")) +
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
  
  ggsave("../article/pictures/20241107_reorganize/mlr/dot_whisker_plot_three_models_FA.png", 
         width = 16, height =6.6, units = c("in"), dpi = 300, bg='white')
  
  dtplot
  ggsave("../article/pictures/20241107_reorganize/mlr/dot_whisker_plot_three_models_no_number.png", 
         width = 10, height =4, units = c("in"), dpi = 300, bg='white')
}

{
  # Regression after removing some samples
  lower_bound <- quantile(merge_data$FSIQ, 0.05, na.rm = TRUE)
  upper_bound <- quantile(merge_data$FSIQ, 0.95, na.rm = TRUE)
  
  rm_merge_data <- merge_data[merge_data$FSIQ >= lower_bound & merge_data$FSIQ <= upper_bound, ]
  imp_ini <- mice(rm_merge_data, m = 1, print=F, maxit = 1, method = "pmm")
  pred <- imp_ini$pred
  pred
  pred[, 1:500] <- 0
  meth <- imp_ini$meth
  meth
  meth[c('chsex_new', 'yzhdc_FOLIC1', 'yzhdc_SMKN3')] <- 'logreg'
  meth[c('yqyzdc_educton', 'yqyzdc_income', 'yqyzdc_folv4_1')] <- 'polr'
  
  imp <- mice(rm_merge_data, meth = meth, pred=pred, print=F, m=10, seed=2024)
  
  compound_names <- names(rm_merge_data)[2:500]
  confounder_names <- names(rm_merge_data)[c(501:504,506,507, 516)]
  
  dependent_vars <- c("VCI", "VSI", "FRI", "WMI", "PSI", "FSIQ")
  pol_sheets <- list()
  
  for (var in dependent_vars) {
    pol_sheet <- data.frame()
    for (compound_var in compound_names) {
      if (var == "PSI") {
        formula <- paste("PSIRevised", "~", paste(c(compound_var, confounder_names), collapse = " + "))
      }
      else {
        formula <- paste(var, "~", paste(c(compound_var, confounder_names), collapse = " + "))
      }
      fit <- with(imp, lm(as.formula(formula)))
      pooled_result <- pool(fit)
      pol <- summary(pooled_result, conf.int=TRUE) %>%
        filter(term==compound_var)
      pol_sheet <- rbind(pol_sheet, pol)
    }
    pol_sheets[[var]] <- pol_sheet
  }
  
  merge_pol_sheet <- pol_sheets[['VCI']] %>%
    select(term, estimate, p.value, `2.5 %`, `97.5 %`) %>%
    rename('VCI estimate' = estimate, 'VCI p' = p.value, 'VCI lower' = `2.5 %`, 'VCI upper' = `97.5 %`) %>%
    inner_join(pol_sheets[['VSI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
               %>% rename('VSI estimate' = estimate, 'VSI p' = p.value, 'VSI lower' = `2.5 %`, 'VSI upper' = `97.5 %`), by='term') %>%
    inner_join(pol_sheets[['FRI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
               %>% rename('FRI estimate' = estimate, 'FRI p' = p.value, 'FRI lower' = `2.5 %`, 'FRI upper' = `97.5 %`), by='term') %>%
    inner_join(pol_sheets[['WMI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
               %>% rename('WMI estimate' = estimate, 'WMI p' = p.value, 'WMI lower' = `2.5 %`, 'WMI upper' = `97.5 %`), by='term') %>%
    inner_join(pol_sheets[['PSI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
               %>% rename('PSI estimate' = estimate, 'PSI p' = p.value, 'PSI lower' = `2.5 %`, 'PSI upper' = `97.5 %`), by='term') %>%
    inner_join(pol_sheets[['FSIQ']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`) 
               %>% rename('FSIQ estimate' = estimate, 'FSIQ p' = p.value, 'FSIQ lower' = `2.5 %`, 'FSIQ upper' = `97.5 %`), by='term') %>%
    inner_join(cpd_info %>% select(unique_features, Name), by = c('term' = 'unique_features')) %>%
    mutate(`VCI FDR` = p.adjust(`VCI p`, method='fdr'),
           `VSI FDR` = p.adjust(`VSI p`, method='fdr'),
           `FRI FDR` = p.adjust(`FRI p`, method='fdr'),
           `WMI FDR` = p.adjust(`WMI p`, method='fdr'),
           `PSI FDR` = p.adjust(`PSI p`, method='fdr'),
           `FSIQ FDR` = p.adjust(`FSIQ p`, method='fdr'))
  write.csv(merge_pol_sheet, './20241107_reanalysis/mlr_compounds_batch2_remove_extreme.csv', row.names = FALSE)
}

{
  # Plot of significant compounds in different models (Sensitivity analysis)
  model_1_data <- read.csv('./20241107_reanalysis/mlr_compounds_batch2_20250307.csv', check.names = FALSE) %>%
    mutate(Model = 'Model_1')
  model_2_data <- read.csv('./20241107_reanalysis/mlr_compounds_batch2_no_adjust_bmi.csv', check.names = FALSE) %>%
    mutate(Model = 'Model_2')
  model_3_data <- read.csv('./20241107_reanalysis/mlr_compounds_batch2_no_adjust_income_education.csv', check.names = FALSE) %>%
    mutate(Model = 'Model_3')
  model_4_data <- read.csv('./20241107_reanalysis/mlr_compounds_batch2_no_adjustment.csv', check.names = FALSE) %>%
    mutate(Model = 'Model_4')
  model_5_data <- read.csv('./20241107_reanalysis/mlr_compounds_batch2_remove_extreme.csv', check.names = FALSE) %>%
    mutate(Model = 'Model_5')
  
  id_var <- c('VCI', 'VSI', 'FRI', 'WMI', 'PSI', 'FSIQ')
  for (var in id_var) {
    var_p <- paste(var, 'p')
    var_estimate <- paste(var, 'estimate')
    var_lower <- paste(var, 'lower')
    var_upper <- paste(var, 'upper')
    sig_cpds <- model_1_data[model_1_data[[var_p]] < 0.05, ]$term
    
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
      mutate(Significance = ifelse(!!sym(var_p) < 0.05, "Significant", "Not Significant")) %>%
      arrange(!!sym(var_estimate), Class) %>%
      mutate(Name = factor(Name, levels = unique(Name)))
    
    ggplot(plot_data, aes(x = Name, y = !!sym(var_estimate), shape = Model, color = Significance)) +
      geom_pointrange(aes(ymin =!!sym(var_lower), ymax = !!sym(var_upper)), 
                      position = position_dodge(width = 0.4)) +  # Separate dependent variables to avoid overlap
      geom_point(size = 0.01, position = position_dodge(width = 0.4)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      scale_color_manual(values = c("Not Significant" = "black", "Significant" = "#7570b3")) +
      scale_shape_manual(values = c(16, 17, 18, 15, 14)) + 
      scale_x_discrete(expand = c(0.1, 0.1)) +
      # coord_flip() +
      theme_minimal() +
      labs(title = "",
           y = "Beta Estimate (95% CI)", 
           x = "") +
      theme(panel.grid.major.x = element_blank(),
            axis.ticks.x = element_blank(), 
            axis.line.x = element_blank(),
            axis.text.x = element_text(color = "black", size=9, angle=45, hjust=1, vjust=1),
            axis.title.y = element_text(color = "black", size=9),
            axis.title.x = element_blank(),
            strip.text = element_blank(),
            legend.position = 'bottom')
    ggsave(paste0("../article/pictures/20241107_reorganize/mlr/sentivity_", var, '.png'),
           width = 15, height=8, units = 'in', dpi = 300, bg='white')
    
  }
}

# Dot-whisper plot for linear regression
merge_pol_sheet <- read.csv("./20241107_reanalysis/mlr_compounds_batch2_20250307.csv", check.names = FALSE)
# 使用FDR计算没有小于0.05的代谢物
# 两种方法展示： 1. bar plot展示智力相关代谢物的分布； 2. 展示验证化合物的dot-whisper plot

{
  # Bar plot to display the related classes
  merge_pol_sheet <- read.csv("./20241107_reanalysis/mlr_compounds_batch2_20250307.csv", check.names = FALSE) %>%
    inner_join(cpd_info %>% select(unique_features, Class_pie), by = c('term' = 'unique_features'))
  
  iq_data <- merge_pol_sheet %>%
    select(term, Class_pie, `VCI p`, `VSI p`, `FRI p`, `WMI p`, `PSI p`, `FSIQ p`) %>%
    pivot_longer(cols = ends_with("p"), 
                 names_to = "IQ", 
                 values_to = "p_value") %>%
    mutate(IQ = sub("\\ p", "", IQ)) # Remove ".p" suffix from IQ names
  
  # Filter significant features (p_value < 0.05)
  significant_features <- iq_data %>%
    filter(p_value < 0.05)
  
  # Count significant features by IQ and Class
  significant_counts <- significant_features %>%
    group_by(IQ, Class_pie) %>%
    summarize(count = n(), .groups = "drop") %>%
    mutate(IQ = factor(IQ, levels=c('VCI', 'VSI', 'FRI', 'WMI', 'PSI', 'FSIQ')))
  
  # Calculate total counts for each IQ to place annotations
  total_counts <- significant_counts %>%
    group_by(IQ) %>%
    summarize(total = sum(count), .groups = "drop")
  
  # Create the bar plot
  ggplot(significant_counts, aes(x = IQ, y = count, fill = Class_pie)) +
    geom_bar(stat = "identity", position = "stack", width=0.7) +
    scale_fill_manual(
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
    geom_text(data = total_counts, aes(x = IQ, y = total + 1.3, label = total), inherit.aes = FALSE) +
    labs(x = '',
         y = "Number of Significant Compounds",
         fill = "Class") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black"))
  ggsave("../article/pictures/20241107_reorganize/mlr/batch2_sig_cpd_classes.png",
         width = 5,
         height =4, 
         units = c("in"), 
         dpi = 300, bg='white')
}

{
  # Upset plot
  merge_pol_sheet <- read.csv("./20241107_reanalysis/mlr_compounds_batch2_20241118.csv", check.names = FALSE) %>%
    inner_join(cpd_info %>% select(unique_features, Class_pie), by = c('term' = 'unique_features'))
  
  {
  # Upset plot of different iq indexes (313电脑)
  cpd_info <- read_excel('./data/identified_compounds_information_batch2.xlsx', sheet=1)
  merge_pol_sheet <- read.csv("./data/mlr_compounds_batch2_20241118.csv", check.names = FALSE) %>%
    inner_join(cpd_info %>% select(unique_features, Class_pie), by = c('term' = 'unique_features')) %>%
    mutate(VCI = if_else(`VCI p` < 0.05, TRUE, FALSE),
           VSI = if_else(`VSI p` < 0.05, TRUE, FALSE),
           FRI = if_else(`FRI p` < 0.05, TRUE, FALSE),
           WMI = if_else(`WMI p` < 0.05, TRUE, FALSE),
           PSI = if_else(`PSI p` < 0.05, TRUE, FALSE),
           FSIQ = if_else(`FSIQ p` < 0.05, TRUE, FALSE)) %>%
    filter(VCI == TRUE | VSI == TRUE |FRI == TRUE|WMI == TRUE |PSI ==TRUE|FSIQ==TRUE) %>%
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
  
  cols <-  colnames(merge_pol_sheet)[28:33]
  
  upset(merge_pol_sheet, 
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
        min_size=5,
        sort_sets=FALSE,
        sort_intersections='descending',
        name='Significant compounds (p < 0.05) in different scores'
  )
  ggsave('./upset_IQ_indexes_comparison.png', width = 6, height = 5, units = 'in',
         dpi=300)
  
  
}
}

{
  # Dot plot of significant correlated compounds
  batch2_sig <- read.csv('./20241107_reanalysis/mlr_compounds_batch2_20250307.csv', check.names = FALSE)
  combined_df <- read.csv('./20241107_reanalysis/correlation_metab_iq_batch2.csv')
  sig_cpd <- combined_df %>%
    filter(VCI_ap < 0.05 | VSI_ap < 0.05 | PSI_ap < 0.05 |FRI_ap < 0.05|WMI_ap < 0.05|FSIQ_ap < 0.05) 
  
  cpd_info <- read_excel('../DATA/statistics/identified_compounds_information_batch2.xlsx', sheet=1)
  class_palette <- c("Amino acids" = '#a6cee3',
                     "Bile acids" = '#1f78b4',
                     "Carnitines" = '#b2df8a',
                     "Fatty acids" = '#33a02c',
                     "Fatty acyls" = '#fb9a99',
                     "Glycerophosphates" = '#d9d9d9',
                     "Glycerophosphocholines" = '#003c30',
                     "Glycerophosphoethanolamines" = '#ffff99',
                     "Glycerophosphoglycerols"= '#8e0152',
                     "Glycerophosphoinositols" = '#313695',
                     "Nucleotides" = '#e31a1c',
                     "Organic acids" = '#fdbf6f',
                     "Others" = '#ff7f00',
                     "Sphingolipids" = '#cab2d6',
                     "Steroids" = '#6a3d9a',
                     "Xenobiotics"  = '#b15928')
  filtered_data <- batch2_sig %>%
    filter(term %in% sig_cpd$unique_features) %>%
    mutate(effect = if_else(`FSIQ estimate` < 0, -1, 1))
  
  plot_data <- filtered_data %>%
    pivot_longer(
      cols = starts_with("VCI") | starts_with("VSI") | starts_with("FRI") | starts_with("WMI") | starts_with("PSI") | starts_with("FSIQ"),
      names_to = c("index", ".value"), 
      names_pattern = "(.*) (.*)"
    ) %>%
    mutate(significance = ifelse(p < 0.05, "Significant", "Not Significant")) %>%
    rename('Index' = 'index', 'Significance' = 'significance') %>%
    inner_join(cpd_info %>% select(Name, Class), by=c("Name" = "Name")) %>%
    mutate(Name = ifelse(Name == 'Perfluorononanoic acid', 'PFNA', Name)) %>%
    mutate(Name = ifelse(Name == 'Perfluorodecanoic acid', 'PFDA', Name)) %>%
    mutate(Name = ifelse(Name == 'Perfluoroundecanoic acid', 'PFUnDA', Name)) %>%   # Perfluorooctanesulfonic acid, PFOS
    mutate(Name = ifelse(Name == 'Linoleic acid', 'FA 18:2', Name)) %>%
    group_by(Name) %>%  # Group by the compound
    mutate(y = ifelse(row_number() == 1, 1, 0)) %>% # Add 'y', 1 for the first row, 0 for others
    ungroup() %>%
    arrange(effect, Class) %>%  # First sort by 'effect' and then by 'Class'
    mutate(Name = factor(Name, levels = unique(Name)))  %>%
    mutate(Index = factor(Index, levels = c("VCI", "VSI", "FRI", "WMI", "PSI", "FSIQ")))
  
  dtplot <- ggplot(plot_data, aes(x = Name, y = estimate, shape = Index, color = Significance)) +
    geom_pointrange(aes(ymin = lower, ymax = upper), 
                    position = position_dodge(width = 0.6)) +  # Separate dependent variables to avoid overlap
    geom_point(size = 0.01, position = position_dodge(width = 0.6)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("Significant" = "#8856a7", "Not Significant" = "black")) +  # Red for p < 0.05
    scale_shape_manual(values = c(21, 22, 23, 24, 25, 8)) +  # Assign different shapes for different dependent variables
    coord_flip() +
    theme_minimal() +
    labs(title = "",
         x = "", 
         y = "Beta Estimate (95% CI)") +
    theme(panel.grid.major.y = element_blank(),
          axis.ticks.y = element_blank(), 
          axis.line.y = element_blank(),
          axis.text.y = element_blank(), 
          axis.title.y = element_blank())
  
  p_y <-
    ggplot(plot_data) +
    geom_col(aes(x = Name, y = y, fill = Class)) +
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
          legend.position = 'none')
  
  p_y + dtplot +
    plot_layout(widths = c(.07, 1), guides='collect') &
    theme()
  ggsave("../article/pictures/20241107_reorganize/mlr/dot_whisker_plot_batch2_sig_correlated_compounds.png", width = 10, height =11, units = c("in"), dpi = 300, bg='white')
}

{
  # 画出显著代谢物的趋势分布
  d_vars <- c("VCI", "VSI", "FRI", "WMI", "PSI", "FSIQ")
  for (d_var in d_vars) {
    d_var_p <- paste(d_var, 'p')
    d_var_est <- paste(d_var, 'estimate')
    batch2_sig <- read.csv('./20241107_reanalysis/mlr_compounds_batch2_20250307.csv', check.names = FALSE) %>%
      filter(!!sym(d_var_p) < 0.05) %>%
      inner_join(cpd_info %>% select(Name, Class), by=c("Name" = "Name")) %>%
      mutate(Effect = if_else(!!sym(d_var_est) > 0, 'Positive', 'Negative')) %>%
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
                                        "Sphingolipids",
                                        "Amino acids",
                                        "Steroids",
                                        "Carnitines",
                                        "Bile acids",
                                        "Organic acids",
                                        "Nucleotides",
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
  ggsave(paste0("../article/pictures/20241107_reorganize/mlr/batch2_sig_compounds_", d_var,".png"), 
         width = 4, height = 4, units = c("in"), dpi = 300, bg='white')
  }
  
}

# Visualize the linear trend
colnames(merge_data)[colnames(merge_data) == "PSIRevised"] <- "PSI"
for (i in 1:nrow(validated)) {
  x <- validated[i, 'cpd']
  x_name <- cpd_info[cpd_info['unique_features'] == x, 'Name']
  y <- validated[i, 'var']
  
  # # Perform the linear regression to extract statistics
  # lm_model <- lm(as.formula(paste(y, "~", x)), data = merge_data)
  # lm_summary <- summary(lm_model)
  # slope <- round(lm_summary$coefficients[2, 1], 3) # Coefficient for x
  # intercept <- round(lm_summary$coefficients[1, 1], 3) # Intercept
  # r_squared <- round(lm_summary$r.squared, 3) # R-squared
  
  plot <- ggplot(merge_data, aes_string(x = x, y = y)) +
    geom_point(size = 0.5) +
    geom_smooth(method = 'lm', se = TRUE) + # Add regression line with CI
    labs(x = x_name, y = y) +
    geom_rug(linewidth=0.2) + # Add rug plot
    # annotate("text", x = Inf, y = Inf, label = paste0("y = ", intercept, " + ", slope, "x\nR² = ", r_squared),
    #          hjust = 1.1, vjust = 2, size = 5, color = "blue", parse = FALSE) + # Add regression info
    theme_minimal()
  
  # Save the plot
  ggsave(filename = paste0('../article/pictures/20241107_reorganize/mlr/lm_', x, '_', y, '.png'),
         plot = plot, width = 4, height = 4, units = "in", dpi = 300, bg = 'white')
}




