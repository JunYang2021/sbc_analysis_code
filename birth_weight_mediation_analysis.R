library(tidyverse)
library(readxl)
library(mediation)

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
}

{
  # mediation: metabolites -> gestational age -> birth weight
  merge_data <- merge_data %>%
    filter(!is.na(yunzhou))
  metabolites <- colnames(merge_data)[2:505]
  y <- 'brwe'
  mediator <- 'yunzhou'
  df <- data.frame(
    X_var = character(),
    ACME_Estimate = numeric(),
    ACME_Lower = numeric(),
    ACME_Upper = numeric(),
    ACME_p_value = numeric(),
    ADE_Estimate = numeric(),
    ADE_Lower = numeric(),
    ADE_Upper = numeric(),
    ADE_p_value = numeric(),
    Total_Effect = numeric(),
    TE_Lower = numeric(),
    TE_Upper = numeric(),
    TE_p_value = numeric(),
    Prop_Mediated = numeric(),
    PM_Lower = numeric(),
    PM_Upper = numeric(),
    PM_p_value = numeric(),
    stringsAsFactors = FALSE
  )
  for (x in metabolites) {
    form_mediator <- paste(mediator, '~', x)
    form_outcome <- paste(y, '~', x, '+', mediator)
    mediator_model <- lm(as.formula(form_mediator), data= merge_data)
    outcome_model <- lm(as.formula(form_outcome), data = merge_data)
    mediation_result <- mediate(mediator_model, outcome_model, treat = x, mediator = mediator, sims = 100)
    
    acme_est <- mediation_result$d0
    acme_lower <- mediation_result$d0.ci[1]
    acme_upper <- mediation_result$d0.ci[2]
    acme_p <- mediation_result$d0.p
    
    ade_est <- mediation_result$z0
    ade_lower <- mediation_result$z0.ci[1]
    ade_upper <- mediation_result$z0.ci[2]
    ade_p <- mediation_result$z0.p
    
    total_effect <- mediation_result$tau.coef
    total_lower <- mediation_result$tau.ci[1]
    total_upper <- mediation_result$tau.ci[2]
    total_p <- mediation_result$tau.p
    
    prop_mediated <- mediation_result$n0
    prop_lower <- mediation_result$n0.ci[1]
    prop_upper <- mediation_result$n0.ci[2]
    prop_p <- mediation_result$n0.p
    
    df <- rbind(df, data.frame(
      X_var = x,
      ACME_Estimate = acme_est,
      ACME_Lower = acme_lower,
      ACME_Upper = acme_upper,
      ACME_p_value = acme_p,
      ADE_Estimate = ade_est,
      ADE_Lower = ade_lower,
      ADE_Upper = ade_upper,
      ADE_p_value = ade_p,
      Total_Effect = total_effect,
      TE_Lower = total_lower,
      TE_Upper = total_upper,
      TE_p_value = total_p,
      Prop_Mediated = prop_mediated,
      PM_Lower = prop_lower,
      PM_Upper = prop_upper,
      PM_p_value = prop_p
    ))
  }
  df <- df %>%
    inner_join(cpd_info %>% dplyr::select(Name, unique_features), by=c("X_var" = "unique_features"))
  write.csv(df, '20241107_reanalysis/brwe_modeling/metabolites_ga_bw_mediation.csv', row.names = FALSE)
  filtered_df <- df %>%
    filter(ACME_p_value < 0.05)
}

{
  # mediation: pre-pregnancy BMI -> metabolites -> birth weight
  merge_data <- merge_data %>%
    filter(!is.na(before_bmi)) # 1103
  metabolites <- colnames(merge_data)[2:500]
  y <- 'brwe'
  x <- 'before_bmi'
  df <- data.frame(
    Mediator = character(),
    ACME_Estimate = numeric(),
    ACME_Lower = numeric(),
    ACME_Upper = numeric(),
    ACME_p_value = numeric(),
    ADE_Estimate = numeric(),
    ADE_Lower = numeric(),
    ADE_Upper = numeric(),
    ADE_p_value = numeric(),
    Total_Effect = numeric(),
    TE_Lower = numeric(),
    TE_Upper = numeric(),
    TE_p_value = numeric(),
    Prop_Mediated = numeric(),
    PM_Lower = numeric(),
    PM_Upper = numeric(),
    PM_p_value = numeric(),
    stringsAsFactors = FALSE
  )
  for (mediator in metabolites) {
    form_mediator <- paste(mediator, '~', x)
    form_outcome <- paste(y, '~', x, '+', mediator)
    mediator_model <- lm(as.formula(form_mediator), data= merge_data)
    outcome_model <- lm(as.formula(form_outcome), data = merge_data)
    mediation_result <- mediate(mediator_model, outcome_model, treat = x, mediator = mediator, sims = 1000)
    
    acme_est <- mediation_result$d0
    acme_lower <- mediation_result$d0.ci[1]
    acme_upper <- mediation_result$d0.ci[2]
    acme_p <- mediation_result$d0.p
    
    ade_est <- mediation_result$z0
    ade_lower <- mediation_result$z0.ci[1]
    ade_upper <- mediation_result$z0.ci[2]
    ade_p <- mediation_result$z0.p
    
    total_effect <- mediation_result$tau.coef
    total_lower <- mediation_result$tau.ci[1]
    total_upper <- mediation_result$tau.ci[2]
    total_p <- mediation_result$tau.p
    
    prop_mediated <- mediation_result$n0
    prop_lower <- mediation_result$n0.ci[1]
    prop_upper <- mediation_result$n0.ci[2]
    prop_p <- mediation_result$n0.p
    
    df <- rbind(df, data.frame(
      Mediator = mediator,
      ACME_Estimate = acme_est,
      ACME_Lower = acme_lower,
      ACME_Upper = acme_upper,
      ACME_p_value = acme_p,
      ADE_Estimate = ade_est,
      ADE_Lower = ade_lower,
      ADE_Upper = ade_upper,
      ADE_p_value = ade_p,
      Total_Effect = total_effect,
      TE_Lower = total_lower,
      TE_Upper = total_upper,
      TE_p_value = total_p,
      Prop_Mediated = prop_mediated,
      PM_Lower = prop_lower,
      PM_Upper = prop_upper,
      PM_p_value = prop_p
    ))
  }
  com_df <- df %>%
    inner_join(cpd_info %>% dplyr::select(Name, unique_features), by=c("Mediator" = "unique_features")) %>%
    mutate(ACME_fdr = p.adjust(ACME_p_value, method='fdr'),
           ADE_fdr = p.adjust(ADE_p_value, method='fdr'),
           TE_fdr = p.adjust(TE_p_value, method='fdr'),
           PM_fdr = p.adjust(PM_p_value, method='fdr'))
  write.csv(com_df, '20241107_reanalysis/brwe_modeling/bmi_metabolites_bw_mediation.csv', row.names = FALSE)
  filtered_df <- com_df %>%
    filter(ACME_p_value < 0.05)
}
