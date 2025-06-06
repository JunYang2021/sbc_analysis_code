library(tidyverse)
library(readxl)
library(mediation)

{ 
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
  colnames(merge_data)[513] <- 'PSI'
}

{
  # mediation: folic acid intake -> metabolites -> VSI
  c_merge_data <- merge_data %>%
    filter(!is.na(yzhdc_FOLIC1)) # 867
  metabolites <- colnames(c_merge_data)[2:500]
  # y_list <- c('VCI', 'VSI', 'FRI', 'WMI', 'PSI', 'FSIQ')
  # y_list <- c('VCI','FRI', 'WMI', 'PSI', 'FSIQ')
  y_list <- c('PSI', 'FSIQ')
  
  for (y in y_list) {
    x <- 'yzhdc_FOLIC1'
    cat(x, y)
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
      mediator_model <- lm(as.formula(form_mediator), data= c_merge_data)
      outcome_model <- lm(as.formula(form_outcome), data = c_merge_data)
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
    write.csv(com_df, paste0('20241107_reanalysis/iq_modeling/folic_acid_metabolites_', y, '_mediation.csv'), row.names = FALSE)
  }
  
  filtered_df <- com_df %>%
    filter(ACME_p_value < 0.05)
}

{
  # mediation: pre-pregnancy BMI-> metabolites -> VSI
  c_merge_data <- merge_data %>%
    filter(!is.na(before_bmi)) # 1122
  metabolites <- colnames(c_merge_data)[2:500]
  # y_list <- c('VCI', 'VSI', 'FRI', 'WMI', 'PSI', 'FSIQ')
  y_list <- c('VCI','FRI', 'WMI', 'PSI', 'FSIQ')
  for (y in y_list){
    x <- 'before_bmi'
    cat(x, y)
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
      mediator_model <- lm(as.formula(form_mediator), data= c_merge_data)
      outcome_model <- lm(as.formula(form_outcome), data = c_merge_data)
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
    write.csv(com_df, paste0('20241107_reanalysis/iq_modeling/bmi_metabolites_', y, '_mediation.csv'), row.names = FALSE)
  }
  
  filtered_df <- com_df %>%
    filter(ACME_p_value < 0.05)
}

{
  # mediation: multivitamin intake-> metabolites -> VCI
  c_merge_data <- merge_data %>%
    filter(!is.na(yqyzdc_folv4_1)) # 966
  c_merge_data$yqyzdc_folv4_1 <- factor(c_merge_data$yqyzdc_folv4_1,
                                     levels = c(3, 2, 1),
                                     labels =c('Sometimes','Sometimes', 'Everyday'))
  metabolites <- colnames(c_merge_data)[2:500]
  # y_list <- c('VCI', 'VSI', 'FRI', 'WMI', 'PSI', 'FSIQ')
  y_list <- c('VSI','FRI', 'WMI', 'PSI', 'FSIQ')
  for (y in y_list) {
    x <- 'yqyzdc_folv4_1'
    cat(x, y)
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
      mediator_model <- lm(as.formula(form_mediator), data= c_merge_data)
      outcome_model <- lm(as.formula(form_outcome), data = c_merge_data)
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
    write.csv(com_df, paste0('20241107_reanalysis/iq_modeling/multivitamin_metabolites_',y, '_mediation.csv'), row.names = FALSE)
  }
  
  filtered_df <- com_df %>%
    filter(ACME_p_value < 0.05)
}


{
  # Sankey plot: Organize data
  sankey_data <- data.frame(
    exposure = character(),
    mediator = character(),
    outcome = character(),
    `ACME (95% CI)` = character(),
    ACME_p = numeric(),  # float
    `TE (95% CI)` = character(),
    TE_p = numeric(),  # float
    `TE (95% CI)` = character(),
    TE_p = numeric(),  # float
    `PM (95% CI)` = character(),
    PM_p = numeric(),
    stringsAsFactors = FALSE
  )
  
  x_list <- c('bmi', 'multivitamin', 'folic_acid')
  y_list <- c('VCI', 'VSI', 'FRI', 'WMI', 'PSI', 'FSIQ')
  for (x in x_list) {
    for (y in y_list) {
      data_file <- paste0('20241107_reanalysis/iq_modeling/', x, '_metabolites_', y, '_mediation.csv')
      
      res_df <- read.csv(data_file)
      cur_m <- res_df[res_df$ACME_p_value < 0.05 & res_df$TE_p_value < 0.05, 'Mediator']
      cur_acme <- res_df[res_df$ACME_p_value < 0.05 & res_df$TE_p_value < 0.05, 'ACME_Estimate']
      cur_acme_l <- res_df[res_df$ACME_p_value < 0.05 & res_df$TE_p_value < 0.05, 'ACME_Lower']
      cur_acme_u <- res_df[res_df$ACME_p_value < 0.05 & res_df$TE_p_value < 0.05, 'ACME_Upper']
      cur_acme <- paste0(round(cur_acme, 2), 
                         " (", round(cur_acme_l, 2), ", ", round(cur_acme_u, 2), ")")
      cur_acmep <- res_df[res_df$ACME_p_value < 0.05 & res_df$TE_p_value < 0.05, 'ACME_p_value']
      
      cur_te <- res_df[res_df$ACME_p_value < 0.05 & res_df$TE_p_value < 0.05, 'Total_Effect']
      cur_te_l <- res_df[res_df$ACME_p_value < 0.05 & res_df$TE_p_value < 0.05, 'TE_Lower']
      cur_te_u <- res_df[res_df$ACME_p_value < 0.05 & res_df$TE_p_value < 0.05, 'TE_Upper']
      cur_te <- paste0(round(cur_te, 2), 
                       " (", round(cur_te_l, 2), ", ", round(cur_te_u, 2), ")")
      cur_tep <- res_df[res_df$ACME_p_value < 0.05 & res_df$TE_p_value < 0.05, 'TE_p_value']
      
      cur_pm <- res_df[res_df$ACME_p_value < 0.05 & res_df$TE_p_value < 0.05, 'Prop_Mediated']
      cur_pm_l <- res_df[res_df$ACME_p_value < 0.05 & res_df$TE_p_value < 0.05, 'PM_Lower']
      cur_pm_u <- res_df[res_df$ACME_p_value < 0.05 & res_df$TE_p_value < 0.05, 'PM_Upper']
      cur_pm <- paste0(round(cur_pm, 2), 
                       " (", round(cur_pm_l, 2), ", ", round(cur_pm_u, 2), ")")
      cur_pmp <- res_df[res_df$ACME_p_value < 0.05 & res_df$TE_p_value < 0.05, 'PM_p_value']
      
      if (length(cur_m) > 0) {
        cur_x <- rep(x, length(cur_m)) 
        cur_y <- rep(y, length(cur_m)) 
        sankey_data <- rbind(sankey_data, data.frame(
          exposure = cur_x,
          mediator = cur_m,
          outcome = cur_y,
          `ACME (95% CI)` = cur_acme,
          ACME_p = cur_acmep,
          `TE (95% CI)` = cur_te,
          TE_p = cur_tep,
          `PM (95% CI)` = cur_pm,
          PM_p = cur_pmp,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  sankey_data <- sankey_data %>%
    inner_join(cpd_info %>% dplyr::select(Name, unique_features, Class_med), by=c("mediator" = "unique_features")) %>%
    rename(m_name=Name, m_class = Class_med)
  write.csv(sankey_data,
            file="20241107_reanalysis/iq_modeling/mediation_analysis_intelligence_outcome_sankey.csv", row.names = FALSE)
}

{
  # Sankey plot: Make plot
  library(ggsankey)
  sankey_data <- read.csv("20241107_reanalysis/iq_modeling/mediation_analysis_intelligence_outcome_sankey.csv")
  
  sankey_data <- sankey_data %>%
    mutate(exposure = factor(exposure, levels = c('bmi', 'multivitamin', 'folic_acid'),
                             labels=c('Pre-pregnancy BMI',
                                      'Multivitamin intake',
                                      'Folic acid intake')),
           m_class = factor(m_class, levels = c('Amino acids',
                                                'Bile acids',
                                                'Carnitines',
                                                'Fatty acids',
                                                'Fatty acyls',
                                                'NAEs',
                                                'Nucleotides',
                                                'Organic acids',
                                                'Steroids',
                                                'PAs',
                                                'PCs',
                                                'PEs',
                                                'PIs',
                                                'Sphingolipids',
                                                'Xenobiotics',
                                                'Others')),
           outcome = factor(outcome, levels= c('FSIQ', 'VCI', 'VSI',
                                               'FRI', 'WMI', 'PSI'))) %>%
    rename(`Prenatal factors` = exposure, Metabolites = m_class, `Cognitive development` = outcome)
  
  # Prepare the data for Sankey using ggsankey::make_long()
  sankey_data_long <- sankey_data %>%
    ggsankey::make_long(`Prenatal factors`, Metabolites, `Cognitive development`)
  
  # Inspect the long data structure (for understanding)
  head(sankey_data_long)
  
  # Create a Sankey plot using ggplot2
  ggplot(sankey_data_long, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
    geom_sankey(flow.alpha = 0.5, node.color = "gray30") +  # Adjust flow and node colors
    geom_sankey_label(size = 3, color = "black",fill = "white", hjust = 0) +  # Labels for nodes with text box
    theme_sankey(base_size = 30) +  # Theme customization for better readability
    theme(legend.position = "none",  # Remove the legend for node color
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 11,color='black',face='bold', vjust=12), 
          axis.ticks.x = element_blank()) +
    scale_fill_viridis_d(option = "D")
  ggsave('../article/pictures/20241107_reorganize/sankey_mediation.png',
         width = 8, height =6, units = c("in"), dpi = 300, bg='white')
}
