library(tidyverse)
library(readxl)
library(mice)
library(patchwork)
library(ggrepel) # For smart labeling

# Perform linear regression in Batch 2 (Add gender interaction term)
cpd_data <- read.csv('../DATA/peakmat_features/batch2_reorganize_20241107/batch2_compound_data.csv')
cpd_info <- read_excel('../DATA/statistics/identified_compounds_information_batch2.xlsx', sheet=1)
cl_data <- read.csv('../DATA/statistics/clinical_data_batch2.csv')

cl_data$VSI[is.na(cl_data$VSI)] <- mean(cl_data$VSI, na.rm=TRUE)
cl_data$WMI[is.na(cl_data$WMI)] <- mean(cl_data$WMI, na.rm=TRUE)
cl_data$PSIRevised[is.na(cl_data$PSIRevised)] <- mean(cl_data$PSIRevised, na.rm=TRUE)

# Gender interaction term
gender_interaction_term <- data.frame(term=sapply(cpd_info$unique_features, function(x) paste0(x, ":chsex_new2")),
                                          term_name=sapply(cpd_info$Name, function(x) paste0(x, ":gender")))

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

imp <- mice(merge_data, meth = meth, pred=pred, print=F, m=10, seed=2025)

compound_names <- names(merge_data)[2:500]
confounder_names <- names(merge_data)[c(501:504,506,507, 516)]

dependent_vars <- c("VCI", "VSI", "FRI", "WMI", "PSI", "FSIQ")
pol_sheets <- list()

for (var in dependent_vars) {
  pol_sheet <- data.frame()
  for (compound_var in compound_names) {
    interaction_term <- paste(compound_var, '*', 'chsex_new')
    if (var == "PSI") {
      formula <- paste("PSIRevised", "~", paste(c(compound_var,interaction_term, confounder_names), collapse = " + "))
    }
    else {
      formula <- paste(var, "~", paste(c(compound_var, interaction_term, confounder_names), collapse = " + "))
    }
    fit <- with(imp, lm(as.formula(formula)))
    pooled_result <- pool(fit)
    pol <- summary(pooled_result, conf.int=TRUE) %>%
      filter(term==paste(compound_var, ":chsex_new2", sep = ''))
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
  inner_join(gender_interaction_term, by='term')
write.csv(merge_pol_sheet, './20241107_reanalysis/mlr_compounds_batch2_gender_interacted.csv', row.names = FALSE)




{
  # Stratified analysis on gender
  cpd_data <- read.csv('../DATA/peakmat_features/batch2_reorganize_20241107/batch2_compound_data.csv')
  cpd_info <- read_excel('../DATA/statistics/identified_compounds_information_batch2.xlsx', sheet=1)
  cl_data <- read.csv('../DATA/statistics/clinical_data_batch2.csv')
  
  cl_data$VSI[is.na(cl_data$VSI)] <- mean(cl_data$VSI, na.rm=TRUE)
  cl_data$WMI[is.na(cl_data$WMI)] <- mean(cl_data$WMI, na.rm=TRUE)
  cl_data$PSIRevised[is.na(cl_data$PSIRevised)] <- mean(cl_data$PSIRevised, na.rm=TRUE)
  
  boy_cl_data <- cl_data %>%
    filter(chsex_new == 1)
  summary(boy_cl_data)
  girl_cl_data <- cl_data %>%
    filter(chsex_new == 2)
  summary(girl_cl_data)
  
  merge_data <- cpd_data %>%
    inner_join(cl_data, by='sample_name') %>%
    inner_join(birth_data, by='sample_name') %>%
    select(-c(fmbs_CBH46, yqyzdc_income))   # 收入缺失值太多
  cols_to_factor <- c('chsex_new', 'yqyzdc_educton', 'yqyzdc_folv4_1', 'yzhdc_FOLIC1',
                      'yzhdc_SMKN3')
  merge_data[cols_to_factor] <- lapply(merge_data[cols_to_factor], as.factor)
  merge_data$yqyzdc_educton <- factor(merge_data$yqyzdc_educton,
                                      levels = c(2, 3, 4, 5, 6),
                                      labels = c('Low', 'Low','Low','Medium', 'High'))
  
  str(merge_data)
  # gender stratification
  boy_merge_data <- merge_data %>%
    filter(chsex_new == 1)
  girl_merge_data <- merge_data %>%
    filter(chsex_new == 2)
  
  # Multiple imputation for confounders
  imp_ini <- mice(boy_merge_data, m = 1, print=F, maxit = 1, method = "pmm")
  pred <- imp_ini$pred
  pred
  pred[, 1:500] <- 0
  meth <- imp_ini$meth
  meth
  meth[c('yzhdc_FOLIC1', 'yzhdc_SMKN3')] <- 'logreg'
  meth[c('yqyzdc_educton', 'yqyzdc_folv4_1')] <- 'polr'
  
  imp_boy <- mice(boy_merge_data, meth = meth, pred=pred, print=F, m=10, seed=2026)
  plot(imp_boy)
  stripplot(imp_boy, yqyzdc_educton~.imp, pch=5, cex=1)
  
  imp_girl <- mice(girl_merge_data, meth = meth, pred=pred, print=F, m=10, seed=2026)
  plot(imp_girl)
  
  compound_names <- names(merge_data)[2:500]
  confounder_names <- names(merge_data)[c(502:503,505,506, 515)]
  
  dependent_vars <- c("VCI", "VSI", "FRI", "WMI", "PSI", "FSIQ")
  boy_pol_sheets <- list()
  
  for (var in dependent_vars) {
    pol_sheet <- data.frame()
    for (compound_var in compound_names) {
    if (var == "PSI") {
      formula <- paste("PSIRevised", "~", paste(c(compound_var, confounder_names), collapse = " + "))
    }
    else {
      formula <- paste(var, "~", paste(c(compound_var, confounder_names), collapse = " + "))
    }
    fit <- with(imp_boy, lm(as.formula(formula)))
    pooled_result <- pool(fit)
    pol <- summary(pooled_result, conf.int=TRUE) %>%
      filter(term==compound_var) %>%
      mutate(z = estimate / std.error)
    pol_sheet <- rbind(pol_sheet, pol)
    }
    boy_pol_sheets[[var]] <- pol_sheet
  }
  
  boy_merge_pol_sheet <- boy_pol_sheets[['VCI']] %>%
    select(term, estimate, p.value, `2.5 %`, `97.5 %`, z) %>%
    rename('VCI estimate' = estimate, 'VCI p' = p.value, 'VCI lower' = `2.5 %`, 'VCI upper' = `97.5 %`, 'VCI z' = z) %>%
    inner_join(boy_pol_sheets[['VSI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`, z) 
               %>% rename('VSI estimate' = estimate, 'VSI p' = p.value, 'VSI lower' = `2.5 %`, 'VSI upper' = `97.5 %`, 'VSI z' = z), by='term') %>%
    inner_join(boy_pol_sheets[['FRI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`, z) 
               %>% rename('FRI estimate' = estimate, 'FRI p' = p.value, 'FRI lower' = `2.5 %`, 'FRI upper' = `97.5 %`, 'FRI z' = z), by='term') %>%
    inner_join(boy_pol_sheets[['WMI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`, z) 
               %>% rename('WMI estimate' = estimate, 'WMI p' = p.value, 'WMI lower' = `2.5 %`, 'WMI upper' = `97.5 %`, 'WMI z' = z), by='term') %>%
    inner_join(boy_pol_sheets[['PSI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`, z) 
               %>% rename('PSI estimate' = estimate, 'PSI p' = p.value, 'PSI lower' = `2.5 %`, 'PSI upper' = `97.5 %`, 'PSI z' = z), by='term') %>%
    inner_join(boy_pol_sheets[['FSIQ']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`, z) 
               %>% rename('FSIQ estimate' = estimate, 'FSIQ p' = p.value, 'FSIQ lower' = `2.5 %`, 'FSIQ upper' = `97.5 %`, 'FSIQ z' = z), by='term') %>%
    inner_join(cpd_info %>% select(unique_features, Name), by = c('term' = 'unique_features'))

  
  girl_pol_sheets <- list()
  
  for (var in dependent_vars) {
    pol_sheet <- data.frame()
    for (compound_var in compound_names) {
      if (var == "PSI") {
        formula <- paste("PSIRevised", "~", paste(c(compound_var, confounder_names), collapse = " + "))
      }
      else {
        formula <- paste(var, "~", paste(c(compound_var, confounder_names), collapse = " + "))
      }
      fit <- with(imp_girl, lm(as.formula(formula)))
      pooled_result <- pool(fit)
      pol <- summary(pooled_result, conf.int=TRUE) %>%
        filter(term==compound_var) %>%
        mutate(z = estimate / std.error)
      pol_sheet <- rbind(pol_sheet, pol)
    }
    girl_pol_sheets[[var]] <- pol_sheet
  }
  
  girl_merge_pol_sheet <- girl_pol_sheets[['VCI']] %>%
    select(term, estimate, p.value, `2.5 %`, `97.5 %`, z) %>%
    rename('VCI estimate' = estimate, 'VCI p' = p.value, 'VCI lower' = `2.5 %`, 'VCI upper' = `97.5 %`, 'VCI z' = z) %>%
    inner_join(girl_pol_sheets[['VSI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`, z) 
               %>% rename('VSI estimate' = estimate, 'VSI p' = p.value, 'VSI lower' = `2.5 %`, 'VSI upper' = `97.5 %`, 'VSI z' = z), by='term') %>%
    inner_join(girl_pol_sheets[['FRI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`, z) 
               %>% rename('FRI estimate' = estimate, 'FRI p' = p.value, 'FRI lower' = `2.5 %`, 'FRI upper' = `97.5 %`, 'FRI z' = z), by='term') %>%
    inner_join(girl_pol_sheets[['WMI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`, z) 
               %>% rename('WMI estimate' = estimate, 'WMI p' = p.value, 'WMI lower' = `2.5 %`, 'WMI upper' = `97.5 %`, 'WMI z' = z), by='term') %>%
    inner_join(girl_pol_sheets[['PSI']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`, z) 
               %>% rename('PSI estimate' = estimate, 'PSI p' = p.value, 'PSI lower' = `2.5 %`, 'PSI upper' = `97.5 %`, 'PSI z' = z), by='term') %>%
    inner_join(girl_pol_sheets[['FSIQ']] %>% select(term, estimate, p.value, `2.5 %`, `97.5 %`, z) 
               %>% rename('FSIQ estimate' = estimate, 'FSIQ p' = p.value, 'FSIQ lower' = `2.5 %`, 'FSIQ upper' = `97.5 %`, 'FSIQ z' = z), by='term') %>%
    inner_join(cpd_info %>% select(unique_features, Name), by = c('term' = 'unique_features'))
  
  boy_girl_str_sheet <- bind_rows(
    boy_merge_pol_sheet %>% mutate(`Fetal sex` = 'Boy'),
    girl_merge_pol_sheet %>% mutate(`Fetal sex` = 'Girl')
  )
  write.csv(boy_girl_str_sheet, file="./20241107_reanalysis/multiple_regression_batch2_gender_stratification_20241128.csv", row.names = FALSE)
}



{
  # Plot scatter plots of stratified analysis
  boy_girl_str_sheet <- read.csv("./20241107_reanalysis/multiple_regression_batch2_gender_stratification_20241128.csv", check.names = FALSE)
    
  
  iq_items <- c('VCI', 'VSI', 'FRI', 'WMI', 'PSI', 'FSIQ')
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
  for (iq_item in iq_items) {
    iq_p <- paste(iq_item, 'p')
    iq_z <- paste(iq_item, 'z')
    boy_iq_p <- paste('Boy', iq_item, 'p')
    boy_iq_z <- paste('Boy', iq_item, 'z')
    girl_iq_p <- paste('Girl', iq_item, 'p')
    girl_iq_z <- paste('Girl', iq_item, 'z')
    
    data_stratified <- boy_girl_str_sheet %>%
      select(term, Name, !!sym(iq_p), !!sym(iq_z), `Fetal sex`) %>%
      pivot_wider(
        names_from = `Fetal sex`,
        values_from = c(iq_p, iq_z),
        names_glue = "{`Fetal sex`} {.value}"
      ) %>%
      mutate(
        Shape = case_when(
          !!sym(boy_iq_p) < 0.05 & !!sym(girl_iq_p) < 0.05 ~ "Both significant",
          !!sym(boy_iq_p) < 0.05 ~ "Significant in males",
          !!sym(girl_iq_p) < 0.05 ~ "Significant in females",
          TRUE ~ "Not significant"
        )
      ) %>%
        mutate(
          Quadrant = factor(case_when(
            !!sym(boy_iq_z) > 0 & !!sym(girl_iq_z) > 0 ~ 'I',
            !!sym(boy_iq_z) < 0 & !!sym(girl_iq_z) > 0 ~ 'II',
            !!sym(boy_iq_z) < 0 & !!sym(girl_iq_z) < 0 ~ 'III',
            !!sym(boy_iq_z) > 0 & !!sym(girl_iq_z) < 0 ~ 'IV',
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
                                      'Others'))) %>%
      mutate(Name = if_else(Name == '(3b,4b,11b,14b)-11-Ethoxy-3,4-epoxy-14-hydroxy-12-cyathen-15-al 14-xyloside',
                            'CHEBI:168534', Name))
  print(iq_item)
  print(summary(data_stratified$Quadrant))
  
  ggplot(data_stratified, aes(x = !!sym(boy_iq_z), y = !!sym(girl_iq_z), shape = Shape, color=Class_pie)) +
    geom_point(size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "gray") +
    geom_hline(yintercept = 0, color = "black") +
    geom_vline(xintercept = 0, color = "black") +
    geom_text_repel(aes(label = Name), size = 2, max.overlaps = 9, color='black',
                    min.segment.length = Inf) + # Add labels
    # geom_text(check_overlap = TRUE, color='black', size=2)+
    scale_shape_manual(values = c('Both significant' = 16,
                                  'Significant in males'=17,
                                  'Significant in females'=15,
                                  'Not significant'=18)) +
    scale_color_manual(values = class_palette) +
    labs(
      x = paste("Z-scores of",  iq_item, "(male offspring)"),
      y = paste("Z-scores of",  iq_item, "(female offspring)"),
      shape = "Significance",
      color='Class'
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      text = element_text(size = 12),
      panel.border = element_blank()
    )
  ggsave(paste("../article/pictures/20241107_reorganize/stratified/", iq_item, '.png', sep=''),
         width = 6,
         height =5,
         units = c("in"),
         dpi = 300,
         bg='white')
  } 
  
}


{
  # Gender interacted variables
  sig_mlr_estimates <- read.csv('./20241107_reanalysis/mlr_compounds_batch2_20250307.csv',
                                check.names = FALSE) %>%
    filter(`VCI p` < 0.05 | `VSI p` < 0.05 | `PSI p` < 0.05 |`FRI p` < 0.05|`WMI p` < 0.05|`FSIQ p` < 0.05)
  
  sig_interacted <- read.csv('./20241107_reanalysis/mlr_compounds_batch2_gender_interacted.csv',
                              check.names = FALSE) %>%
    filter(`VCI p` < 0.05 | `VSI p` < 0.05 | `PSI p` < 0.05 |`FRI p` < 0.05|`WMI p` < 0.05|`FSIQ p` < 0.05) %>%
    mutate(feature = gsub(':chsex_new2', "", term)) %>%
    filter(feature %in% sig_mlr_estimates$term)
  
  boy_girl_str_sheet <- read.csv("./20241107_reanalysis/multiple_regression_batch2_gender_stratification_20241128.csv",
                                 check.names = FALSE)
  
  # Dot-whisper plot
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
  
  plot_data <- boy_girl_str_sheet %>%
    filter(term %in% sig_interacted$feature) %>%
    pivot_longer(
      cols = starts_with("VCI") | starts_with("VSI") | starts_with("FRI") | starts_with("WMI") | starts_with("PSI") | starts_with("FSIQ"),
      names_to = c("index", ".value"), 
      names_pattern = "(.*) (.*)"
    ) %>%
    mutate(significance = ifelse(p < 0.05, "Significant", "Not Significant")) %>%
    rename('Index' = 'index', 'Significance' = 'significance') %>%
    inner_join(cpd_info %>% select(unique_features, Class), by=c("term" = "unique_features")) %>%
    group_by(Name) %>%  # Group by the compound
    mutate(y = ifelse(row_number() == 1, 1, 0)) %>% # Add 'y', 1 for the first row, 0 for others
    ungroup() %>%
    arrange(Class) %>%  # First sort by 'effect' and then by 'Class'
    mutate(Name = factor(Name, levels = unique(Name)))  %>%
    mutate(Index = factor(Index, levels = c("VCI", "VSI", "FRI", "WMI", "PSI", "FSIQ")))
  
  dtplot <- ggplot(plot_data, aes(x = Name, y = estimate, shape = Index, color = `Fetal sex`)) +
    geom_pointrange(aes(ymin = lower, ymax = upper), 
                    position = position_dodge(width = 0.6)) +  # Separate dependent variables to avoid overlap
    geom_point(size = 0.01, position = position_dodge(width = 0.6)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("Boy" = "#1b9e77", "Girl" = "#7570b3")) +  # Red for p < 0.05
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
          axis.title.y = element_blank(),
          strip.text = element_blank())
  
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
          plot.background = element_blank())
  
  p_y + dtplot +
    plot_layout(widths = c(.07, 1), guides='collect') &
    theme()
  ggsave("../article/pictures/20241107_reorganize/stratified/dot_whisker_plot_batch2_gender_significant_cpd_stra.png", width = 10, height =11, units = c("in"), dpi = 300, bg='white')
}

{
  # Gender interacted variables （VCI）
  sig_interacted <- read.csv('./20241107_reanalysis/mlr_compounds_batch2_gender_interacted.csv',
                             check.names = FALSE) %>%
    filter(`VCI p` < 0.05) %>%
    mutate(feature = gsub(':chsex_new2', "", term))
  
  boy_girl_str_sheet <- read.csv("./20241107_reanalysis/multiple_regression_batch2_gender_stratification_20241128.csv",
                                 check.names = FALSE)
  
  # Dot-whisper plot
  class_palette <- c("Amino acids" = '#a6cee3',
                     "Bile acids" = '#1f78b4',
                     "Carnitines" = '#b2df8a',
                     "Fatty acids" = '#33a02c',
                     "Fatty acyls" = '#fb9a99',
                     "Glycerophosphates" = '#d9d9d9',
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
  
  plot_data <- boy_girl_str_sheet %>%
    filter(term %in% sig_interacted$feature) %>%
    mutate(Significance = ifelse(`VCI p` < 0.05, "Significant", "Not Significant")) %>%
    inner_join(cpd_info %>% select(unique_features, Class), by=c("term" = "unique_features")) %>%
    group_by(Name) %>%  # Group by the compound
    mutate(effect = if_else(`VCI estimate`[`Fetal sex` == "Boy"] < `VCI estimate`[`Fetal sex` == "Girl"], -1, 1)) %>%
    mutate(y = ifelse(row_number() == 1, 1, 0)) %>% # Add 'y', 1 for the first row, 0 for others
    ungroup() %>%
    arrange(effect, Class, Name) %>%
    mutate(Name = factor(Name, levels = unique(Name))) %>%
    mutate(Class = if_else(Class == 'Glycerophosphoethanolamines', 'PE', Class)) %>%
    mutate(Class = if_else(Class == 'Glycerophosphocholines', 'PC', Class))
  
  dtplot <- ggplot(plot_data, aes(x = Name, y = `VCI estimate`, shape = `Fetal sex`, color = Significance)) +
    geom_pointrange(aes(ymin = `VCI lower`, ymax = `VCI upper`), 
                    position = position_dodge(width = 0.3)) +  # Separate dependent variables to avoid overlap
    geom_point(size = 0.01, position = position_dodge(width = 0.7)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("Not Significant" = "black", "Significant" = "#7570b3")) +
    scale_shape_manual(values = c(16, 17)) + 
    coord_flip() +
    theme_minimal() +
    labs(title = "",
         x = "", 
         y = "Beta Estimate (95% CI)") +
    theme(panel.grid.major.y = element_blank(),
          axis.ticks.y = element_blank(), 
          axis.line.y = element_blank(),
          axis.text.y = element_blank(), 
          axis.title.y = element_blank(),
          strip.text = element_blank())
  
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
          plot.background = element_blank())
  
  p_y + dtplot +
    plot_layout(widths = c(.07, 1), guides='collect') &
    theme()
  ggsave("../article/pictures/20241107_reorganize/stratified/dot_whisker_plot_batch2_gender_stra_VCI.png", width = 10, height =11, units = c("in"), dpi = 300, bg='white')
}

{
  # Gender interacted variables （FSIQ）
  sig_interacted <- read.csv('./20241107_reanalysis/mlr_compounds_batch2_gender_interacted.csv',
                             check.names = FALSE) %>%
    filter(`FSIQ p` < 0.05) %>%
    mutate(feature = gsub(':chsex_new2', "", term))
  
  boy_girl_str_sheet <- read.csv("./20241107_reanalysis/multiple_regression_batch2_gender_stratification_20241128.csv",
                                 check.names = FALSE)
  
  # Dot-whisper plot
  class_palette <- c("Amino acids" = '#a6cee3',
                     "Bile acids" = '#1f78b4',
                     "Carnitines" = '#b2df8a',
                     "Fatty acids" = '#33a02c',
                     "Fatty acyls" = '#fb9a99',
                     "Glycerophosphates" = '#d9d9d9',
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
  
  plot_data <- boy_girl_str_sheet %>%
    filter(term %in% sig_interacted$feature) %>%
    mutate(Significance = ifelse(`FSIQ p` < 0.05, "Significant", "Not Significant")) %>%
    inner_join(cpd_info %>% select(unique_features, Class), by=c("term" = "unique_features")) %>%
    group_by(Name) %>%  # Group by the compound
    mutate(effect = if_else(`FSIQ estimate`[`Fetal sex` == "Boy"] < `FSIQ estimate`[`Fetal sex` == "Girl"], -1, 1)) %>%
    mutate(y = ifelse(row_number() == 1, 1, 0)) %>% # Add 'y', 1 for the first row, 0 for others
    ungroup() %>%
    arrange(effect, Class, Name) %>%
    mutate(Name = factor(Name, levels = unique(Name))) %>%
    mutate(Class = if_else(Class == 'Glycerophosphoethanolamines', 'PE', Class)) %>%
    mutate(Class = if_else(Class == 'Glycerophosphocholines', 'PC', Class))
  
  dtplot <- ggplot(plot_data, aes(x = Name, y = `FSIQ estimate`, shape = `Fetal sex`, color = Significance)) +
    geom_pointrange(aes(ymin = `FSIQ lower`, ymax = `FSIQ upper`), 
                    position = position_dodge(width = 0.3)) +  # Separate dependent variables to avoid overlap
    geom_point(size = 0.01, position = position_dodge(width = 0.7)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("Not Significant" = "black", "Significant" = "#7570b3")) +
    scale_shape_manual(values = c(16, 17)) + 
    coord_flip() +
    theme_minimal() +
    labs(title = "",
         x = "", 
         y = "Beta Estimate (95% CI)") +
    theme(panel.grid.major.y = element_blank(),
          axis.ticks.y = element_blank(), 
          axis.line.y = element_blank(),
          axis.text.y = element_blank(), 
          axis.title.y = element_blank(),
          strip.text = element_blank())
  
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
          plot.background = element_blank())
  
  p_y + dtplot +
    plot_layout(widths = c(.07, 1), guides='collect') &
    theme()
  ggsave("../article/pictures/20241107_reorganize/stratified/dot_whisker_plot_batch2_gender_stra_FSIQ.png", width = 10, height =14, units = c("in"), dpi = 300, bg='white')
}


{
  # Gender interacted variables
  boy_girl_str_sheet <- read.csv("./20241107_reanalysis/multiple_regression_batch2_gender_stratification_20241128.csv",
                                 check.names = FALSE) %>%
    mutate(`Fetal sex` = ifelse(`Fetal sex` == "Boy", "Male",
                                ifelse(`Fetal sex` == "Girl", "Female", `Fetal sex`)))
  
  # Dot-whisper plot
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
  
  for (iq_var in c('VCI', 'VSI', 'WMI', 'FRI', 'PSI', 'FSIQ')) {
    iq_p <- paste(iq_var, 'p')
    iq_estimate <- paste(iq_var, 'estimate')
    iq_l <- paste(iq_var, 'lower')
    iq_u <- paste(iq_var, 'upper')
    
    sig_interacted <- read.csv('./20241107_reanalysis/mlr_compounds_batch2_gender_interacted.csv',
                             check.names = FALSE) %>%
      filter(!!sym(iq_p) < 0.05) %>%
      mutate(feature = gsub(':chsex_new2', "", term))
    
    plot_data <- boy_girl_str_sheet %>%
      filter(term %in% sig_interacted$feature) %>%
      mutate(Significance = ifelse(get(iq_p) < 0.05, "Significant", "Not Significant")) %>%
      inner_join(cpd_info %>% select(unique_features, Class), by=c("term" = "unique_features")) %>%
      group_by(Name) %>%  # Group by the compound
      mutate(effect = if_else(get(iq_estimate)[`Fetal sex` == "Male"] < get(iq_estimate)[`Fetal sex` == "Female"], -1, 1)) %>%
      mutate(y = ifelse(row_number() == 1, 1, 0)) %>% # Add 'y', 1 for the first row, 0 for others
      ungroup() %>%
      arrange(effect, Class, Name) %>%
      mutate(Name = factor(Name, levels = unique(Name))) %>%
      mutate(Class = if_else(Class == 'Glycerophosphates', 'PA', Class),
             Class = if_else(Class == 'Glycerophosphocholines', 'PC', Class),
             Class = if_else(Class == 'Glycerophosphoethanolamines', 'PE', Class),
             Class = if_else(Class == 'Glycerophosphoglycerols', 'PG', Class),
             Class = if_else(Class == 'Glycerophosphoinositols', 'PI', Class))
    
    dtplot <- ggplot(plot_data, aes(x = Name, y = !!sym(iq_estimate), shape = `Fetal sex`, color = Significance)) +
      geom_pointrange(aes(ymin = !!sym(iq_l), ymax = !!sym(iq_u)), 
                      position = position_dodge(width = 0.3)) +  # Separate dependent variables to avoid overlap
      geom_point(size = 0.01, position = position_dodge(width = 0.7)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      scale_color_manual(values = c("Not Significant" = "black", "Significant" = "#7570b3")) +
      scale_shape_manual(values = c("Male"=16, "Female"=17)) +  # Assign different shapes for different dependent variables
      coord_flip() +
      theme_minimal() +
      labs(title = "",
           x = "", 
           y = "Beta Estimate (95% CI)") +
      theme(panel.grid.major.y = element_blank(),
            axis.ticks.y = element_blank(), 
            axis.line.y = element_blank(),
            axis.text.y = element_blank(), 
            axis.title.y = element_blank(),
            strip.text = element_blank())
    
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
            plot.background = element_blank())
    
    p_y + dtplot +
      plot_layout(widths = c(.07, 1), guides='collect') &
      theme()
    ggsave(paste0("../article/pictures/20241107_reorganize/stratified/dot_whisker_plot_batch2_gender_interacted_cpd_stra_",
                  iq_var,
                  ".png"),
           width = 10,
           height =11, 
           units = c("in"), 
           dpi = 300, 
           bg='white')
  }
}

{
  # Gender interacted variables for VCI (FIGURE 5)
  boy_girl_str_sheet <- read.csv("./20241107_reanalysis/multiple_regression_batch2_gender_stratification_20241128.csv",
                                 check.names = FALSE)
  
  # Dot-whisper plot
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
  
  iq_var <- 'VCI'
  iq_p <- paste(iq_var, 'p')
  iq_estimate <- paste(iq_var, 'estimate')
  iq_l <- paste(iq_var, 'lower')
  iq_u <- paste(iq_var, 'upper')
  
  sig_interacted <- read.csv('./20241107_reanalysis/mlr_compounds_batch2_gender_interacted.csv',
                             check.names = FALSE) %>%
    filter(!!sym(iq_p) < 0.05) %>%
    mutate(feature = gsub(':chsex_new2', "", term))
  
  plot_data <- boy_girl_str_sheet %>%
    filter(term %in% sig_interacted$feature) %>%
    mutate(Significance = ifelse(get(iq_p) < 0.05, "Significant", "Not Significant")) %>%
    inner_join(cpd_info %>% select(unique_features, Class), by=c("term" = "unique_features")) %>%
    group_by(Name) %>%  # Group by the compound
    mutate(effect = if_else(get(iq_estimate)[`Fetal sex` == "Boy"] < get(iq_estimate)[`Fetal sex` == "Girl"], -1, 1)) %>%
    mutate(y = ifelse(row_number() == 1, 1, 0)) %>% # Add 'y', 1 for the first row, 0 for others
    ungroup() %>%
    arrange(effect, Class, desc(Name)) %>%
    mutate(Name = factor(Name, levels = unique(Name))) %>%
    mutate(Class = if_else(Class == 'Glycerophosphates', 'PA', Class),
           Class = if_else(Class == 'Glycerophosphocholines', 'PC', Class),
           Class = if_else(Class == 'Glycerophosphoethanolamines', 'PE', Class),
           Class = if_else(Class == 'Glycerophosphoglycerols', 'PG', Class),
           Class = if_else(Class == 'Glycerophosphoinositols', 'PI', Class))
  
  dtplot <- ggplot(plot_data, aes(x = Name, y = !!sym(iq_estimate), shape = `Fetal sex`, color = Significance)) +
    geom_pointrange(aes(ymin = !!sym(iq_l), ymax = !!sym(iq_u)), 
                    position = position_dodge(width = 0.3)) +  # Separate dependent variables to avoid overlap
    geom_point(size = 0.01, position = position_dodge(width = 0.7)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("Not Significant" = "black", "Significant" = "#7570b3")) +
    scale_shape_manual(values = c(16, 17)) +  # Assign different shapes for different dependent variables
    theme_minimal() +
    labs(title = "",
         x = "", 
         y = "Beta Estimate (95% CI)") +
    theme(panel.grid.major.x = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color = "black", size=9)
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
  ggsave(paste0("../article/pictures/20241107_reorganize/stratified/dot_whisker_plot_batch2_gender_interacted_cpd_stra_main",
                iq_var,
                ".png"),
         width = 11,
         height =5, 
         units = c("in"), 
         dpi = 300, 
         bg='white')
}