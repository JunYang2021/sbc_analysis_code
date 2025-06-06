# Construct models for IQ
library(tidyverse)
library(readxl)
library(caret)  

{ 
  cpd_data <- read.csv('../DATA/peakmat_features/batch2_reorganize_20241107/batch2_compound_data.csv')
  cpd_info <- read_excel('../DATA/statistics/identified_compounds_information_batch2.xlsx', sheet=1)
  cl_data <- read.csv('../DATA/statistics/clinical_data_batch2.csv')
  
  cl_data$VSI[is.na(cl_data$VSI)] <- mean(cl_data$VSI, na.rm=TRUE)
  cl_data$WMI[is.na(cl_data$WMI)] <- mean(cl_data$WMI, na.rm=TRUE)
  cl_data$PSIRevised[is.na(cl_data$PSIRevised)] <- mean(cl_data$PSIRevised, na.rm=TRUE)
  names(cl_data)[names(cl_data) == 'PSIRevised'] <- 'PSI'
  birth_data <- read.csv('../DATA/statistics/clinical_birth_data_batch2.csv')
  
  
  batch2_data <- read.csv('./20241107_reanalysis/mlr_compounds_batch2_20250307.csv', check.names = FALSE)
  
  merge_data <- cpd_data %>%
    inner_join(cl_data, by='sample_name') %>%
    inner_join(birth_data, by='sample_name')
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
  # Construct machine learning model to get variance explained
  library(caret)
  
  # Get Rsquared from all the factors
  dependent_vars <- c("VCI", "VSI", "FRI", "WMI", "PSI", "FSIQ")
  influence_factors <- c('metabolomics','chsex_new', 'yqyzdc_educton', 'yqyzdc_income',
                         'yqyzdc_folv4_1', 'maternal_age', 'before_bmi', 'yzhdc_FOLIC1',
                         'yzhdc_SMKN3', 'yunzhou')
  
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
    if (cur_fac == 'metabolomics') {
      # Perform elastic net regression
      for (d_var in dependent_vars) {
        print(d_var)
        d_var_p <- paste(d_var, 'p')
        sig_cpds <- batch2_data %>%
          filter(!!sym(d_var_p) < 0.05)
        print(nrow(sig_cpds))
        formula <- as.formula(paste(d_var, "~", paste(sig_cpds$term, collapse = " + ")))
        cv_model <- train(formula,
                          data = merge_data, 
                          method = "glmnet", 
                          trControl = train_control)
        valid_r2 <- mean(cv_model$resample$Rsquared)
        
        n <- nrow(merge_data)
        k <- nrow(sig_cpds)
        adj_r2 <- 1 - (1 - valid_r2) * (n - 1) / (n - k - 1)
        
        train_predictions <- predict(cv_model, newdata = merge_data)
        
        # Actual values
        y_actual <- merge_data[[d_var]]
        
        # Compute R-squared for the training set
        sst <- sum((y_actual - mean(y_actual))^2)  # Total sum of squares
        sse <- sum((y_actual - train_predictions)^2)  # Residual sum of squares
        train_r2 <- 1 - (sse / sst)
        
        # Save results
        results <- rbind(results, data.frame(Factor = "metabolomics", 
                                             DependentVar = d_var,
                                             TrainR2 = train_r2,
                                             ValidR2 = valid_r2,
                                             AdjR2 = adj_r2,
                                             pvalue='-'))
      }
    }
    else {
      # Perform which model (only one IV)
      for (d_var in dependent_vars) {
        print(d_var)
        merge_data_com <- merge_data %>%
          filter(!is.na(get(cur_fac)))
        formula <- as.formula(paste(d_var, "~", cur_fac))
        cv_model <- train(formula, data = merge_data_com, 
                          method = "lm", trControl = train_control)
        valid_r2 <- mean(cv_model$resample$Rsquared)
        
        # Fit the model on the full training data (non-CV)
        full_model <- lm(formula, data = merge_data_com)
        # print(summary(full_model))
        
        # Calculate R2 on training set
        train_predictions <- predict(full_model, merge_data_com)
        train_actuals <- merge_data_com[[d_var]]
        ss_total <- sum((train_actuals - mean(train_actuals))^2)
        ss_residual <- sum((train_actuals - train_predictions)^2)
        train_r2 <- 1 - (ss_residual / ss_total)
        
        p_value <- coef(summary(full_model))[2, "Pr(>|t|)"]
        
        n <- nrow(merge_data_com)
        k <- 1
        adj_r2 <- 1 - (1 - valid_r2) * (n - 1) / (n - k - 1)
        
        # Save results
        results <- rbind(results, data.frame(Factor = cur_fac, 
                                             DependentVar = d_var,
                                             TrainR2 = train_r2,
                                             ValidR2 = valid_r2,
                                             AdjR2 = adj_r2,
                                             pvalue=p_value))
      }
    }
  }
  
  write.csv(results, "./20241107_reanalysis/iq_modeling/adjusted_r2_results_elastic_net.csv",
            row.names = FALSE)
}

{
  # Construct lm model to get correlations (for covarites)
  dependent_vars <- c("VCI", "VSI", "FRI", "WMI", "PSI", "FSIQ")
  influence_factors <- c('chsex_new', 'yqyzdc_educton', 'yqyzdc_income',
                         'yqyzdc_folv4_1', 'maternal_age', 'before_bmi', 'yzhdc_FOLIC1',
                         'yzhdc_SMKN3', 'yunzhou', 'fmbs_CBH46', 'brwe')
  results_list <- list()
  id <- 1
  for (d_var in dependent_vars) {
    for (cur_fac in influence_factors) {
      print(cur_fac)
      # Perform which model (only one IV)
      merge_data_com <- merge_data %>%
        filter(!is.na(get(cur_fac)))
      formula <- as.formula(paste(d_var,"~", cur_fac))
      full_model <- lm(formula, data = merge_data_com)
      print(summary(full_model))
      model_summary <- summary(full_model)
      coefficients <- model_summary$coefficients
      conf_int <- confint(full_model)
      
      # Remove intercept row
      coefficients <- coefficients[-1, , drop = FALSE]
      conf_int <- conf_int[-1, , drop = FALSE]
      
      # Store results
      results_list[[id]] <- data.frame(
        Variable = cur_fac,
        DV = d_var,
        Term = rownames(coefficients),
        Estimate = coefficients[, 1],
        Conf.Low = conf_int[, 1],
        Conf.High = conf_int[, 2],
        P.Value = coefficients[, 4],
        N = nrow(merge_data_com)
      )
      id <- id+1
    
  }
  }
  results <- do.call(rbind, results_list) %>%
    select(-Variable) %>%
    pivot_wider(
      names_from = DV,
      values_from = c(Estimate, Conf.Low, Conf.High, P.Value, N)
    ) %>%
    mutate(VCI_combined = paste0(sprintf("%.3f", Estimate_VCI), " (", sprintf("%.3f", Conf.Low_VCI), ", ", sprintf("%.3f", Conf.High_VCI), ")")) %>%
    mutate(VSI_combined = paste0(sprintf("%.3f", Estimate_VSI), " (", sprintf("%.3f", Conf.Low_VSI), ", ", sprintf("%.3f", Conf.High_VSI), ")")) %>%
    mutate(FRI_combined = paste0(sprintf("%.3f", Estimate_FRI), " (", sprintf("%.3f", Conf.Low_FRI), ", ", sprintf("%.3f", Conf.High_FRI), ")")) %>%
    mutate(WMI_combined = paste0(sprintf("%.3f", Estimate_WMI), " (", sprintf("%.3f", Conf.Low_WMI), ", ", sprintf("%.3f", Conf.High_WMI), ")")) %>%
    mutate(PSI_combined = paste0(sprintf("%.3f", Estimate_PSI), " (", sprintf("%.3f", Conf.Low_PSI), ", ", sprintf("%.3f", Conf.High_PSI), ")")) %>%
    mutate(FSIQ_combined = paste0(sprintf("%.3f", Estimate_FSIQ), " (", sprintf("%.3f", Conf.Low_FSIQ), ", ", sprintf("%.3f", Conf.High_FSIQ), ")"))

  write.csv(results, "./20241107_reanalysis/iq_modeling/covariates_iq_lm.csv",
            row.names = FALSE)
  
  
}

{
  # Plot the variance explained
  formal_names <- c("chsex_new" = "Child gender", 
                    "yqyzdc_educton" = "Maternal education", 
                    "yqyzdc_income" = "Family income", 
                    "yqyzdc_folv4_1" = "Multivitamin intake", 
                    "maternal_age" = "Maternal age", 
                    "before_bmi" = "Pre-pregnancy BMI", 
                    "yzhdc_FOLIC1" = "Folic acid intake", 
                    "yzhdc_SMKN3" = "Passive smoking", 
                    "yunzhou" = "Gestational age", 
                    "metabolomics" = "Metabolomics")
  r2_results <- read.csv("./20241107_reanalysis/iq_modeling/adjusted_r2_results_elastic_net.csv")
  r2_results$Factor <- formal_names[r2_results$Factor]
  
  r2_results$Factor <- factor(r2_results$Factor,
                              levels = c("Metabolomics","Maternal education",
                                         "Family income", "Child gender", 
                                         "Multivitamin intake","Folic acid intake", 
                                         "Passive smoking",
                                         "Maternal age","Gestational age","Pre-pregnancy BMI"))
  r2_results$DependentVar <- factor(r2_results$DependentVar, levels = c("VCI", "VSI",
                                                                        "WMI", "FRI",
                                                                        "PSI", "FSIQ"))
  
  ggplot(r2_results, aes(x = Factor, y = ValidR2 * 100, fill = DependentVar)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
    geom_text(aes(label = round(ValidR2 * 100, 1)), 
              position = position_dodge(width = 0.8), vjust = -0.5, size = 2.2) +
    scale_fill_brewer(palette='Dark2') +
    labs(x = NULL, y = "R-squared (%)", fill = NULL) + 
    theme_minimal() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = 'black'),
          legend.position = c(0.85, 0.75),       # Move legend to upper middle
          legend.justification = "center")   # Center legend
  ggsave('../article/pictures/20241107_reorganize/cv/intelligence_elastic_net.png',
         width = 9, height=5, units = 'in', dpi = 300)
}


{
  dependent_vars <- c("VCI", "VSI", "FRI", "WMI", "PSI", "FSIQ")
  for (d_var in dependent_vars) {
    train_control <- trainControl(method = "repeatedcv",      # Cross-validation method
                                number = 10,                # Number of folds
                                repeats = 20,              # Number of repeats
                                savePredictions = TRUE)
    d_var_p <- paste(d_var, 'p')
    sig_cpds <- batch2_data %>%
      filter(!!sym(d_var_p) < 0.05)
    formula <- as.formula(paste(d_var, "~", paste(sig_cpds$term, collapse = " + ")))
    cv_model <- train(formula,
                      data = merge_data, 
                      method = "glmnet", 
                      trControl = train_control)
    r2_values <- cv_model$resample$Rsquared
    ggplot(data = data.frame(R2 = r2_values), aes(x = R2)) +
      geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black", alpha = 0.7) +
      labs(x = "R-squared", y = "Frequency") +
      theme_classic()
    
    imp<- varImp(cv_model)
    plot(imp)
    
    imp_data <- as.data.frame(varImp(cv_model)$importance) %>%
      mutate(Name = sig_cpds$Name)
    
    imp_data_fig <- imp_data %>%
      filter(Overall != 0) %>%
      arrange(desc(Overall))
    
    
    # Create the plot
    ggplot(imp_data_fig, aes(x = reorder(Name, Overall), y = Overall)) +
      geom_segment(aes(xend = Name, y = 0, yend = Overall), color = "steelblue", linewidth=1) +
      geom_point(color = "red", size = 3) +  # Add a red dot at the end of each stick
      coord_flip() +  # Flip the axes for better visualization
      labs(
        x = "",
        y = "Importance"
      ) +
      theme_classic() +
      theme(axis.text.x = element_text(color = 'black'),
            axis.text.y = element_text(color = 'black'))
    
    if (d_var == 'VCI') {
      ggsave(paste0('../article/pictures/20241107_reorganize/cv/', d_var, '_elastic_net.png'),
           width = 7, height=9, units = 'in', dpi = 300)
    }
    else {
      ggsave(paste0('../article/pictures/20241107_reorganize/cv/', d_var, '_elastic_net.png'),
             width = 7, height=6, units = 'in', dpi = 300)
    }
    
  }
  
}


{
  # Plot the importance and estimate
  library(patchwork)
  dependent_vars <- c("VCI", "VSI", "FRI", "WMI", "PSI", "FSIQ")
  for (d_var in dependent_vars) {
    train_control <- trainControl(method = "repeatedcv",      # Cross-validation method
                                  number = 10,                # Number of folds
                                  repeats = 1,              # Number of repeats
                                  savePredictions = TRUE)
    d_var_p <- paste(d_var, 'p')
    d_var_est <- paste(d_var, 'estimate')
    d_var_low <- paste(d_var, 'lower')
    d_var_up <- paste(d_var, 'upper')
    sig_cpds <- batch2_data %>%
      filter(!!sym(d_var_p) < 0.05)
    formula <- as.formula(paste(d_var, "~", paste(sig_cpds$term, collapse = " + ")))
    cv_model <- train(formula,
                      data = merge_data, 
                      method = "glmnet", 
                      trControl = train_control)
    
    imp_data <- as.data.frame(varImp(cv_model)$importance) %>%
      mutate(Name = sig_cpds$Name)
    
    imp_data_fig <- imp_data %>%
      filter(Overall != 0) %>%
      arrange(desc(Overall)) %>%
      inner_join(batch2_data, by='Name') %>%
      mutate(Name = if_else(Name == 'Sorbic Acid', 'FA 6:2', Name),
             Name = if_else(Name == 'PC O-38:6 (isomer)', 'PC O-38:6', Name),
             Name = if_else(Name == 'PC O-36:5 (isomer)', 'PC O-36:5', Name),
             Name = if_else(Name == 'trans-Dodec-2-enoic acid', 'FA 12:1', Name),
             Name = if_else(Name == 'Perfluorohexanesulfonic acid', 'PFHxS', Name),
             Name = if_else(Name == '14-Deoxy-11,12-didehydroandrographolide', '14-dehydro Andrographolide', Name)) %>%
      pivot_longer(
        cols = starts_with("VCI") | starts_with("VSI") | starts_with("FRI") | starts_with("WMI") | starts_with("PSI") | starts_with("FSIQ"),
        names_to = c("Index", ".value"), 
        names_pattern = "(.*) (.*)"
      ) %>%
      mutate(Significance = ifelse(p < 0.05, "Significant", "Not Significant")) %>%
      mutate(Name = factor(Name, levels = unique(Name)))  %>%
      mutate(Index = factor(Index, levels = c("VCI", "VSI", "FRI", "WMI", "PSI", "FSIQ")))
    
    dtplot <- ggplot(imp_data_fig, aes(x = Name, y = estimate, shape = Index, color = Significance)) +
      geom_pointrange(aes(ymin = lower, ymax = upper), 
                      position = position_dodge(width = 0.5)) +  # Separate dependent variables to avoid overlap
      geom_point(size = 0.1, position = position_dodge(width = 0.5)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      scale_color_manual(values = c("Significant" = "#8856a7", "Not Significant" = "black")) +  # Red for p < 0.05
      scale_shape_manual(values = c(21, 22, 23, 24, 25, 8)) +  # Assign different shapes for different dependent variables
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
            )
    
    p_y <- ggplot(imp_data_fig, aes(x = Name, y = Overall)) +
      geom_segment(aes(xend = Name, y = 0, yend = Overall), color = "steelblue", linewidth=1) +
      geom_point(color = "red", size = 1) +  # Add a red dot at the end of each stick
      labs(
        y = "Importance",
        x = ""
      ) +
      scale_y_continuous(breaks = c(0, 50, 100), limits = c(0, 105)) +
      theme_classic() +
      theme(axis.text.y = element_text(color = 'black'),
            axis.text.x = element_text(color = "black", size=10, angle=45, hjust=1, vjust=1),
            axis.title.y = element_text(color = "black", size=9))
    
    # dtplot + p_y +
    #   plot_layout(heights = c(1, .2), guides='collect') &
    #   theme(legend.position = "bottom")
    dtplot + p_y +
        plot_layout(heights = c(1, .2), guides='collect') &
        theme(legend.position = "right")
    
    if (d_var == 'VCI') {
      ggsave(paste0('../article/pictures/20241107_reorganize/cv/', d_var, '_elastic_net_estimate.png'),
             width = 10, height=4, units = 'in', dpi = 300)
    }
    else {
      if (d_var == 'VSI') {
        dtplot <- ggplot(imp_data_fig, aes(x = Name, y = estimate, shape = Index, color = Significance)) +
          geom_pointrange(aes(ymin = lower, ymax = upper), 
                          position = position_dodge(width = 0.5)) +  # Separate dependent variables to avoid overlap
          geom_point(size = 0.1, position = position_dodge(width = 0.5)) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
          scale_color_manual(values = c("Significant" = "#8856a7", "Not Significant" = "black")) +  # Red for p < 0.05
          scale_shape_manual(values = c(21, 22, 23, 24, 25, 8)) +  # Assign different shapes for different dependent variables
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
          )
        
        p_y <- ggplot(imp_data_fig, aes(x = Name, y = Overall)) +
          geom_segment(aes(xend = Name, y = 0, yend = Overall), color = "steelblue", linewidth=1) +
          geom_point(color = "red", size = 1) +  # Add a red dot at the end of each stick
          labs(
            y = "Importance",
            x = ""
          ) +
          scale_y_continuous(breaks = c(0, 50, 100), limits = c(0, 105)) +
          theme_classic() +
          theme(axis.text.y = element_text(color = 'black'),
                axis.text.x = element_text(color = "black", size=8, angle=60, hjust=1, vjust=1),
                axis.title.y = element_text(color = "black", size=9))
        
        dtplot + p_y +
          plot_layout(heights = c(1, .2), guides='collect') &
          theme(legend.position = "right")
        ggsave(paste0('../article/pictures/20241107_reorganize/cv/', d_var, '_elastic_net_estimate.png'),
               width = 12, height=8, units = 'in', dpi = 300)
      }
      else {
        ggsave(paste0('../article/pictures/20241107_reorganize/cv/', d_var, '_elastic_net_estimate.png'),
             width = 12, height=8, units = 'in', dpi = 300)
      }
    }
  }
}