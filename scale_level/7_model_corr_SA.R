## Sezgi Ercan - 20.06.2025
## Script 7 for sum score analysis

##-------------------- Sum Score OpenMx Model for Trend ----------------------##
## This model is corrected for parental age

# Necessary Libraries
library(gridExtra)
library(tidyverse)
library(openxlsx)
library(OpenMx)
library(ggpubr)
library(scales)
library(psych)
library(MASS)
library(grid)

rm(list = ls()) # clean memory

# This script takes "sumscore_Data.RDS" file from first script as input. !!! Not "sumscore_RawData.RDS" !!!!
# Please set the working directory to where the above-mentioned file is.
setwd("./sum_score_analyses/data_SA")


# Please indicate your desired output directory. Outputs will be automatically saved there.
output_dir <- "./sum_score_analyses/model_results_SA/corr"


#
#-----------------------------Code Starts---------------------------#
#


# Quantile to ordinal function
quantile_to_ordinal <- function(x) {
  
  # Compute default quantiles
  probs <- c(0, 1/3, 2/3, 1)
  q <- quantile(x, probs = probs, na.rm = TRUE)
  
  # Check if breaks are unique
  if (length(unique(q)) < length(q)) {
    
    q <- c(-Inf, 0, median(x[x > 0], na.rm = TRUE), Inf)
    
  }
  
  # Return cut result as numeric
  as.numeric(as.character(cut(x, 
                              breaks = q, 
                              labels = FALSE, 
                              include.lowest = TRUE) - 1))
}


df <- readRDS("sumscore_Data.rds")
df <- df %>%
  mutate(sex_chr = sex,
         sex = ifelse(sex == "male", 0, 1))


df <- df %>%
  group_by(yob, respondent, assigned_age) %>%
  mutate(no_data = ifelse(all(is.na(score)), TRUE, FALSE)) %>%
  ungroup()

df <- df %>%
  filter(no_data == FALSE)


# loop specs
raters <- c("m", "f", "t")

ages <- c("7", "10", "12")

scales <- c("emp", "dsm")

sig_thresh <- 0.05/6


rater_names <- list("m" = "Mother", "f" = "Father", "t" = "Teacher")
scale_names <- list("emp" = "AP", "dsm" = "ADHP")


# Dataframes and plot list for the results
trend_plot_data <- data.frame()
openmx_results <- data.frame()
hist_list <- list()
hist_list_ord <- list()
ordinal_info <- data.frame()


for (scale_i in scales){

  scale_name <- scale_names[[scale_i]]

  for (rater in raters) {
      
    rater_name <- rater_names[[rater]]
      
    for (rated_age in ages) {
      
      # filter data
      data <- df %>%
        filter(scale == scale_i) %>%
        filter(respondent == rater) %>%
        filter(assigned_age == rated_age)
      
      # remove non-paired subjects
      data <- data %>%
        group_by(twin_id) %>%
        filter(n() == 2) %>%
        ungroup()
      
      
      hist_p <- ggplot(data, aes(x = score)) +
        geom_histogram(colour = "black", fill = "#1f77b4", binwidth = 1) +
        labs(
          x = "Sum Score",
          y = "Subjects (in Thousands)",
          title = paste0(rater_name, "-Rated ", rated_age, "-years-old Subjects ")
        ) +
        theme_classic() +
        theme(
          plot.margin = unit(c(0.5, 0.25, 0, 0.25), "cm"),
          plot.title = element_text(hjust = 0.5, size = 21, family = "Arial"),
          axis.title = element_text(size = 20, family = "Arial"),
          axis.text = element_text(size = 18, family = "Arial"),
          panel.grid.major.x = element_line(size = 0.5, color = "gray90"),
          panel.grid.major.y = element_line(size = 0.5, color = "gray90")
        ) +
        scale_y_continuous(
          labels = scales::label_number(scale = 1/1000, suffix = "K"))
      
      plot_name <- paste0(scale_i, "_", rater, rated_age)
      hist_list[[plot_name]] <- hist_p # Store the plot in the list
      
        
      ## Converting to wide format
      # Separate shared and twin-specific variables
      shared_vars <- data %>%
        group_by(twin_id) %>%
        slice(1) %>%  # Take one row per twin_id
        ungroup() %>%
        dplyr::select(twin_id, zyg, yob, assigned_age, parent_Age)
      
      
      # Reshape twin-specific variables
      twin_specific <- data %>%
        dplyr::select(twin_id, birthorder, sex, score) %>%
        pivot_wider(
          names_from = birthorder,
          values_from = c(sex, score),
          names_sep = "_twin"
        )
      
      # Merge the shared and twin-specific parts
      data_wide <- left_join(shared_vars, twin_specific, by = "twin_id") 
      
      data_wide <- as.data.frame(data_wide)

      
      names(data_wide) <- c("twin_id", "zyg", "cohort", "assigned_age", "parent_Age", "sex1", "sex2",
                            "cph1", "cph2")
      
      min_coh <- min(data_wide$cohort)
      min_parent_age <- min(data_wide$parent_Age, na.rm = TRUE)
      
      # MZ and DZ datasets  
      mz_data <- data_wide %>%
        filter(zyg == "MZ")
      dz_data <- data_wide %>%
        filter(zyg == "DZ")
      
      
      ## Continous to ordinal scale (0,1,2)
      # Cutting data
      mz_data$oph1 <- quantile_to_ordinal(mz_data$cph1)
      
      mz_data$oph2 <- quantile_to_ordinal(mz_data$cph2)
      
      dz_data$oph1 <- quantile_to_ordinal(dz_data$cph1)
      
      dz_data$oph2 <- quantile_to_ordinal(dz_data$cph2)
      
      
      quantiles <- data.frame("scale" = scale_i,
                              "rater" = rater,
                              "age" = rated_age,
                              "mz1" = table(mz_data$oph1),
                              "mz2" = table(mz_data$oph2),
                              "dz1" = table(dz_data$oph1),
                              "dz2" = table(dz_data$oph2))
      
      
      ordinal_info <- rbind(ordinal_info, quantiles)
      
      
      # Ordinal Distribution Plot
      ord_hist_df <- data.frame(c(mz_data$oph1, mz_data$oph2, dz_data$oph1, dz_data$oph2))
      names(ord_hist_df) <- c("sumscores_ordinal")
      ord_hist_df <- na.omit(ord_hist_df)
      
      hist_p_ord <- ggplot(ord_hist_df, aes(x = as.factor(sumscores_ordinal))) +
        geom_bar(colour = "black", fill = "#1f77b4") +
        labs(
          x = "Ordinal Sum Scores",
          y = "Number of Subjects",
          title = paste0(rater_name, "-Rated ", rated_age, "-years-old Subjects ")
        ) +
        theme_classic() +
        theme(
          plot.margin = unit(c(0.5, 0.25, 0, 0.25), "cm"),
          plot.title = element_text(hjust = 0.5, size = 18, family = "Arial"),
          axis.title = element_text(size = 18, family = "Arial"),
          axis.text = element_text(size = 15, family = "Arial"),
          panel.grid.major.x = element_line(size = 0.5, color = "gray90"),
          panel.grid.major.y = element_line(size = 0.5, color = "gray90")
        )
      
      hist_p_ord
      hist_list_ord[[plot_name]] <- hist_p_ord # Store the plot in the list
      
      
      # standardize cohort
      mz_data <- mz_data %>%
        mutate(cohort = cohort - min_coh) %>%
        mutate(parent_Age = parent_Age - min_parent_age)
      
      dz_data <- dz_data %>%
        mutate(cohort = cohort - min_coh) %>%
        mutate(parent_Age = parent_Age - min_parent_age)
      
      
      ddatmz <- mz_data %>%
        dplyr::select(c(oph1, oph2, cohort, sex1, sex2, parent_Age))
      
      ddatdz <- dz_data %>%
        dplyr::select(c(oph1, oph2, cohort, sex1, sex2, parent_Age))
      
      
      
      ##---------------------------------OpenMx-------------------------------------##
      mx_intervals = TRUE # TRUE if you want to have conf_intervarls in model results
      
      
      ## starting values
      # regression models for obtaining starting values
      
      mz1_coef <- summary(polr(as.factor(oph1) ~ cohort + sex1 + parent_Age, data = ddatmz, Hess = TRUE))$coeff
      mz2_coef <- summary(polr(as.factor(oph2) ~ cohort + sex2 + parent_Age, data = ddatmz, Hess = TRUE))$coeff
      
      dz1_coef <- summary(polr(as.factor(oph1) ~ cohort + sex1 + parent_Age, data = ddatdz, Hess = TRUE))$coeff
      dz2_coef <- summary(polr(as.factor(oph2) ~ cohort + sex2 + parent_Age, data = ddatdz, Hess = TRUE))$coeff
      
      
      # main effect of polynomial cohort on ph
      svpcoh1 <- round((mz1_coef[1] + mz2_coef[1] + dz1_coef[1] + dz2_coef[1]) / 4, 3)
      
      # main effect of sex on ph
      svb1_sex <- round((mz1_coef[2] + mz2_coef[2] + dz1_coef[2] + dz2_coef[2]) / 4, 3)
      
      # main effect of parental Age on ph
      svb1_Age <- round((mz1_coef[3] + mz2_coef[3] + dz1_coef[3] + dz2_coef[3]) / 4, 3)
      
      
      # starting values for thresholds
      qmz = polychoric(ddatmz[, c("oph1", "oph2")])$tau  # estimated thresholds
      qdz = polychoric(ddatdz[, c("oph1", "oph2")])$tau  # estimated thresholds
      
      q2 = t(qmz +qdz)/2 # mz dz average
      thr1 = round(mean(q2[1,1:2]),2) # twin 1 twin 2 average
      thr2 = round(mean(q2[2,1:2]),2) #
      svthr = matrix(c(thr1,thr1,thr2,thr2),2,2,byrow=T) # 2x2 matrix for threshold starting values
      
      # Starting value for covariance matrix
      mz_cor <- polychoric(ddatmz[, c("oph1", "oph2")])$rho[1,2]
      dz_cor <- polychoric(ddatdz[, c("oph1", "oph2")])$rho[1,2]
      
      
      #-------------------------------Model Starts-----------------------------------#
      #-------------------------------Model Starts-----------------------------------#
      # remove the rows with NA values for the definition variables
      cols <- c("parent_Age","sex1","sex2","cohort")
      keepMZ <- complete.cases(ddatmz[ , cols])
      keepDZ <- complete.cases(ddatdz[ , cols])
      
      ddatmz <- ddatmz[keepMZ, , drop=FALSE]
      ddatdz <- ddatdz[keepDZ, , drop=FALSE]
      
      
      selVars=c('oph1', 'oph2')
      #
      # SAT Model
      # Data objects for Multiple Groups .. define as ordinal
      #
      ddatmz$oph1 <- mxFactor(ddatmz$oph1, levels=c(0:2) )
      ddatmz$oph2 <- mxFactor(ddatmz$oph2, levels=c(0:2) )
      ddatdz$oph1 <- mxFactor(ddatdz$oph1, levels=c(0:2) )
      ddatdz$oph2 <- mxFactor(ddatdz$oph2, levels=c(0:2) )
      #
      dataMZ   	<- mxData( observed=ddatmz, type="raw" )
      dataDZ   	<- mxData( observed=ddatdz, type="raw" )  
        
      #
      covMZ     	<- mxMatrix(type='Symm',nrow=2, ncol=2, free=T, values=c(1,mz_cor,1), labels=c('vphmz','cmz','vphmz'), name='SMZ')
      covDZ    	<- mxMatrix(type='Symm',nrow=2, ncol=2, free=T,  values=c(1,dz_cor,1), labels=c('vphdz','cdz','vphdz'), name='SDZ') 
      
      
      # Matrix & Algebra for expected means vector and expected thresholds
      b0    	<- mxMatrix( type="Full", nrow=1, ncol=2, labels=c('b0','b0'),free=TRUE, name="B0" )
      
      threshold 	<- mxMatrix( type="Full", nrow=2, ncol=2, free=FALSE, 
                              values = svthr, name="threshold" )
      
      # Covariate matrices
      defct    	<- mxMatrix( type="Full", nrow=1, ncol=2, free=FALSE, 
                             labels=c("data.cohort","data.cohort"), name="Coh" )

      defsex    	<- mxMatrix( type="Full", nrow=1, ncol=2, free=FALSE, 
                              labels=c("data.sex1","data.sex2"), name="Sex" )
      
      defAge    	<- mxMatrix( type="Full", nrow=1, ncol=2, free=FALSE, 
                              labels=c("data.parent_Age","data.parent_Age"), name="Age" )
      
      
      # Betas
      Bsex     	<- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, 
                             values=svb1_sex, label="bsex", name="Bsex" )
      
      B_Age     	<- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, 
                              values=svb1_Age, label="b_Age", name="B_Age" )
      
      Bcoh     	<- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, 
                             values=svpcoh1, label="bcoh", name="Bcoh" )
      
      # 
      expm <- mxAlgebra(expression = B0 + Coh*Bcoh + Sex*Bsex + Age*B_Age, name='Mean') 
      
      #
      # Create Algebra for expected Threshold Matrices
      #
      # Algebra to compute total variances and standard deviations (diagonal only)
      #
      # Algebra for expected Mean and Variance/Covariance Matrices in MZ & DZ twins
      covMZ_    	<- mxAlgebra( expression= SMZ, name="expCovMZ" )
      covDZ_    	<- mxAlgebra( expression= SDZ, name="expCovDZ" )
      #
      ##
      # Expectation objects for Multiple Groups
      expMZ     	<- 	mxExpectationNormal( covariance="expCovMZ", means="Mean", dimnames=selVars,
                                          thresholds="threshold" )
      expDZ     	<- 	mxExpectationNormal( covariance="expCovDZ", means="Mean", dimnames=selVars,
                                          thresholds="threshold" )
      #
      #
      bits 		<-c(expm, b0,defct, defAge, defsex, Bsex, Bcoh, B_Age)
      funML     	<- 	mxFitFunctionML()
      #
      modelMZ  	<- mxModel( bits,threshold, covMZ, covMZ_, dataMZ, funML, expMZ, name="MZ" )
      modelDZ  	<- mxModel( bits,threshold, covDZ, covDZ_, dataDZ, funML,expDZ, name="DZ" )
      #
      #
      #
      # Create Confidence Interval Objects
      ciCoh     <- mxCI("bcoh", interval = 1 - sig_thresh)
      #
      #
      # Combine Groups
      multi     	<- mxFitFunctionMultigroup( c("MZ","DZ") )
      satModel 	<- mxModel( "SAT", modelMZ, modelDZ, funML, multi,ciCoh  )
      #
      # ---
      out1 = mxTryHard(satModel, intervals = T)
      
      # test of no trend
      no_trend <- mxModel(out1, name="no_trend" )
      no_trend <- omxSetParameters(no_trend, labels='bcoh', values=0, free=F)
      fit_no_trend <- mxTryHard(no_trend, intervals = mx_intervals)
      
      trend_test <- data.frame(mxCompare(out1, fit_no_trend))   # ns
      
      # --- RESCUE BLOCK FOR NEGATIVE diffLL --- #
      if (!is.na(trend_test[2, "diffLL"]) && trend_test[2, "diffLL"] < -0.1) {
        
        cat("\nNegative diffLL detected! Rescuing SAT model with nested parameters...\n")
        
        # 1. Extract the good parameters from the nested model
        good_params <- omxGetParameters(fit_no_trend)
        
        # 2. Inject them into the SAT model as starting values
        sat_rescue <- omxSetParameters(satModel, 
                                       labels = names(good_params), 
                                       values = good_params)
        
        # 3. Secure the global minimum FIRST 
        sat_rescued_fit <- mxTryHard(sat_rescue, intervals = FALSE, extraTries = 50)
        
        # Now compute the Confidence Intervals
        out1 <- mxRun(sat_rescued_fit, intervals = TRUE)
        
        # Re-run the nested model
        no_trend <- mxModel(out1, name="no_trend" )
        no_trend <- omxSetParameters(no_trend, labels='bcoh', values=0, free=F)
        fit_no_trend <- mxTryHard(no_trend, intervals = mx_intervals)
        
        # Overwrite the trend_test with the corrected comparison
        trend_test <- data.frame(mxCompare(out1, fit_no_trend))
      }
      # ---------------------------------------- #
      
      model_info <- data.frame("scale" = c("", scale_i),
                               "rater" = c("", rater),
                               "age" = c("", rated_age)
      )
      
      trend_test <- cbind(trend_test, model_info)
      print(c(rater_name, scale_name, rated_age))
      print(trend_test)
      
      ##--------------------------------End of OpenMx Model-------------------------##
      
      # Adding all the results in to a df
      openmx_results <- rbind(openmx_results, trend_test)
      
      
      # dataframe for plotting the trend
      out_paramters <- data.frame(out1$output$estimate)
      
      # obtaining outbput parameters
      b0 = out_paramters["b0",1]
      bsex = out_paramters["bsex",1]
      bcoh = out_paramters["bcoh",1]
      bAge = out_paramters["b_Age",1]
      vphmz = out_paramters["vphmz",1]
      vphdz = out_paramters["vphdz",1]
      thresh_1 = svthr[1,1]
      thresh_2 = svthr[2,2]
      l_conf_int = tryCatch(round(out1$output$confidenceIntervals["bcoh",1], 4), error = function(e) NA)
      u_conf_int = tryCatch(round(out1$output$confidenceIntervals["bcoh",3], 4), error = function(e) NA)                    
      
      # If CI fails, fallback to CI (Standard Error based)
      if (is.na(l_conf_int) || is.na(u_conf_int)) {
        cat(" -> Profile CI missing: Calculating Standard Error based CI instead...")
        sum_out <- summary(out1)
        b_se <- sum_out$parameters[sum_out$parameters$name == "bcoh", "Std.Error"]
        
        # Calculate Z-value based onsig_thresh (0.05/6)
        z_val <- qnorm(1 - (sig_thresh / 2))
        
        l_conf_int <- round(bcoh - (z_val * b_se), 4)
        u_conf_int <- round(bcoh + (z_val * b_se), 4)
      }                 
      
      
      p_val = trend_test[2, "p"]
      
      
      significance <- ifelse(p_val < sig_thresh, TRUE, FALSE)
      
      # Creating years dataset for plotting the trend
      combined_zyg <- rbind(mz_data, dz_data)
      
      years_df <- data.frame(unique(combined_zyg$cohort))
      colnames(years_df) <- "cohort"
      years_df$cohort <- years_df[order(years_df$cohort), ]
      
      
      
      trend_plot_df <- years_df %>%
        mutate(
          mean_liability = (bcoh * cohort),
          beta_0 = b0,
          beta_bc = bcoh,
          beta_sex = bsex,
          beta_Age = bAge
          
        ) %>%
        mutate(est_prob_0 = pnorm(thresh_1, mean = mean_liability, sd = sqrt((vphmz + vphdz)/2), lower.tail = T, log.p = FALSE)) %>%
        mutate(est_prob_2 = pnorm(thresh_2, mean = mean_liability, sd = sqrt((vphmz + vphdz)/2), lower.tail = F, log.p = FALSE)) %>%
        mutate(est_prob_1 = 1- (est_prob_0 + est_prob_2))
      
      # Observed probabilities-------------------------
      
      prob_lookUp <- combined_zyg %>%
        pivot_longer(cols = c("oph1", "oph2"), 
                     names_to = "birthorder", 
                     values_to = "ph") %>%
        filter(!is.na(ph))
      
      
      observed_probs <- prob_lookUp %>%
        group_by(cohort) %>%
        summarise(count_0 = sum(ph == 0, na.rm = T),
                  count_1 = sum(ph == 1, na.rm = T),
                  count_2 = sum(ph == 2, na.rm = T),
                  n = n()
        ) %>%
        ungroup() %>%
        mutate(obs_prob_0 = count_0/n,
               obs_prob_1 = count_1/n,
               obs_prob_2 = count_2/n)
      
      # Observed probabilities-------------------------
      
      trend_plot_df <- trend_plot_df %>%
        left_join(observed_probs, by = c("cohort"))
      
      
      trend_plot_df <- trend_plot_df %>%
        mutate(lower_thresh = thresh_1,
               upper_thresh = thresh_2,
               up_conf = u_conf_int,
               low_conf = l_conf_int,
               scale = scale_i,
               rater = rater,
               rated_age = rated_age,
               significance = significance) %>%
        dplyr::select(-count_0, -count_1, -count_2)
      
      trend_plot_df <- trend_plot_df %>%
        mutate(coded_cohort = cohort, 
          cohort = cohort + min_coh)
      
      
      trend_plot_data <- rbind(trend_plot_data, trend_plot_df)

      # printing the progress
      cat("Processed:", scale_i, rater, rated_age, "\n")
  
  
    } # age loop
  } # rater loop    
} # scale loop


# Format output
trend_plot_data_formatted <- trend_plot_data %>%
  mutate(
    `99.17% CI for b_bc` = paste0("[", low_conf, ", ", up_conf, "]")
    # Rename scale and rater variables
    # scale = case_when(
    #   scale == "emp" ~ "AP",
    #   scale == "dsm" ~ "ADHP",
    #   TRUE ~ scale
    # ),
    # rater = case_when(
    #   rater == "m" ~ "Mother",
    #   rater == "f" ~ "Father",
    #   rater == "t" ~ "Teacher",
    #   TRUE ~ rater
    # )
  ) %>%
  dplyr::select(
    scale = scale,
    rater = rater,
    rated_age = rated_age,
    `Birth Cohort` = cohort,
    `Birth Cohort as Coded` = coded_cohort,
    `N (sample size)` = n,
    `lower_thresh` = lower_thresh,
    `upper_thresh` = upper_thresh,
    `Mean Liability` = mean_liability,
    `b_0` = beta_0,
    `beta_bc` = beta_bc,
    `beta_sex` = beta_sex,
    `beta_Age` = beta_Age,
    `99.17% CI for b_bc`,
    `significance` = significance,
    `Estimated Prob(0)` = est_prob_0,
    `Estimated Prob(1)` = est_prob_1,
    `Estimated Prob(2)` = est_prob_2,
    `Observed Prob(0)` = obs_prob_0,
    `Observed Prob(1)` = obs_prob_1,
    `Observed Prob(2)` = obs_prob_2,
    low_conf,
    up_conf
  )


# Save excel results
file_path_plotDf <- file.path(output_dir, "openmx_parameters_corr_SA.xlsx")
write.xlsx(trend_plot_data_formatted, file = file_path_plotDf, rowNames = FALSE)

file_path_res <- file.path(output_dir, "openmx_results_corr_SA.xlsx")
write.xlsx(openmx_results, file = file_path_res, rowNames = FALSE)

file_path_res <- file.path(output_dir, "ordinal_info_corr_SA.xlsx")
write.xlsx(ordinal_info, file = file_path_res, rowNames = FALSE)


# Saving the histograms
emp_list <- hist_list[grepl("^emp", names(hist_list))]
dsm_list <- hist_list[grepl("^dsm", names(hist_list))]


emp_plot <- grid.arrange(grobs = emp_list, ncol = 3, nrow = 3)


plot_path_emp <- file.path(output_dir, paste0("emp", "_hist_corr.png"))
ggsave(plot_path_emp, plot = emp_plot, width = 1920 / 100, height = 1440 / 100, dpi = 300, units = "in", bg = "white")


dsm_plot <- grid.arrange(grobs = dsm_list, ncol = 3, nrow = 3)


plot_path_dsm <- file.path(output_dir, paste0("dsm", "_hist_corr.png"))
ggsave(plot_path_dsm, plot = dsm_plot, width = 1920 / 100, height = 1440 / 100, dpi = 300, units = "in", bg = "white")


# Ordinal histograms
# Saving the histograms
emp_list_ord <- hist_list_ord[grepl("^emp", names(hist_list_ord))]
dsm_list_ord <- hist_list_ord[grepl("^dsm", names(hist_list_ord))]


emp_plot_ord <- grid.arrange(grobs = emp_list_ord, ncol = 3, nrow = 3,
                         top = textGrob("Distribution of the Ordinal Sum Scores - AP Scale",
                                        gp = gpar(fontsize = 20)))


plot_path_emp_ord <- file.path(output_dir, paste0("emp_ord", "_hist_corr.png"))
ggsave(plot_path_emp_ord, plot = emp_plot_ord, width = 1920 / 100, height = 1440 / 100, dpi = 300, units = "in", bg = "white")


dsm_plot_ord <- grid.arrange(grobs = dsm_list_ord, ncol = 3, nrow = 3,
                         top = textGrob("Distribution of the Ordinal Sum Scores - ADHP Scale",
                                        gp = gpar(fontsize = 20)))


plot_path_dsm_ord <- file.path(output_dir, paste0("dsm_ord", "_hist_corr.png"))
ggsave(plot_path_dsm_ord, plot = dsm_plot_ord, width = 1920 / 100, height = 1440 / 100, dpi = 300, units = "in", bg = "white")




