## Sezgi Ercan - 20.06.2025
## Script 6 for sum score analysis

##-------------------- Sum score ADE Model Fitting ----------------------##

# Necessary Libraries
library(tidyverse)
library(OpenMx)
library(MASS)
library(psych)
library(polycor)
library(gridExtra)
library(grid)
library(openxlsx)

rm(list = ls()) # clean memory


# This script takes "sumscore_Data.RDS" file from first script as input. !!! Not "sumscore_RawData.RDS" !!!!
# Please set the working directory to where the above-mentioned file is.
setwd("./sum_score_analyses/data_SA")


# Please indicate your desired output directory. ADE model outputs will be automatically saved there.
output_dir <- "./sum_score_analyses/ADE_model_Results_SA"


#
#-----------------------------Code Starts---------------------------#
#


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


scales <- c("emp", "dsm")
ages <- c(7,10,12)

# Empty df for the results
results_df <- data.frame()
trend_test <- data.frame()
corr_df <- data.frame()

scale_names <- list("emp" = "AP", "dsm" = "ADHP")


for (scale_i in scales){
  for (rated_age in ages){
    
    scale_name <- scale_names[[scale_i]]
    
    data <- df %>%
      filter(respondent == "m") %>%
      filter(assigned_age == rated_age) %>%
      filter(scale == scale_i)
    
    data <- data %>%
      group_by(twin_id) %>%
      filter(n() == 2) %>%
      ungroup()
    
    
    ## Converting to wide format
    # Separate shared and twin-specific variables
    shared_vars <- data %>%
      group_by(twin_id) %>%
      slice(1) %>%  # Take one row per twin_id
      ungroup() %>%
      dplyr::select(twin_id, zyg, yob, assigned_age)
    
    
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
    
    names(data_wide) <- c("twin_id", "zyg", "cohort", "assigned_age", "sex1", "sex2",
                          "cph1", "cph2")
    
    min_coh <- min(data_wide$cohort)
    
    # MZ and DZ datasets
    mz_data <- data_wide %>%
      filter(zyg == "MZ") %>%
      mutate(cohort = cohort-min_coh)
    
    dz_data <- data_wide %>%
      filter(zyg == "DZ") %>%
      mutate(cohort = cohort-min_coh)
    
    
    ## Continous to ordinal scale (0,1,2)
    # Cutting data
    mz_data$oph1 <- quantile_to_ordinal(mz_data$cph1)
    
    mz_data$oph2 <- quantile_to_ordinal(mz_data$cph2)
    
    dz_data$oph1 <- quantile_to_ordinal(dz_data$cph1)
    
    dz_data$oph2 <- quantile_to_ordinal(dz_data$cph2)
    
    
    #---------------------------------OpenMx--------------------------------------##
    
    mx_intervals = TRUE # TRUE if you want to have conf_intervarls in model results
    selVars <- c("oph1", "oph2")
    
    
    modeltype="ADE"  # the model type is ADE or ACE. here: ADE
    nv=1 #
    ntv=nv*2 #
    
    
    ## starting values
    # regression models for obtaining starting values
    
    mz1_coef <- summary(polr(as.factor(oph1) ~ cohort + sex1, data = mz_data, Hess = TRUE))$coeff
    mz2_coef <- summary(polr(as.factor(oph2) ~ cohort + sex2, data = mz_data, Hess = TRUE))$coeff
    
    dz1_coef <- summary(polr(as.factor(oph1) ~ cohort + sex1, data = dz_data, Hess = TRUE))$coeff
    dz2_coef <- summary(polr(as.factor(oph2) ~ cohort + sex2, data = dz_data, Hess = TRUE))$coeff
    
    
    # main effect of polynomial cohort on ph
    svpcoh1 <- round((mz1_coef[1] + mz2_coef[1] + dz1_coef[1] + dz2_coef[1]) / 4, 3)
    
    
    # B0
    svb0 = (mean(mz_data$cph1, na.rm = T) + mean(mz_data$cph2, na.rm = T) +
              mean(dz_data$cph1, na.rm = T) + mean(dz_data$cph2, na.rm = T)) / 4
    
    
    # main effect of sex on ph
    svb1_sex <- round((mz1_coef[4] + mz2_coef[4] + dz1_coef[4] + dz2_coef[4]) / 4, 3)
    
    
    # correlations and confidence intervals
    mz_cor <- polychor(mz_data$oph1, mz_data$oph2, ML=TRUE, std.err=TRUE)
    se_mz <- sqrt(mz_cor$var[1, 1])
    
    comp_cases_mz <- sum(complete.cases(mz_data[, c("oph1", "oph2")]))
    
    lower_mz <- round(mz_cor$rho - 1.96 * se_mz, 2)
    upper_mz <- round(mz_cor$rho + 1.96 * se_mz,2)
    
    cat(nrow(mz_data), comp_cases_mz, round(mz_cor$rho, 2), paste0("(",lower_mz, "-", upper_mz, ")"))
    
    
    dz_cor <- polychor(dz_data$oph1, dz_data$oph2, ML=TRUE, std.err=TRUE)
    se_dz <- sqrt(dz_cor$var[1, 1])
    
    comp_cases_dz <- sum(complete.cases(dz_data[, c("oph1", "oph2")]))
    
    lower_dz <- round(dz_cor$rho - 1.96 * se_dz, 2)
    upper_dz <- round(dz_cor$rho + 1.96 * se_dz,2)
    
    cat(nrow(dz_data), comp_cases_dz, round(dz_cor$rho, 2), paste0("(",lower_dz, "-", upper_dz, ")"))
    
    corr_df <- rbind(corr_df, data.frame(
      Scale = scale_name,
      Age = rated_age,
      MZ_total = nrow(mz_data),
      MZ_comp = comp_cases_mz,
      MZ_Correlation = mz_cor$rho,
      MZ_StdError = se_mz,
      MZ_CI_Lower = lower_mz,
      MZ_CI_Upper = upper_mz,
      DZ_total = nrow(dz_data),
      DZ_comp = comp_cases_dz,
      DZ_Correlation = dz_cor$rho,
      DZ_StdError = se_dz,
      DZ_CI_Lower = lower_dz,
      DZ_CI_Upper = upper_dz
      ))
    
    
    # unstandardized (raw) variances
    rmz = polychoric(mz_data[, c("oph1", "oph2")])$rho[2]
    rdz = polychoric(dz_data[, c("oph1", "oph2")])$rho[2]
    
    svVA = (4*rdz - rmz) * 1
    svVD = (2*rmz - 4*rdz) * 1
    svVE = (1 - svVA - svVD) * 1
    
    
    # starting values for the path coefficients
    svsdA = sqrt(svVA)
    svsdD = sqrt(svVD)
    svsdE = sqrt(svVE)
    
    # starting values for thresholds
    qmz = polychoric(mz_data[, c("oph1", "oph2")])$tau  # estimated thresholds
    qdz = polychoric(dz_data[, c("oph1", "oph2")])$tau  # estimated thresholds
    
    q2 = t(qmz +qdz)/2 # mz dz average
    thr1 = round(mean(q2[1,1:2]),2) # twin 1 twin 2 average
    thr2 = round(mean(q2[2,1:2]),2) #
    thfix = matrix(c(thr1,thr1,thr2,thr2),2,2,byrow=T) # 2x2 matrix for threshold starting values
    
    # Starting value for covariance matrix
    mz_cor <- polychoric(mz_data[, c("oph1", "oph2")])$rho[1,2]
    dz_cor <- polychoric(dz_data[, c("oph1", "oph2")])$rho[1,2]
    
    
    # ------------------- Openmx model starts --------------------- ##
    
    
    # Data objects
    mz_data$oph1 <- mxFactor(mz_data$oph1, levels=c(0:2) )
    mz_data$oph2 <- mxFactor(mz_data$oph2, levels=c(0:2) )
    dz_data$oph1 <- mxFactor(dz_data$oph1, levels=c(0:2) )
    dz_data$oph2 <- mxFactor(dz_data$oph2, levels=c(0:2) )
    #
    #
    # Create regression model for expected Mean Matrices
    B0_     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, 
                         values=svb0, labels=c("b0","b0"), name="B0" )
    
    Bcoh     	<- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE,
                           values=svpcoh1, label="bcoh", name="Bcoh" )
    
    Bsex     	<- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE,
                           values=svb1_sex, label="bsex", name="Bsex" )
    
    
    ThreshFixed = mxMatrix( type = "Full", nrow=2, ncol=2, free=F, 
                            labels = c('t0','t1','t0','t1'),values=thfix, name='Th') # the fixed thresholds based on the data
    #
    #
    # define sex and cohort
    Sex <- mxMatrix(type="Full", nrow=1, ncol=ntv, free=FALSE, 
                    labels=c('data.sex1','data.sex2'), name="sex")
    #
    Cohort    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=FALSE,
                             labels=c("data.cohort","data.cohort"), name="coh" )
    
    #
    # define cohorts for moderation effect
    Coh1 <- mxMatrix(type="Full", nrow=1, ncol=1, free=FALSE, labels=c('data.cohort'),name="coh1")
    
    #
    ExpMean <- mxAlgebra(expression = B0 + coh*Bcoh + sex*Bsex, name='expMean')
    
    #
    #
    #ACDE model
    # Create Matrices for Variance Components
    pathA      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svsdA, label="pA11", name="pthA" )
    #
    modA      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=0, label="mA11", name="mA" )
    
    #
    if (modeltype=="ACE") {
      pathC      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svsdC, label="pC11", name="pthC" )
      pathD      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=FALSE, values=0, label="pD11", name="pthD" )
      #
      modC      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=0, label="mC11", name="mC" )
      modD      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=FALSE, values=0, label="mD11", name="mD" )
      #
    }
    #
    if (modeltype=="ADE") {
      pathC      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=FALSE, values=0, label="pC11", name="pthC" )
      pathD      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svsdD, label="pD11", name="pthD" )
      modC      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=FALSE, values=0, label="mC11", name="mC" )
      modD      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=0, label="mD11", name="mD" )
    }
    #
    pathE      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svsdE, label="pE11", name="pthE" )
    modE      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=0, label="mE11", name="mE" )
    #
    # ACDE model
    # Create Matrices for Variance Components
    covE      <- mxAlgebra( expression = (pthE+mE%*%coh1)%*%(pthE+mE%*%coh1), name="VE" )
    covC      <- mxAlgebra( expression = (pthC+mC%*%coh1)%*%(pthC+mC%*%coh1), name="VC" )
    covD      <- mxAlgebra( expression = (pthD+mD%*%coh1)%*%(pthD+mD%*%coh1), name="VD" )
    covA      <- mxAlgebra( expression = (pthA+mA%*%coh1)%*%(pthA+mA%*%coh1), name="VA" )
    #ACDE model
    # Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
    covP      <- mxAlgebra( expression= VA+VC+VD+VE, name="V" )
    covMZ     <- mxAlgebra( expression= VA+VC+VD, name="cMZ" )
    covDZ     <- mxAlgebra( expression= 0.5%x%VA+VC+.25%x%VD, name="cDZ" )
    expCovMZ  <- mxAlgebra( expression= rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
    expCovDZ  <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )
    
    #ACE model
    # Create Data Objects for Multiple Groups
    #ACE model
    dataMZ    <- mxData( observed=mz_data, type="raw" )
    dataDZ    <- mxData( observed=dz_data, type="raw" )
    #
    #
    # Create Expectation Objects for Multiple Groups
    expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="expMean", threshold="Th", dimnames=selVars,threshnames=selVars )
    expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="expMean", threshold="Th", dimnames=selVars,threshnames=selVars )
    funML     <- mxFitFunctionML()
    #ACDE model
    # Create Model Objects for Multiple Groups
    pars      <- list(B0_, Bcoh, Bsex, ThreshFixed, pathA, pathC, pathD, pathE,
                      modA, modD, modC, modE)
    
    modelMZ   <- mxModel( pars, covA, covC, covD, covE, covP, ExpMean, Sex, Cohort, Coh1, covMZ, expCovMZ, ExpMean, dataMZ, expMZ, funML, name="MZ" )
    modelDZ   <- mxModel( pars, covA, covC, covD, covE, covP, ExpMean, Sex, Cohort, Coh1, covDZ, expCovDZ, ExpMean, dataDZ, expDZ, funML, name="DZ" )
    multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )
    # ACDE model
    # Create Confidence Interval Objects
    ciACDE     <- mxCI(c("mA","mD","mE"))
    # ACDE model
    # Build Model with Confidence Intervals
    omodelACDEmod  <- mxModel(name=modeltype, pars, modelMZ, modelDZ, multi, ciACDE )
    #ACDE model
    # Run ACDE Model
    ofitACDEmod2    <- mxTryHard( omodelACDEmod, intervals = T)
    osumACDEmod2    <- summary( ofitACDEmod2 )
    
    
    # test A mod
    otestA=omxSetParameters(ofitACDEmod2, labels=c('mA11'), value=0, free=FALSE, name = "DE")
    otestA_out    <- mxTryHard( otestA, intervals=F)
    test_A <- data.frame(mxCompare(ofitACDEmod2, otestA_out))
    
    # test D mod
    otestD=omxSetParameters(ofitACDEmod2, labels=c('mD11'), value=0, free=FALSE, name = "AE")
    otestD_out    <- mxTryHard( otestD, intervals=F)
    test_D <- data.frame(mxCompare(ofitACDEmod2, otestD_out))

    # test E mod
    otestE=omxSetParameters(ofitACDEmod2, labels=c('mE11'), value=0, free=FALSE, name = "AD")
    otestE_out    <- mxTryHard( otestE, intervals=F)
    test_E <- data.frame(mxCompare(ofitACDEmod2, otestE_out))
    
    # A and D dropped
    otestAD=omxSetParameters(ofitACDEmod2, labels=c("mA11", "mD11"), value=c(0,0), free=FALSE, name = "E")
    otestAD_out    <- mxTryHard( otestAD, intervals=F)
    test_AD <- data.frame(mxCompare(ofitACDEmod2, otestAD_out))
    
    # omnibus test
    om_test=omxSetParameters(ofitACDEmod2, labels=c("mA11", "mD11", "mE11"), value=c(0,0,0), free=FALSE, name = "OM")
    om_test_out    <- mxTryHard( om_test, intervals=F)
    test_OM <- data.frame(mxCompare(ofitACDEmod2, om_test_out))



    
    model_info <- data.frame("scale" = c("", scale_i),
                             "rater" = c("", "m"),
                             "age" = c("", rated_age)
    )
    
    test_res <- rbind(test_A, test_D, test_E, test_AD, test_OM)
    test_res <- cbind(test_res, model_info)
    
    trend_test <- rbind(trend_test, test_res)
    
    
    # Years and Birth cohorts for the change in variance
    combined_zyg <- rbind(mz_data, dz_data)
    years <- sort(unique(combined_zyg$cohort))
    years_df <- data.frame(years = years,
                           birth_cohort = years + min_coh,
                           birth_cohort_short = sprintf("%02d", (years + min_coh) %% 100))
    
    
    
    models <- list(ofitACDEmod2, otestA_out, otestD_out, otestE_out, otestAD_out, om_test_out)
    
    
    # Loop through the list
    for (each_model in models) {
      
      model_name <- each_model@name
      
      # Get the named vector of estimates
      out_paramters <- data.frame(each_model$output$estimate)
      
      # obtaining output parameters
      b0 = out_paramters["b0",1]
      bsex = out_paramters["bsex",1]
      bcoh = out_paramters["bcoh",1]
      
      path_A = out_paramters["pA11",1]
      path_D = out_paramters["pD11",1]
      path_E = out_paramters["pE11",1]
      
      m_A1 = ifelse(is.na(out_paramters["mA11",1]), 0, out_paramters["mA11",1])
      m_D1 = ifelse(is.na(out_paramters["mD11",1]), 0, out_paramters["mD11",1])
      m_E1 = ifelse(is.na(out_paramters["mE11",1]), 0, out_paramters["mE11",1])
      
      var_A = (path_A + years_df$years*m_A1)^2
      var_D = (path_D + years_df$years*m_D1)^2
      var_E = (path_E + years_df$years*m_E1)^2
      total_var = var_A + var_D + var_E
      
      var_df <- cbind(years_df, var_A, var_D, var_E, total_var)
      
      var_df <- var_df %>%
        mutate(std_var_A = var_A / total_var,
               std_var_D = var_D / total_var,
               std_var_E = var_E / total_var)
      
      B_mA = paste0(round(m_A1, 4), " (", round(ofitACDEmod2$output$confidenceIntervals["mA11",1], 4), 
                   " - ", round(ofitACDEmod2$output$confidenceIntervals["mA11",3], 4), ")")
      B_mD = paste0(round(m_D1, 4), " (", round(ofitACDEmod2$output$confidenceIntervals["mD11",1], 4), 
                   " - ", round(ofitACDEmod2$output$confidenceIntervals["mD11",3], 4), ")")
      B_mE = paste0(round(m_E1, 4), " (", round(ofitACDEmod2$output$confidenceIntervals["mE11",1], 4), 
                   " - ", round(ofitACDEmod2$output$confidenceIntervals["mE11",3],4), ")")
      
      
      # Correlation by birth cohorts
      mz_desc <- mz_data %>%
        group_by(cohort) %>%
        summarise(
          MZ_poly_cor = possibly(~ polychor(., oph2, ML=TRUE, std.err=TRUE)$rho, otherwise = NA)(oph1),
          n_MZ = n(),
          .groups = 'drop'
        )
      
      dz_desc <- dz_data %>%
        group_by(cohort) %>%
        summarise(
          Dz_poly_cor = possibly(~ polychor(., oph2, ML=TRUE, std.err=TRUE)$rho, otherwise = NA)(oph1),
          n_DZ = n(),
          .groups = 'drop'
        )
      
      res_df <- data.frame(
        "Scale" = scale_i,
        "Age" = rated_age,
        "Model" = model_name,
        "Birth cohort" = years_df$years + min_coh,
        "cohort" = years_df$years,
        "B_mA" = B_mA,
        "B_mD" = B_mD,
        "B_mE" = B_mE,
        "var_phen" = total_var,
        "varA" = var_A,
        "varD" = var_D,
        "varE" = var_E,
        "std_varA" = var_df$std_var_A,
        "std_varD" = var_df$std_var_D,
        "std_varE" = var_df$std_var_E
      )
      
      res_df <- res_df %>%
        left_join(mz_desc, by = "cohort") %>%
        left_join(dz_desc, by = "cohort")
      
      results_df <- rbind(results_df, res_df)

    }
    
    #-----------------------------------------------------------------
    
   
     # Print progress
    cat("Processed:", scale_i, rated_age, "\n")
    
  } 
  
}

write.xlsx(results_df, file = paste0(output_dir, "/ADE_parameters_SA.xlsx"))
write.xlsx(trend_test, file = paste0(output_dir, "/ADE_ModelResults_SA.xlsx"))
write.xlsx(corr_df, file = paste0(output_dir, "/polychor_SA.xlsx"))

  
