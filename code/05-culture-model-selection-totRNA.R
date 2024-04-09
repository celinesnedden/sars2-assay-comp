# This runs the full model selection procedure for the culture models,
#     using totRNA as the primary predictor

# Prep environment -------------------------------------------------------------

## Install CmdStanR, CmdStan, and RStan for model compilation ------------------

#  We recommend you review the installation instructions at:
#     https://mc-stan.org/cmdstanr/articles/cmdstanr.html 
#     https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started 

library("cmdstanr")
library("rstan")


## Install other packages ------------------------------------------------------

req_pkgs <- c("loo", "stringr", "tidyverse", "gtools", "ggpubr", "matrixStats")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependences = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))


## Small convenience function --------------------------------------------------

`%notin%` <- Negate(`%in%`)  # for convenience


## Set chains, iterations, cores, & seed ----------------------------------------

options(mc.cores = parallel::detectCores())

debug <- FALSE
if (debug == FALSE) {
  n_chains <- 6
  n_iter <- 4000
} else if (debug == TRUE) {
  n_chains <- 1
  n_iter <- 500
}

seed <- 5923


# Load data --------------------------------------------------------------------

dat <- read.csv("./data/pred-sg-data.csv")


# Prep data --------------------------------------------------------------------

dat.full <- subset(dat, !is.na(pos_inf) & cens_total %in% c(0, 1)) # for calculations below
dat.culture <- subset(dat, cens_total == 0 & !is.na(pos_inf)) # for fitting

# Convert columns to numerics
dat.culture$val_sg_pred <- as.numeric(dat.culture$val_sg_pred)
dat.culture$val_total <- as.numeric(dat.culture$val_total)

# Set articles as numeric factors
dat.culture$article_idx <-  as.numeric(factor(dat.culture$article))

# Set TG indicator based on total RNA protocol: ordered by decreasing expected RNA quantity
dat.culture$tg_idx[dat.culture$tg_total %in% c("N")] <- 1
dat.culture$tg_idx[dat.culture$tg_total %in% c("E")] <- 2
dat.culture$tg_idx[dat.culture$tg_total %in% c("S")] <- 3


# Evaluate RNA negative samples ------------------------------------------------

dat.neg <- subset(dat, cens_total == 1 & !is.na(pos_inf))

# total RNA 
TN_T_neg <- nrow(subset(dat.neg, pos_inf == 0 & cens_total == 1))
FN_T_neg <- nrow(subset(dat.neg, pos_inf == 1 & cens_total == 1))



# Run model --------------------------------------------------------------------

## Load model ------------------------------------------------------------------

model.culture.xval <- cmdstan_model('./code/stan-model-culture-xval.stan')


## Assign folds ----------------------------------------------------------------

n_folds <- 10
set.seed(16669)
dat.culture$fold <- kfold_split_stratified(K = n_folds, x = dat.culture$article_idx)


## Prep results table ----------------------------------------------------------

# probably read in table from culture base models code ?
tbl_cols <- c("model", "predictors", "elpd", "elpd_se", "elpd_loo", "elpd_loo_se", 
              "p_loo_loo", "p_loo_loo_se", "pareto_k_max",
              "F_score", "MCC",
              "n_correct_test", "n_correct_train",
              "percent_total", "percent_pos", "percent_neg")

tblS4 <- as.data.frame(matrix(nrow = 0, ncol = length(tbl_cols)))
colnames(tblS4) <- tbl_cols
log_lik <- as.data.frame(matrix(nrow = n_chains * n_iter / 2, 
                                ncol = nrow(dat.culture)))


## Run simplest model ----------------------------------------------------------

pred.df <- list(true_val = c(),
                pred_val = c())


### K-fold cross validation ---------------------------------------------------

n_corr_test <- 0
n_corr_train <- 0

for (k in 1:n_folds) {
  
  cat("Running fold ", k, "for the quantitative base model: T \n \n")
  
  # For ease of formatting data for stan
  dat.train <- subset(dat.culture, fold != k)
  dat.test <- subset(dat.culture, fold == k)
  
  # Formatting for stan
  dat.stan <- c(list(
    
    # Number of observations, articles, & levels 
    N_test = dim(dat.test)[1],
    N_train = dim(dat.train)[1],
    N_article = max(dat.train$article_idx),
    L_sp = max(dat.train$sp_idx),
    L_st = max(dat.train$st_idx),
    L_tg = max(dat.train$tg_idx),
    L_dpi = max(dat.train$dpi_idx),
    L_age = max(dat.train$age_idx),
    L_cell = max(dat.train$cell_idx),
    
    # Toggle predictors on/off
    t_inc = 1,
    sg_inc = 0,
    dose_inc = 0,
    dpi_inc = 0,
    sp_inc = 0,
    tg_inc = 0,
    st_inc = 0,
    age_inc = 0,
    sex_inc = 0,
    cell_inc = 0,
    assay_inc = 0,
    
    # Response variables
    pos_inf_train = dat.train$pos_inf,
    pos_inf_test = dat.test$pos_inf,
    
    # Predictor variables
    t_train = dat.train$val_total,
    t_test = dat.test$val_total,
    sg_train = dat.train$val_sg_pred,
    sg_test = dat.test$val_sg_pred,
    dose_train = dat.train$log10_dose_pfu,
    dose_test = dat.test$log10_dose_pfu,
    tg_train = dat.train$tg_idx,
    tg_test = dat.test$tg_idx,
    dpi_train = dat.train$dpi_idx,
    dpi_test = dat.test$dpi_idx,
    age_train = dat.train$age_idx,
    age_test = dat.test$age_idx,
    sex_train = dat.train$sex_idx,
    sex_test = dat.test$sex_idx,
    st_train = dat.train$st_idx,
    st_test = dat.test$st_idx,
    sp_train = dat.train$sp_idx,
    sp_test = dat.test$sp_idx,
    cell_train = dat.train$cell_idx,
    cell_test = dat.test$cell_idx,
    assay_train = dat.train$assay_idx,
    assay_test = dat.test$assay_idx,
    
    # Article hierarchy
    article_train = dat.train$article_idx, 
    article_test = dat.test$article_idx))
  
  # Run current model
  fit.current <- model.culture.xval$sample(data = dat.stan,
                                           chains = n_chains,
                                           parallel_chains = n_chains,
                                           iter_warmup = n_iter / 2,
                                           iter_sampling = n_iter / 2,
                                           refresh = 1000,
                                           seed = seed)
  
  # Store log likelihood for all test data
  log_lik_pars <- paste0("log_lik[", 1:nrow(dat.test), "]")
  log_lik[, dat.culture$fold == k] <- as.data.frame(fit.current$draws(log_lik_pars, 
                                                                    format = "matrix"))
  
  # Determine number of correctly classified samples in test data
  p_pos_test <- paste0("p_pos_test[", 1:nrow(dat.test), "]") 
  p_test_df <- fit.current$draws(p_pos_test, format = "df")
  p_test_mean <- as.data.frame(colMeans(p_test_df))
  p_test_mean <- p_test_mean[rownames(p_test_mean) %notin% c(".chain",
                                                             ".iteration",
                                                             ".draw"), ]
  pred_test <- p_test_mean
  pred_test[p_test_mean >= 0.5] <- 1
  pred_test[p_test_mean < 0.5] <- 0
  n_correct <- sum(dat.test$pos_inf == pred_test)
  n_corr_test <- n_corr_test + n_correct
  
  # Add to prediction df
  pred.df$true_val <- c(pred.df$true_val, dat.test$pos_inf)
  pred.df$pred_val <- c(pred.df$pred_val, pred_test)
  
  # Determine number of correctly classified samples in training data
  p_pos_train <- paste0("p_pos_train[", 1:nrow(dat.train), "]")
  p_train_df <- fit.current$draws(p_pos_train, format = "df")
  p_train_mean <- as.data.frame(colMeans(p_train_df))
  p_train_mean <- p_train_mean[rownames(p_train_mean) %notin% c(".chain",
                                                                ".iteration",
                                                                ".draw"), ]
  pred_train <- p_train_mean
  pred_train[p_train_mean >= 0.5] <- 1
  pred_train[p_train_mean < 0.5] <- 0
  n_correct <- sum(dat.train$pos_inf == pred_train)
  n_corr_train <- n_corr_train + n_correct
  
}

# Generate statistics based on equations
elpd_current <- sum(log(colMeans(exp(as.matrix(log_lik)))))
elpd_se_current <- sqrt(ncol(log_lik)) * sd(log(colMeans(exp(as.matrix(log_lik)))))
looic_current <- -2 * elpd_current 
looic_se_current <- 2 * elpd_se_current
ploo_current <- sum(colVars(as.matrix(log_lik)))
ploo_se_current <- sqrt(ncol(log_lik)) * sd(colVars(as.matrix(log_lik)))
pareto_k_max <- NA

# Generate F-score & MCC
TP <- sum(pred.df$true_val == 1 & pred.df$pred_val == 1)
FP <- sum(pred.df$true_val == 0 & pred.df$pred_val == 1)
TN <- sum(pred.df$true_val == 0 & pred.df$pred_val == 0) + TN_T_neg
FN <- sum(pred.df$true_val == 1 & pred.df$pred_val == 0) + FN_T_neg
F_score <- (2 * TP) / (2 * TP + FP + FN)
MCC <- ((TP * TN) - (FP * FN))/(sqrt(TP + FP)*sqrt(TP + FN)*sqrt(TN + FP)*sqrt(TN + FN))

# Percent of samples correctly classified, stratified by Pos/Neg/All 
percent_pos <- TP / (TP + FN)  * 100
percent_neg <- TN / (TN + FP) * 100
percent_total <- (TN + TP) / (TP + FP + TN + FN) * 100

# Store log_lik data frame for future reference
loglik_name <- paste0("log_lik_", "c1") # model name for future ref
assign(loglik_name, log_lik)


### PSIS-LOO method ------------------------------------------------------------

# Will run slightly slower, but prevents needing an additional model file
dat.train <- dat.culture
dat.test <- dat.culture

# Formatting for stan
dat.stan <- c(list(
  
  # Number of observations, articles, & levels 
  N_test = dim(dat.test)[1],
  N_train = dim(dat.train)[1],
  N_article = max(dat.train$article_idx),
  L_sp = max(dat.train$sp_idx),
  L_st = max(dat.train$st_idx),
  L_tg = max(dat.train$tg_idx),
  L_dpi = max(dat.train$dpi_idx),
  L_age = max(dat.train$age_idx),
  L_cell = max(dat.train$cell_idx),
  
  # Toggle predictors on/off
  t_inc = 1,
  sg_inc = 0, 
  dose_inc = 0,
  dpi_inc = 0,
  sp_inc = 0,
  tg_inc = 0,
  st_inc = 0,
  age_inc = 0,
  sex_inc = 0,
  cell_inc = 0,
  assay_inc = 0,
  
  # Response variables
  pos_inf_train = dat.train$pos_inf,
  pos_inf_test = dat.test$pos_inf,
  
  # Predictor variables
  t_train = dat.train$val_total,
  t_test = dat.test$val_total,
  sg_train = dat.train$val_sg_pred,
  sg_test = dat.test$val_sg_pred,
  dose_train = dat.train$log10_dose_pfu,
  dose_test = dat.test$log10_dose_pfu,
  tg_train = dat.train$tg_idx,
  tg_test = dat.test$tg_idx,
  dpi_train = dat.train$dpi_idx,
  dpi_test = dat.test$dpi_idx,
  age_train = dat.train$age_idx,
  age_test = dat.test$age_idx,
  sex_train = dat.train$sex_idx,
  sex_test = dat.test$sex_idx,
  st_train = dat.train$st_idx,
  st_test = dat.test$st_idx,
  sp_train = dat.train$sp_idx,
  sp_test = dat.test$sp_idx,
  cell_train = dat.train$cell_idx,
  cell_test = dat.test$cell_idx,
  assay_train = dat.train$assay_idx,
  assay_test = dat.test$assay_idx,
  
  # Article hierarchy
  article_train = dat.train$article_idx, 
  article_test = dat.test$article_idx))

# Run current model
fit.current <- model.culture.xval$sample(data = dat.stan,
                                         chains = n_chains,
                                         parallel_chains = n_chains,
                                         iter_warmup = n_iter / 2,
                                         iter_sampling = n_iter / 2,
                                         refresh = 1000)

# Store log likelihood
log_lik_pars <- paste0("log_lik[", 1:nrow(dat.culture), "]")
log_lik <- as.data.frame(fit.current$draws(log_lik_pars,format = "matrix"))

# Use loo package to extract elpd, ic, p_loo, MC_se, and pareto_k statistics
loo_current <- loo(as.matrix(log_lik), r_eff = NA)
elpd_loo_current <- loo_current$estimates[1, 1]
elpd_loo_se_current <- loo_current$estimates[1, 2]
ploo_current <- loo_current$estimates[2, 1]
ploo_se_current <- loo_current$estimates[2, 2]
pareto_k_max <- max(pareto_k_values(loo_current))

# Append statistics to performance table
next_row <- c(list(model = "c1",
                   predictors = "T",
                   elpd = elpd_current,
                   elpd_se = elpd_se_current,
                   elpd_loo = elpd_loo_current,
                   elpd_loo_se = elpd_loo_se_current,
                   p_loo_loo = ploo_current,
                   p_loo_loo_se = ploo_se_current,
                   pareto_k_max = pareto_k_max,
                   F_score = F_score,
                   MCC = MCC,
                   n_correct_test = n_corr_test,
                   n_correct_train = n_corr_train,
                   percent_pos = percent_pos,
                   percent_neg = percent_neg,
                   percent_total = percent_total))
tblS4 <- rbind(tblS4, next_row)
View(tblS4)

# Convert log likelihoods to loo object for future reference
loo_name <- paste0("loo_", "c1") # model name for future ref
assign(loo_name, loo_current)



## Run stepwise forward regression ---------------------------------------------

all_predictors <- c("cell_idx", "assay_idx", "log10_dose_pfu", "st_idx", "sp_idx", 
                    "age_idx", "sex_idx", "dpi_idx", "tg_idx")
req_predictors <- c()

for (n_predictors in 1:length(all_predictors)) {
  # Iterate over the number of predictors that can be in the model
  remaining_predictors <- all_predictors[all_predictors %notin% req_predictors]
  cat("Now running models with --", n_predictors, "-- predictors", "\n")
  
  for (ii_pred in 1:length(remaining_predictors)) {
    # Iterate over the possible combinations of remaining predictors
    current_predictors <- c(req_predictors, remaining_predictors[ii_pred]) # list of predictors to use
    mdl_name <- paste0("c", length(current_predictors) + 1, ".", ii_pred) # model name for future ref
    
    # Toggle on/off parameters based on inclusion
    if ("log10_dose_pfu" %in% current_predictors) {dose_pred <- 1}  else {dose_pred <- 0}
    if ("dpi_idx" %in% current_predictors) {dpi_pred <- 1}  else {dpi_pred <- 0}
    if ("sp_idx" %in% current_predictors) {sp_pred <- 1}  else {sp_pred <- 0}
    if ("st_idx" %in% current_predictors) {st_pred <- 1}  else {st_pred <- 0}
    if ("tg_idx" %in% current_predictors) {tg_pred <- 1}  else {tg_pred <- 0}
    if ("age_idx" %in% current_predictors) {age_pred <- 1}  else {age_pred <- 0}
    if ("sex_idx" %in% current_predictors) {sex_pred <- 1}  else {sex_pred <- 0}
    if ("cell_idx" %in% current_predictors) {cell_pred <- 1}  else {cell_pred <- 0}
    if ("assay_idx" %in% current_predictors) {assay_pred <- 1}  else {assay_pred <- 0}
    
    # To calculate number of correctly classified values in train/test data
    n_corr_test <- 0 
    n_corr_train <- 0
    
    # DF for predictions 
    pred.df <- list(true_val = c(), pred_val = c())
    
    # Log likelihood (prevent any assignments percolating across predictors)
    log_lik <- as.data.frame(matrix(nrow = n_chains * n_iter / 2, 
                                    ncol = nrow(dat.culture)))
    
    ### K-fold cross validation ------------------------------------------------
    
    for (k in 1:n_folds) {
      
      cat("Running fold ", k, "for predictor", remaining_predictors[ii_pred],
          "and model", mdl_name, "\n")
      cat("Note: the chosen predictors so far are:", paste0(req_predictors, collapse = ", "), "\n\n")
      
      # For ease of formatting data for stan
      dat.train <- subset(dat.culture, fold != k)
      dat.test <- subset(dat.culture, fold == k)
      
      # Formatting for stan
      dat.stan <- c(list(
        
        # Number of observations, articles, & levels 
        N_test = dim(dat.test)[1],
        N_train = dim(dat.train)[1],
        N_article = max(dat.train$article_idx),
        L_sp = max(dat.culture$sp_idx),
        L_st = max(dat.culture$st_idx),
        L_tg = max(dat.culture$tg_idx),
        L_dpi = max(dat.culture$dpi_idx),
        L_age = 3,
        L_cell = max(dat.culture$cell_idx),
        
        # Toggle predictors on/off
        t_inc = 1,
        sg_inc = 0,
        dose_inc = dose_pred,
        dpi_inc = dpi_pred,
        sp_inc = sp_pred,
        tg_inc = tg_pred,
        st_inc = st_pred,
        age_inc = age_pred,
        sex_inc = sex_pred,
        cell_inc = cell_pred,
        assay_inc = assay_pred,
        
        # Response variables
        pos_inf_train = dat.train$pos_inf,
        pos_inf_test = dat.test$pos_inf,
        
        # Predictor variables
        t_train = dat.train$val_total,
        t_test = dat.test$val_total,
        sg_train = dat.train$val_sg_pred,
        sg_test = dat.test$val_sg_pred,
        dose_train = dat.train$log10_dose_pfu,
        dose_test = dat.test$log10_dose_pfu,
        tg_train = dat.train$tg_idx,
        tg_test = dat.test$tg_idx,
        dpi_train = dat.train$dpi_idx,
        dpi_test = dat.test$dpi_idx,
        age_train = dat.train$age_idx,
        age_test = dat.test$age_idx,
        sex_train = dat.train$sex_idx,
        sex_test = dat.test$sex_idx,
        st_train = dat.train$st_idx,
        st_test = dat.test$st_idx,
        sp_train = dat.train$sp_idx,
        sp_test = dat.test$sp_idx,
        cell_train = dat.train$cell_idx,
        cell_test = dat.test$cell_idx,
        assay_train = dat.train$assay_idx,
        assay_test = dat.test$assay_idx,
        
        # Article hierarchy
        article_train = dat.train$article_idx, 
        article_test = dat.test$article_idx))
      
      # Run current model
      fit.current <- model.culture.xval$sample(data = dat.stan,
                                               chains = n_chains,
                                               parallel_chains = n_chains,
                                               iter_warmup = n_iter / 2,
                                               iter_sampling = n_iter / 2,
                                               refresh = 1000,
                                               seed = seed)
      
      # Store log likelihood for all test data
      log_lik_pars <- paste0("log_lik[", 1:nrow(dat.test), "]")
      log_lik[, dat.culture$fold == k] <- as.data.frame(fit.current$draws(log_lik_pars, 
                                                                          format = "matrix"))
      
      # Determine number of correctly classified samples in test data
      p_pos_test <- paste0("p_pos_test[", 1:nrow(dat.test), "]") 
      p_test_df <- fit.current$draws(p_pos_test, format = "df")
      p_test_mean <- as.data.frame(colMeans(p_test_df))
      p_test_mean <- p_test_mean[rownames(p_test_mean) %notin% c(".chain",
                                                                 ".iteration",
                                                                 ".draw"), ]
      pred_test <- p_test_mean
      pred_test[p_test_mean >= 0.5] <- 1
      pred_test[p_test_mean < 0.5] <- 0
      n_correct <- sum(dat.test$pos_inf == pred_test)
      n_corr_test <- n_corr_test + n_correct
      
      # Add to prediction df
      pred.df$true_val <- c(pred.df$true_val, dat.test$pos_inf)
      pred.df$pred_val <- c(pred.df$pred_val, pred_test)
      
      
      # Determine number of correctly classified samples in training data
      p_pos_train <- paste0("p_pos_train[", 1:nrow(dat.train), "]")
      p_train_df <- fit.current$draws(p_pos_train, format = "df")
      p_train_mean <- as.data.frame(colMeans(p_train_df))
      p_train_mean <- p_train_mean[rownames(p_train_mean) %notin% c(".chain",
                                                                    ".iteration",
                                                                    ".draw"), ]
      pred_train <- p_train_mean
      pred_train[p_train_mean >= 0.5] <- 1
      pred_train[p_train_mean < 0.5] <- 0
      n_correct <- sum(dat.train$pos_inf == pred_train)
      n_corr_train <- n_corr_train + n_correct
    }
    
    # Generate statistics based on equations in Vehtari 2017
    elpd_current <- sum(log(colMeans(exp(as.matrix(log_lik)))))
    elpd_se_current <- sqrt(ncol(log_lik)) * sd(log(colMeans(exp(as.matrix(log_lik)))))
    ploo_current <- sum(colVars(as.matrix(log_lik)))
    ploo_se_current <- sqrt(ncol(log_lik)) * sd(colVars(as.matrix(log_lik)))
    pareto_k_max <- NA
    
    # Generate F-score & MCC
    TP <- sum(pred.df$true_val == 1 & pred.df$pred_val == 1)
    FP <- sum(pred.df$true_val == 0 & pred.df$pred_val == 1)
    TN <- sum(pred.df$true_val == 0 & pred.df$pred_val == 0) + TN_T_neg
    FN <- sum(pred.df$true_val == 1 & pred.df$pred_val == 0) + FN_T_neg
    F_score <- (2 * TP) / (2 * TP + FP + FN)
    MCC <- ((TP * TN) - (FP * FN))/(sqrt(TP + FP)*sqrt(TP + FN)*sqrt(TN + FP)*sqrt(TN + FN))
    
    # Percent of samples correctly classified, stratified by Pos/Neg/All 
    percent_pos <- TP / (TP + FN)  * 100
    percent_neg <- TN / (TN + FP) * 100
    percent_total <- (TN + TP) / (TP + FP + TN + FN) * 100
    
    
    # Store log_lik dataframe for future reference
    loglik_name <- paste0("log_lik_", mdl_name) # model name for future ref
    assign(loglik_name, log_lik)
    
    
    ## PSIS-LOO method ---------------------------------------------------------
    
    cat("Now running PSIS-LOO method for predictor", remaining_predictors[ii_pred],
        "and model", mdl_name, "\n")
    
    # Will run slower, but prevents needing an additional model file
    dat.train <- dat.culture
    dat.test <- dat.culture
    
    # Formatting for stan
    dat.stan <- c(list(
      
      # Number of observations, articles, & levels 
      N_test = dim(dat.test)[1],
      N_train = dim(dat.train)[1],
      N_article = max(dat.culture$article_idx),
      L_sp = max(dat.culture$sp_idx),
      L_st = max(dat.culture$st_idx),
      L_tg = max(dat.culture$tg_idx),
      L_dpi = max(dat.culture$dpi_idx),
      L_age = 3,
      L_cell = max(dat.culture$cell_idx),
      
      # Toggle predictors on/off
      t_inc = 1,
      sg_inc = 0,
      dose_inc = dose_pred,
      dpi_inc = dpi_pred,
      sp_inc = sp_pred,
      tg_inc = tg_pred,
      st_inc = st_pred,
      age_inc = age_pred,
      sex_inc = sex_pred,
      cell_inc = cell_pred,
      assay_inc = assay_pred,
      
      # Response variables
      pos_inf_train = dat.train$pos_inf,
      pos_inf_test = dat.test$pos_inf,
      
      # Predictor variables
      t_train = dat.train$val_total,
      t_test = dat.test$val_total,
      sg_train = dat.train$val_sg_pred,
      sg_test = dat.test$val_sg_pred,
      dose_train = dat.train$log10_dose_pfu,
      dose_test = dat.test$log10_dose_pfu,
      tg_train = dat.train$tg_idx,
      tg_test = dat.test$tg_idx,
      dpi_train = dat.train$dpi_idx,
      dpi_test = dat.test$dpi_idx,
      age_train = dat.train$age_idx,
      age_test = dat.test$age_idx,
      sex_train = dat.train$sex_idx,
      sex_test = dat.test$sex_idx,
      st_train = dat.train$st_idx,
      st_test = dat.test$st_idx,
      sp_train = dat.train$sp_idx,
      sp_test = dat.test$sp_idx,
      cell_train = dat.train$cell_idx,
      cell_test = dat.test$cell_idx,
      assay_train = dat.train$assay_idx,
      assay_test = dat.test$assay_idx,
      
      # Article hierarchy
      article_train = dat.train$article_idx, 
      article_test = dat.test$article_idx))
    
    # Run current model
    fit.current <- model.culture.xval$sample(data = dat.stan,
                                             chains = n_chains,
                                             parallel_chains = n_chains,
                                             iter_warmup = n_iter / 2,
                                             iter_sampling = n_iter / 2,
                                             refresh = 1000,
                                             seed = seed)
    
    # Store log likelihood
    log_lik_pars <- paste0("log_lik[", 1:nrow(dat.culture), "]")
    log_lik <- as.data.frame(fit.current$draws(log_lik_pars,format = "matrix"))
    
    # Use loo package to extract elpd, ic, p_loo, MC_se, and pareto_k statistics
    loo_current <- loo(as.matrix(log_lik), r_eff = NA)
    elpd_loo_current <- loo_current$estimates[1, 1]
    elpd_loo_se_current <- loo_current$estimates[1, 2]
    ploo_current <- loo_current$estimates[2, 1]
    ploo_se_current <- loo_current$estimates[2, 2]
    pareto_k_max <- max(pareto_k_values(loo_current))
    
    # Append statistics to performance table
    next_row <- c(list(model = mdl_name,
                       predictors = paste0(current_predictors, collapse = ", "),
                       elpd = elpd_current,
                       elpd_se = elpd_se_current,
                       elpd_loo = elpd_loo_current,
                       elpd_loo_se = elpd_loo_se_current,
                       p_loo_loo = ploo_current,
                       p_loo_loo_se = ploo_se_current,
                       pareto_k_max = pareto_k_max,
                       F_score = F_score,
                       MCC = MCC,
                       n_correct_test = n_corr_test,
                       n_correct_train = n_corr_train,
                       percent_pos = percent_pos,
                       percent_neg = percent_neg,
                       percent_total = percent_total))
    tblS4 <- rbind(tblS4, next_row)
    View(tblS4)
    
    # Convert log likelihoods to loo object for future reference
    loo_name <- paste0("loo_", mdl_name) # model name for future ref
    assign(loo_name, loo_current)
    
    # Save fit to compare parameter estimates between models
    fit_name <- paste0("fit.cul.", mdl_name)
    assign(fit_name, fit.current)
    
  }
  
  # Find the best predictor for this step, and set as required moving forward
  npred.subset <- subset(tblS4, str_detect(tblS4$model, paste0("c", length(current_predictors) + 1, ".")))
  best_predictors <- npred.subset$predictors[which(npred.subset$elpd == max(npred.subset$elpd))] 
  req_predictors <- unique(c(req_predictors, unlist(strsplit(best_predictors,", "))))
  print(req_predictors)
  cat("The chosen predictors so far are:", paste0(req_predictors, collapse = ", "), "\n")
  
}

tblS4.raw <- tblS4


# Generate Table S4 ------------------------------------------------------------

## Add in elpd_diff ------------------------------------------------------------

best_model <- tblS4$model[which(tblS4$elpd == max(tblS4$elpd))]
best_elpd <- tblS4$elpd[tblS4$model == best_model]
best_log_lik <- get(paste0("log_lik_", best_model))

tblS4$elpd_diff <- 0
tblS4$elpd_diff_se <- 0

for (row_num in 1:nrow(tblS4)) {
  tblS4$elpd_diff[row_num] <- tblS4$elpd[row_num] - best_elpd
  current_log_lik <- get(paste0("log_lik_", tblS4$model[row_num]))
  tblS4$elpd_diff_se[row_num] <- sqrt(nrow(dat)) * sd(log(colMeans(exp(as.matrix(best_log_lik)))) -
                                                        log(colMeans(exp(as.matrix(current_log_lik)))))
}


## Add in loo-based elpd-diff -------------------------------------------------

comp_loo <- loo_compare(x = list("c1" = loo_c1, 
                                 "c2.1" = loo_c2.1, "c2.2" = loo_c2.2, "c2.3" = loo_c2.3,
                                 "c2.4" = loo_c2.4, "c2.5" = loo_c2.5, "c2.6" = loo_c2.6,
                                 "c2.7" = loo_c2.7, "c2.8" = loo_c2.8, "c2.9" = loo_c2.9,
                                 
                                 "c3.1" = loo_c3.1, "c3.2" = loo_c3.2, "c3.3" = loo_c3.3,
                                 "c3.4" = loo_c3.4, "c3.5" = loo_c3.5, "c3.6" = loo_c3.6,
                                 "c3.7" = loo_c3.7, "c3.8" = loo_c3.8,
                                 
                                 "c4.1" = loo_c4.1, "c4.2" = loo_c4.2, "c4.3" = loo_c4.3,
                                 "c4.4" = loo_c4.4, "c4.5" = loo_c4.5, "c4.6" = loo_c4.6,
                                 "c4.7" = loo_c4.7,
                                 
                                 "c5.1" = loo_c5.1, "c5.2" = loo_c5.2, "c5.3" = loo_c5.3,
                                 "c5.4" = loo_c5.4, "c5.5" = loo_c5.5,"c5.6" = loo_c5.6,
                                 
                                 "c6.1" = loo_c6.1, "c6.2" = loo_c6.2, "c6.3" = loo_c6.3,
                                 "c6.4" = loo_c6.4, "c6.5" = loo_c6.5, 
                                 
                                 "c7.1" = loo_c7.1, "c7.2" = loo_c7.2,
                                 "c7.3" = loo_c7.3, "c7.4" = loo_c7.4,
                                 
                                 "c8.1" = loo_c8.1, "c8.2" = loo_c8.2, "c8.3" = loo_c8.3,
                                 
                                 "c9.1" = loo_c9.1, "c9.2" = loo_c9.2,
                                 
                                 "c10.1" = loo_c10.1))


tblS4$elpd_diff_loo <- 0
tblS4$elpd_diff_se_loo <- 0

for (row_num in 1:nrow(comp_loo)) {
  mdl <- rownames(comp_loo)[row_num]
  tblS4$elpd_diff_loo[tblS4$model == mdl] <- comp_loo[row_num, 1]
  tblS4$elpd_diff_se_loo[tblS4$model == mdl] <- comp_loo[row_num, 2]
}



## Compare to "best" model ----------------------------------------------------

# To add rows comparing ELPD between the "best" model
best_model <- "c8.1"
best_elpd <- tblS4$elpd[tblS4$model == best_model]
best_log_lik <- get(paste0("log_lik_", best_model))

tblS4$elpd_diff_best <- 0
tblS4$elpd_diff_best_se <- 0

for (row_num in 1:nrow(tblS4)) {
  tblS4$elpd_diff_best[row_num] <- tblS4$elpd[row_num] - best_elpd
  current_log_lik <- get(paste0("log_lik_", tblS4$model[row_num]))
  tblS4$elpd_diff_best_se[row_num] <- sqrt(nrow(dat)) * sd(log(colMeans(exp(as.matrix(best_log_lik)))) -
                                                             log(colMeans(exp(as.matrix(current_log_lik)))))
}

# Simple ELPD Diff statistic for LOO comparisons of key models
comp_loo_best <- loo_compare(x = list("c1" = loo_c1,
                                      "c8.1" = loo_c8.1, 
                                      "c10.1" = loo_c10.1))
write.csv(comp_loo_best, file = "./outputs/tables/EA-table-culture-inf-key-model-comparison-loo.csv",
          row.names = TRUE)




## Generate percentages for correctly predicted samples -----------------------

tblS4$percent_test <- tblS4$n_correct_test / dim(dat.culture)[1] * 100
tblS4$percent_train <- tblS4$n_correct_train / (dim(dat.culture)[1] * (n_folds - 1)) * 100

# Add correctly predicted RNA negative samples
tblS4$n_correct_test_adj <- tblS4$n_correct_test + TN_T_neg
tblS4$n_correct_train_adj <- tblS4$n_correct_train + (TN_T_neg * (n_folds - 1))

tblS4$percent_test_adj <- tblS4$n_correct_test_adj / (nrow(dat.full)) * 100
tblS4$percent_train_adj <- tblS4$n_correct_train_adj / (nrow(dat.full) * (n_folds - 1)) * 100

# Save raw version
write.csv(tblS4, file = "./outputs/tables/tblS4-culture-totRNA-model-selection-inf-raw.csv",
          row.names = FALSE)



## Aesthetic changes ----------------------------------------------------------

# Round for visualization in table
tblS4$percent_test <- round(tblS4$percent_test, digits = 1)
tblS4$percent_train <- round(tblS4$percent_train, digits = 1)
tblS4$percent_test_adj <- round(tblS4$percent_test_adj, digits = 1)
tblS4$percent_train_adj <- round(tblS4$percent_train_adj, digits = 1)
tblS4$percent_total <- round(tblS4$percent_total, digits = 1)
tblS4$percent_pos <- round(tblS4$percent_pos, digits = 1)
tblS4$percent_neg <- round(tblS4$percent_neg, digits = 1)

tblS4$elpd <- round(tblS4$elpd, digits = 2)
tblS4$elpd_se <- round(tblS4$elpd_se, digits = 2)
tblS4$elpd_diff <- round(tblS4$elpd_diff, digits = 2)
tblS4$elpd_diff_se <- round(tblS4$elpd_diff_se, digits = 2)

tblS4$elpd <- round(tblS4$elpd, digits = 2)
tblS4$elpd_se <- round(tblS4$elpd_se, digits = 2)

tblS4$elpd_loo <- round(tblS4$elpd_loo, digits = 2)
tblS4$elpd_loo_se <- round(tblS4$elpd_loo_se, digits = 2)
tblS4$elpd_diff_loo <- round(tblS4$elpd_diff_loo, digits = 2)
tblS4$elpd_diff_se_loo <- round(tblS4$elpd_diff_se_loo, digits = 2)

tblS4$p_loo_loo <- round(tblS4$p_loo_loo, digits = 2)
tblS4$p_loo_loo_se <- round(tblS4$p_loo_loo_se, digits = 2)
tblS4$pareto_k_max <- round(tblS4$pareto_k_max, digits = 2)

tblS4$F_score <- round(tblS4$F_score, digits = 2)
tblS4$MCC <- round(tblS4$MCC, digits = 2)

# Consolidate predictor names
tblS4$predictors <- str_replace_all(tblS4$predictors, "log10_dose_pfu", "DOSE")
tblS4$predictors <- str_replace_all(tblS4$predictors, "st_idx", "ST")
tblS4$predictors <- str_replace_all(tblS4$predictors, "sp_idx", "SP")
tblS4$predictors <- str_replace_all(tblS4$predictors, "age_idx", "AGE")
tblS4$predictors <- str_replace_all(tblS4$predictors, "sex_idx", "SEX")
tblS4$predictors <- str_replace_all(tblS4$predictors, "dpi_idx", "DPI")
tblS4$predictors <- str_replace_all(tblS4$predictors, "tg_idx", "TG")
tblS4$predictors <- str_replace_all(tblS4$predictors, "cell_idx", "CELL")
tblS4$predictors <- str_replace_all(tblS4$predictors, "assay_idx", "ASSAY")
tblS4$predictors <- str_replace_all(tblS4$predictors, ", ", " + ")
tblS4$predictors[tblS4$model != "c1"] <- paste0("T + ", tblS4$predictors[tblS4$model != "c1"])

# Generate consolidated columns 
tblS4$"ELPD (SE)" <- paste0(round(tblS4$elpd, digits = 2),
                            " (", 
                            round(tblS4$elpd_se, digits = 2), 
                            ")")
tblS4$"ELPD Difference (SE)" <- paste0(round(tblS4$elpd_diff, digits = 2),
                                       " (", 
                                       round(tblS4$elpd_diff_se, digits = 2), 
                                       ")")
tblS4$"ELPD Difference Best (SE)" <- paste0(round(tblS4$elpd_diff_best, digits = 2),
                                            " (", 
                                            round(tblS4$elpd_diff_best_se, digits = 2), 
                                            ")")
tblS4$"LOO - ELPD (SE)" <- paste0(round(tblS4$elpd_loo, digits = 2),
                                  " (", 
                                  round(tblS4$elpd_loo_se, digits = 2), 
                                  ")")
tblS4$"LOO - ELPD Difference (SE)" <- paste0(round(tblS4$elpd_diff_loo, digits = 2),
                                             " (", 
                                             round(tblS4$elpd_diff_se_loo, digits = 2), 
                                             ")")


tblS4$"LOO - Effective parameters (SE)" <- paste0(round(tblS4$p_loo_loo, digits = 2),
                                                  " (", 
                                                  round(tblS4$p_loo_loo_se, digits = 2), 
                                                  ")")



# Rename columns to be more informative
colnames(tblS4)[colnames(tblS4) == "percent_test_adj"] <- "% correctly predicted (test)"
colnames(tblS4)[colnames(tblS4) == "percent_pos"] <- "% correctly predicted (test, culture positive)"
colnames(tblS4)[colnames(tblS4) == "percent_neg"] <- "% correctly predicted (test, culture negative)"
colnames(tblS4)[colnames(tblS4) == "percent_train_adj"] <- "% correctly predicted (train)"
colnames(tblS4)[colnames(tblS4) == "pareto_k_max"] <- "Max. Pareto k"
colnames(tblS4)[colnames(tblS4) == "model"] <- "Model"
colnames(tblS4)[colnames(tblS4) == "predictors"] <- "Predictors"
colnames(tblS4)[colnames(tblS4) == "F_score"] <- "F score"

# Select only the columns of interest, and in the correct order
tblS4 <- tblS4 %>% select(c("Model", "Predictors", 
                            "ELPD Difference (SE)",
                            "ELPD (SE)",
                            "LOO - ELPD Difference (SE)",
                            "LOO - ELPD (SE)", 
                            "LOO - Effective parameters (SE)",
                            "Max. Pareto k",
                            "F score", 
                            "MCC",
                            "% correctly predicted (train)",
                            "% correctly predicted (test)",
                            "% correctly predicted (test, culture positive)",
                            "% correctly predicted (test, culture negative)"))


# Export -----------------------------------------------------------------------

write.csv(tblS4, file = "./outputs/tables/tblS4-culture-totRNA-model-selection-inf-formatted.csv",
          row.names = FALSE)


