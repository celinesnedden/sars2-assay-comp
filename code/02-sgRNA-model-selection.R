# Runs cross validation for the best sgRNA hurdle model
#    First runs the logistic component, then the linear component (which 
#       includes the best logistic component)
# Note that manual changes in the model files are necessary to alternate between
#    informative and non-informative priors

# Prep environment -------------------------------------------------------------

## Install CmdStanR, CmdStan, and RStan for model compilation ------------------

#  We recommend you review the installation instructions at:
#     https://mc-stan.org/cmdstanr/articles/cmdstanr.html 
#     https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started 

library("cmdstanr")
library("rstan")


## Install other packages ------------------------------------------------------

req_pkgs <- c("loo", "stringr", "tidyverse", "gtools", "matrixStats")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependences = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))

`%notin%` <- Negate(`%in%`)  # for convenience


# Set chains, iterations, cores, & seed ----------------------------------------

options(mc.cores = parallel::detectCores())

debug <- FALSE
if (debug == FALSE) {
  n_chains <- 6
  n_iter <- 4000
} else if (debug == TRUE) {
  n_chains <- 1
  n_iter <- 500
}

seed <- 111


# Load data --------------------------------------------------------------------

dat <- read.csv("./data/clean-data.csv")


# Get totRNA negative statistics -----------------------------------------------

dat.neg <- subset(dat, pos_sg %in% c(0, 1) & cens_total == 1)

TN_neg <- nrow(subset(dat.neg, pos_sg == 0))
FN_neg <- nrow(subset(dat.neg, pos_sg == 1))


# Prep data --------------------------------------------------------------------

# Removes data with no sgRNA data and totRNA negative samples (immediately predicted 
#     to also be sgRNA negative)
dat <- subset(dat, pos_sg %in% c(0, 1) & cens_total == 0)

# Screens sgRNA positive samples with no quantitative information out of fitting
#       for the linear component
dat$val_sg[dat$cens_sg == 3] <- -9

# Convert columns to numerics
dat$val_sg <- as.numeric(dat$val_sg)
dat$val_total <- as.numeric(dat$val_total)

# Set articles as numeric factors
dat$article_idx <-  as.numeric(factor(dat$article))


# Logistic model ---------------------------------------------------------------

## Load model ------------------------------------------------------------------

model.log.xval <- cmdstan_model('./code/stan-model-sgRNA-log-xval.stan')
model.log.loo <- cmdstan_model('./code/stan-model-sgRNA-log-loo.stan')


## Assign folds ----------------------------------------------------------------

set.seed(1235)
n_folds <- 10
dat$fold <- kfold_split_stratified(K = n_folds, x = dat$article_idx)


## Set up results table --------------------------------------------------------

tbl_cols <- c("model", "predictors", "elpd", "elpd_se", 
              "elpd_loo", "elpd_loo_se", "p_loo_loo", "p_loo_loo_se", "pareto_k_max", 
              "F_score", "MCC",
              "n_correct_test", "n_correct_train",
              "percent_pos", "percent_neg")

tblS2 <- as.data.frame(matrix(nrow = 0, ncol = length(tbl_cols)))
colnames(tblS2) <- tbl_cols
log_lik <- as.data.frame(matrix(nrow = n_chains * n_iter / 2, ncol = nrow(dat)))


## Run simplest model ----------------------------------------------------------

# Data frame to store predictions 
pred.df <- list(true_val = c(),
                pred_val = c())


### K-fold cross validation ----------------------------------------------------

n_corr_test <- 0
n_corr_train <- 0

for (k in 1:n_folds) {
  
  cat("Running fold ", k, "for the simplest model \n")
  
  # For ease of formatting data for Stan
  dat.train <- subset(dat, fold != k)
  dat.test <- subset(dat, fold == k)
  
  # Formatting the data for Stan
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
    
    # Toggle predictors on/off
    dose_inc = 0,
    dpi_inc = 0,
    sp_inc = 0,
    tg_inc = 0,
    st_inc = 0,
    age_inc = 0,
    sex_inc = 0,
    article_inc = 0,
    
    # Response variables
    sg_pos_train = dat.train$pos_sg,
    sg_pos_test = dat.test$pos_sg,
    
    # Predictor variables
    t_train = dat.train$val_total,
    t_test = dat.test$val_total,
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
    
    # Article hierarchy
    article_train = dat.train$article_idx, 
    article_test = dat.test$article_idx))
  
  # Run current model
  fit.current <- model.log.xval$sample(data = dat.stan,
                                       chains = n_chains,
                                       parallel_chains = n_chains,
                                       iter_warmup = n_iter / 2,
                                       iter_sampling = n_iter / 2,
                                       refresh = 1000,
                                       seed = seed)
  
  
  # Store log likelihood for all test data
  log_lik_pars <- paste0("log_lik[", 1:nrow(dat.test), "]")
  log_lik[, dat$fold == k] <- as.data.frame(fit.current$draws(log_lik_pars, 
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
  n_correct <- sum(dat.test$pos_sg == pred_test)
  n_corr_test <- n_corr_test + n_correct
  
  # Add to prediction df
  pred.df$true_val <- c(pred.df$true_val, dat.test$pos_sg)
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
  n_correct <- sum(dat.train$pos_sg == pred_train)
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

# Generate F-score
TP <- sum(pred.df$true_val == 1 & pred.df$pred_val == 1)
FP <- sum(pred.df$true_val == 0 & pred.df$pred_val == 1)
TN <- sum(pred.df$true_val == 0 & pred.df$pred_val == 0) + TN_neg
FN <- sum(pred.df$true_val == 1 & pred.df$pred_val == 0) + FN_neg
F_score <- (2 * TP) / (2 * TP + FP + FN)
MCC <- ((TP * TN) - (FP * FN))/(sqrt(TP + FP)*sqrt(TP + FN)*sqrt(TN + FP)*sqrt(TN + FN))

# Percent positives and negatives
percent_pos <- TP / (TP + FN)  * 100
percent_neg <- TN / (TN + FP) * 100

# Store log_lik data frame for future reference
loglik_name <- paste0("log_lik_", "l1") 
assign(loglik_name, log_lik)


### PSIS-LOO method ------------------------------------------------------------

# Formatting the data for Stan
dat.stan <- c(list(
  
  # Number of observations, articles, & levels 
  N = nrow(dat),
  N_lin = nrow(subset(dat, cens_sg == 0)),
  N_article = max(dat$article_idx),
  L_sp = max(dat$sp_idx),
  L_st = max(dat$st_idx),
  L_tg = max(dat$tg_idx),
  L_dpi = max(dat$dpi_idx),
  L_age = max(dat$age_idx),
  
  # Toggle predictors on/off
  dose_inc = 0,
  dpi_inc = 0,
  sp_inc = 0,
  tg_inc = 0,
  st_inc = 0,
  age_inc = 0,
  sex_inc = 0,
  article_inc = 0,
  
  # Response variables
  sg_pos = dat$pos_sg,
  
  # Predictor variables
  t = dat$val_total,
  dose = dat$log10_dose_pfu,
  tg = dat$tg_idx,
  dpi = dat$dpi_idx,
  age = dat$age_idx,
  sex = dat$sex_idx,
  st = dat$st_idx,
  sp = dat$sp_idx,
  
  # Housekeeping for tRNA censoring
  t_cens = dat$cens_total,
  
  # Article hierarchy
  article = dat$article_idx))


# Run current model
fit.current <- model.log.loo$sample(data = dat.stan,
                                    chains = n_chains,
                                    parallel_chains = n_chains,
                                    iter_warmup = n_iter / 2,
                                    iter_sampling = n_iter / 2,
                                    refresh = 1000,
                                    seed = seed)

# Store log likelihood
log_lik_pars <- paste0("log_lik[", 1:nrow(dat), "]")
log_lik <- as.data.frame(fit.current$draws(log_lik_pars,format = "matrix"))

# Use loo package to extract elpd, ic, p_loo, MC_se, and pareto_k statistics
loo_current <- loo(as.matrix(log_lik), r_eff = NA)
elpd_loo_current <- loo_current$estimates[1, 1]
elpd_loo_se_current <- loo_current$estimates[1, 2]
ploo_current <- loo_current$estimates[2, 1]
ploo_se_current <- loo_current$estimates[2, 2]
pareto_k_max <- max(pareto_k_values(loo_current))

# Append statistics to performance table
next_row <- c(list(model = "l1",
                   predictors = "",
                   elpd = elpd_current,
                   elpd_se = elpd_se_current,
                   elpd_loo = elpd_loo_current,
                   elpd_loo_se = elpd_loo_se_current,
                   p_loo_loo = ploo_current,
                   p_loo_loo_se = ploo_se_current,
                   pareto_k_max = pareto_k_max,
                   n_correct_test = n_corr_test,
                   n_correct_train = n_corr_train,
                   percent_pos = percent_pos,
                   percent_neg = percent_neg,
                   F_score = F_score,
                   MCC = MCC))
tblS2 <- rbind(tblS2, next_row)
View(tblS2)

# Convert log likelihoods to loo object for future reference
loo_name <- paste0("loo_", "l1")
assign(loo_name, loo_current)


## Run stepwise forward regression ---------------------------------------------

all_predictors <- c("log10_dose_pfu", "st_idx", "sp_idx", 
                    "age_idx", "sex_idx", "dpi_idx", "tg_idx")
req_predictors <- c()

for (n_predictors in 1:length(all_predictors)) {
  # Iterate over the number of predictors that can be in the model
  remaining_predictors <- all_predictors[all_predictors %notin% req_predictors]
  cat("Now running models with --", n_predictors, "-- predictors", "\n")
  
  for (ii_pred in 1:length(remaining_predictors)) {
    # Iterate over the possible combinations of remaining predictors
    current_predictors <- c(req_predictors, remaining_predictors[ii_pred]) # list of predictors to use
    mdl_name <- paste0("l", length(current_predictors) + 1, ".", ii_pred) # model name for future ref
    
    # Toggle on/off parameters based on inclusion
    if ("article_idx" %in% current_predictors) {article_pred <- 1}  else {article_pred <- 0}
    if ("log10_dose_pfu" %in% current_predictors) {dose_pred <- 1}  else {dose_pred <- 0}
    if ("dpi_idx" %in% current_predictors) {dpi_pred <- 1}  else {dpi_pred <- 0}
    if ("sp_idx" %in% current_predictors) {sp_pred <- 1}  else {sp_pred <- 0}
    if ("st_idx" %in% current_predictors) {st_pred <- 1}  else {st_pred <- 0}
    if ("tg_idx" %in% current_predictors) {tg_pred <- 1}  else {tg_pred <- 0}
    if ("age_idx" %in% current_predictors) {age_pred <- 1}  else {age_pred <- 0}
    if ("sex_idx" %in% current_predictors) {sex_pred <- 1}  else {sex_pred <- 0}
    
    # To calculate number of correctly classified values in train/test data
    n_corr_test <- 0 
    n_corr_train <- 0
    pred.df <- list(true_val = c(),
                    pred_val = c())
    
    
    ### K-fold cross validation ------------------------------------------------
    
    for (k in 1:n_folds) {
      
      cat("Running fold ", k, "for predictor", remaining_predictors[ii_pred],
          "and model", mdl_name, "\n")
      cat("Note: the chosen predictors so far are:", 
          paste0(req_predictors, collapse = ", "), "\n")
      
      # For ease of formatting data for stan
      dat.train <- subset(dat, fold != k)
      dat.test <- subset(dat, fold == k)
      
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
        L_age = max(dat$age_idx),
        
        # Toggle predictors on/off
        dose_inc = dose_pred,
        dpi_inc = dpi_pred,
        sp_inc = sp_pred,
        tg_inc = tg_pred,
        st_inc = st_pred,
        age_inc = age_pred,
        sex_inc = sex_pred,
        article_inc = article_pred,
        
        # Response variables
        sg_pos_train = dat.train$pos_sg,
        sg_pos_test = dat.test$pos_sg,
        
        # Predictor variables
        t_train = dat.train$val_total,
        t_test = dat.test$val_total,
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
        
        # Housekeeping for tRNA censoring
        t_cens_train = dat.train$cens_total,
        t_cens_test = dat.test$cens_total,
        
        # Article hierarchy
        article_train = dat.train$article_idx, 
        article_test = dat.test$article_idx))
      
      # Run current model
      fit.current <- model.log.xval$sample(data = dat.stan,
                                                  chains = n_chains,
                                                  parallel_chains = n_chains,
                                                  iter_warmup = n_iter / 2,
                                                  iter_sampling = n_iter / 2,
                                                  refresh = 1000,
                                                  seed = seed)
      
      # Store log likelihood for all test data
      log_lik_pars <- paste0("log_lik[", 1:nrow(dat.test), "]")
      log_lik[, dat$fold == k] <- as.data.frame(fit.current$draws(log_lik_pars, 
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
      n_correct <- sum(dat.test$pos_sg == pred_test)
      n_corr_test <- n_corr_test + n_correct
      
      # Add to prediction df
      pred.df$true_val <- c(pred.df$true_val, dat.test$pos_sg)
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
      n_correct <- sum(dat.train$pos_sg == pred_train)
      n_corr_train <- n_corr_train + n_correct
    }
    
    # Generate statistics based on equations in Vehtari et al. 2017
    elpd_current <- sum(log(colMeans(exp(as.matrix(log_lik)))))
    elpd_se_current <- sqrt(ncol(log_lik)) * sd(log(colMeans(exp(as.matrix(log_lik)))))
    ploo_current <- sum(colVars(as.matrix(log_lik)))
    ploo_se_current <- sqrt(ncol(log_lik)) * sd(colVars(as.matrix(log_lik)))
    pareto_k_max <- NA
    
    # Generate F-score
    TP <- sum(pred.df$true_val == 1 & pred.df$pred_val == 1)
    FP <- sum(pred.df$true_val == 0 & pred.df$pred_val == 1)
    TN <- sum(pred.df$true_val == 0 & pred.df$pred_val == 0) + TN_neg
    FN <- sum(pred.df$true_val == 1 & pred.df$pred_val == 0) + FN_neg
    F_score <- (2 * TP) / (2 * TP + FP + FN)
    MCC <- ((TP * TN) - (FP * FN))/(sqrt(TP + FP)*sqrt(TP + FN)*sqrt(TN + FP)*sqrt(TN + FN))
    
    # Percent positives and negatives
    percent_pos <- TP / (TP + FN)  * 100
    percent_neg <- TN / (TN + FP) * 100
    
    # Store log_lik dataframe for future reference
    loglik_name <- paste0("log_lik_", mdl_name) # model name for future ref
    assign(loglik_name, log_lik)
    
    
    ## PSIS-LOO method ---------------------------------------------------------
    
    cat("Now running PSIS-LOO method for predictor", remaining_predictors[ii_pred],
        "and model", mdl_name, "\n")
    
    # Formatting for stan
    dat.stan <- c(list(
      
      # Number of observations, articles, & levels 
      N = nrow(dat),
      N_article = max(dat$article_idx),
      L_sp = max(dat$sp_idx),
      L_st = max(dat$st_idx),
      L_tg = max(dat$tg_idx),
      L_dpi = max(dat$dpi_idx),
      L_age = max(dat$age_idx),
      
      # Toggle predictors on/off
      dose_inc = dose_pred,
      dpi_inc = dpi_pred,
      sp_inc = sp_pred,
      tg_inc = tg_pred,
      st_inc = st_pred,
      age_inc = age_pred,
      sex_inc = sex_pred,
      article_inc = article_pred,
      
      # Response variables
      sg_pos = dat$pos_sg,
      
      # Predictor variables
      t = dat$val_total,
      dose = dat$log10_dose_pfu,
      tg = dat$tg_idx,
      dpi = dat$dpi_idx,
      age = dat$age_idx,
      sex = dat$sex_idx,
      st = dat$st_idx,
      sp = dat$sp_idx,
      
      # Housekeeping for tRNA censoring
      t_cens = dat$cens_total,
      
      # Article hierarchy
      article = dat$article_idx))
    
    # Run current model
    fit.current <- model.log.loo$sample(data = dat.stan,
                                           chains = n_chains,
                                           parallel_chains = n_chains,
                                           iter_warmup = n_iter / 2,
                                           iter_sampling = n_iter / 2,
                                           refresh = 1000,
                                           seed = seed)
    
    # Store log likelihood
    log_lik_pars <- paste0("log_lik[", 1:nrow(dat), "]")
    log_lik <- as.data.frame(fit.current$draws(log_lik_pars,format = "matrix"))
    
    # Use loo package to extract statistics
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
                       n_correct_test = n_corr_test,
                       n_correct_train = n_corr_train,
                       percent_pos = percent_pos,
                       percent_neg = percent_neg,
                       F_score = F_score,
                       MCC = MCC))
    tblS2 <- rbind(tblS2, next_row)
    View(tblS2)
    
    # Convert log likelihoods to loo object for future reference
    loo_name <- paste0("loo_", mdl_name) # model name for future ref
    assign(loo_name, loo_current)
    
    }

  # Find the best predictor for this step, and set as required moving forward
  npred.subset <- subset(tblS2, str_detect(tblS2$model, paste0("l", length(current_predictors) + 1, ".")))
  best_predictors <- npred.subset$predictors[which(npred.subset$elpd == max(npred.subset$elpd))] 
  req_predictors <- unique(c(req_predictors, unlist(strsplit(best_predictors,", "))))
  cat("The chosen predictors so far are:", paste0(req_predictors, collapse = ", "), "\n")

}

tblS2.raw <- tblS2



## Generate Table S2 -----------------------------------------------------------

### Add in elpd_diff -----------------------------------------------------------

best_model <- tblS2$model[which(tblS2$elpd == max(tblS2$elpd))]
best_elpd <- tblS2$elpd[tblS2$model == best_model]
best_log_lik <- get(paste0("log_lik_", best_model))

tblS2$elpd_diff <- 0
tblS2$elpd_diff_se <- 0

for (row_num in 1:nrow(tblS2)) {
 tblS2$elpd_diff[row_num] <- tblS2$elpd[row_num] - best_elpd
 current_log_lik <- get(paste0("log_lik_", tblS2$model[row_num]))
 tblS2$elpd_diff_se[row_num] <- sqrt(nrow(dat)) * sd(log(colMeans(exp(as.matrix(best_log_lik)))) -
                                                       log(colMeans(exp(as.matrix(current_log_lik)))))
}


### Add in loo-based elpd-diff -------------------------------------------------

comp_loo <- loo_compare(x = list("l1" = loo_l1,
                                 "l2.1" = loo_l2.1, "l2.2" = loo_l2.2, "l2.3" = loo_l2.3,
                                 "l2.4" = loo_l2.4, "l2.5" = loo_l2.5, "l2.6" = loo_l2.6,
                                 "l2.7" = loo_l2.7, 
                                 "l3.1" = loo_l3.1, "l3.2" = loo_l3.2, "l3.3" = loo_l3.3,
                                 "l3.4" = loo_l3.4, "l3.5" = loo_l3.5, "l3.6" = loo_l3.6,
                                 "l4.1" = loo_l4.1, "l4.2" = loo_l4.2, "l4.3" = loo_l4.3,
                                 "l4.4" = loo_l4.4, "l4.5" = loo_l4.5, 
                                 "l5.1" = loo_l5.1, "l5.2" = loo_l5.2, "l5.3" = loo_l5.3,
                                 "l5.4" = loo_l5.4,
                                 "l6.1" = loo_l6.1, "l6.2" = loo_l6.2, "l6.3" = loo_l6.3,
                                 "l7.1" = loo_l7.1, "l7.2" = loo_l7.2,
                                 "l8.1" = loo_l8.1
                                 ))


tblS2$elpd_diff_loo <- 0
tblS2$elpd_diff_se_loo <- 0

for (row_num in 1:nrow(comp_loo)) {
  mdl <- rownames(comp_loo)[row_num]
  tblS2$elpd_diff_loo[tblS2$model == mdl] <- comp_loo[row_num, 1]
  tblS2$elpd_diff_se_loo[tblS2$model == mdl] <- comp_loo[row_num, 2]
}



### Compare to "best" model ----------------------------------------------------

# To add rows comparing ELPD between the "best" model
best_model <- "l4.2"
best_elpd <- tblS2$elpd[tblS2$model == best_model]
best_log_lik <- get(paste0("log_lik_", best_model))

tblS2$elpd_diff_best <- 0
tblS2$elpd_diff_best_se <- 0

for (row_num in 1:nrow(tblS2)) {
  tblS2$elpd_diff_best[row_num] <- tblS2$elpd[row_num] - best_elpd
  current_log_lik <- get(paste0("log_lik_", tblS2$model[row_num]))
  tblS2$elpd_diff_best_se[row_num] <- sqrt(nrow(dat)) * sd(log(colMeans(exp(as.matrix(best_log_lik)))) -
                                                        log(colMeans(exp(as.matrix(current_log_lik)))))
}

# Simple ELPD Diff statistic for LOO comparisons of key models
comp_loo_best <- loo_compare(x = list("l1" = loo_l1,
                                      "l4.2" = loo_l4.2, 
                                      "l8.1" = loo_l8.1))
write.csv(comp_loo_best, file = "./outputs/tables/EA-table-sgRNA-log-inf-key-model-comparison-loo.csv",
          row.names = TRUE)


### Generate percentages for correctly predicted samples -----------------------

# Add correctly predicted RNA negative samples to statistics
tblS2$n_correct_test_adj <- tblS2$n_correct_test + TN_neg
tblS2$n_correct_train_adj <- tblS2$n_correct_train + (TN_neg * (n_folds - 1))

# Percent correct 
tblS2$percent_test <- tblS2$n_correct_test_adj / (nrow(dat) + nrow(dat.neg)) * 100
tblS2$percent_train <- tblS2$n_correct_train_adj / ((nrow(dat) * (n_folds - 1)) + ((n_folds - 1) * nrow(dat.neg))) * 100


### Update predictor names -----------------------------------------------------

tblS2$predictors <- str_remove_all(tblS2$predictors, "_idx")
tblS2$predictors <- str_remove_all(tblS2$predictors, "log10_|_pfu")
tblS2$predictors <- str_replace_all(tblS2$predictors, ",", " +" )
tblS2$predictors <- toupper(tblS2$predictors)


### Save a raw version ---------------------------------------------------------

write.csv(tblS2, file = "./outputs/tables/tblS2-sgRNA-log-model-selection-inf-raw.csv",
          row.names = FALSE)



### Aesthetic changes for table ------------------------------------------------

# Rounding
tblS2$percent_test <- round(tblS2$percent_test, digits = 2)
tblS2$percent_train <- round(tblS2$percent_train, digits = 2)
tblS2$percent_pos <- round(tblS2$percent_pos, digits = 2)
tblS2$percent_neg <- round(tblS2$percent_neg, digits = 2)

tblS2$elpd <- round(tblS2$elpd, digits = 2)
tblS2$elpd_se <- round(tblS2$elpd_se, digits = 2)
tblS2$elpd_diff <- round(tblS2$elpd_diff, digits = 2)
tblS2$elpd_diff_se <- round(tblS2$elpd_diff_se, digits = 2)
tblS2$elpd_diff_best <- round(tblS2$elpd_diff_best, digits = 2)
tblS2$elpd_diff_best_se <- round(tblS2$elpd_diff_best_se, digits = 2)

tblS2$elpd_loo <- round(tblS2$elpd_loo, digits = 2)
tblS2$elpd_loo_se <- round(tblS2$elpd_loo_se, digits = 2)
tblS2$elpd_diff_loo <- round(tblS2$elpd_diff_loo, digits = 2)
tblS2$elpd_diff_se_loo <- round(tblS2$elpd_diff_se_loo, digits = 2)

tblS2$p_loo_loo <- round(tblS2$p_loo_loo, digits = 2)
tblS2$p_loo_loo_se <- round(tblS2$p_loo_loo_se, digits = 2)
tblS2$pareto_k_max <- round(tblS2$pareto_k_max, digits = 2)

tblS2$F_score <- round(tblS2$F_score, digits = 2)
tblS2$MCC <- round(tblS2$MCC, digits = 2)

# Generate consolidated columns 
tblS2$"ELPD (SE)" <- paste0(round(tblS2$elpd, digits = 2),
                           " (", 
                           round(tblS2$elpd_se, digits = 2), 
                           ")")
tblS2$"ELPD Difference (SE)" <- paste0(round(tblS2$elpd_diff, digits = 2),
                                       " (", 
                                       round(tblS2$elpd_diff_se, digits = 2), 
                                       ")")

tblS2$"LOO - ELPD (SE)" <- paste0(round(tblS2$elpd_loo, digits = 2),
                            " (", 
                            round(tblS2$elpd_loo_se, digits = 2), 
                            ")")
tblS2$"LOO - ELPD Difference (SE)" <- paste0(round(tblS2$elpd_diff_loo, digits = 2),
                                       " (", 
                                       round(tblS2$elpd_diff_se_loo, digits = 2), 
                                       ")")


tblS2$"LOO - Effective parameters (SE)" <- paste0(round(tblS2$p_loo_loo, digits = 2),
                                            " (", 
                                            round(tblS2$p_loo_loo_se, digits = 2), 
                                            ")")

tblS2$"ELPD Difference Best (SE)" <- paste0(round(tblS2$elpd_diff_best, digits = 2),
                                            " (", 
                                            round(tblS2$elpd_diff_best_se, digits = 2), 
                                            ")")

# Rename columns to be more informative
colnames(tblS2)[colnames(tblS2) == "percent_test"] <- "% correctly predicted (test)"
colnames(tblS2)[colnames(tblS2) == "percent_pos"] <- "% correctly predicted (positive, test)"
colnames(tblS2)[colnames(tblS2) == "percent_neg"] <- "% correctly predicted (negative, test)"
colnames(tblS2)[colnames(tblS2) == "percent_train"] <- "% correctly predicted (train)"
colnames(tblS2)[colnames(tblS2) == "pareto_k_max"] <- "Max. Pareto k"
colnames(tblS2)[colnames(tblS2) == "model"] <- "Model"
colnames(tblS2)[colnames(tblS2) == "predictors"] <- "Predictors"
colnames(tblS2)[colnames(tblS2) == "F_score"] <- "F score"

# Order by decreasing elpd value
tblS2 <- tblS2[order(tblS2$Model), ]

# Select only the columns of interest, and in the correct order
tblS2 <- tblS2 %>% dplyr::select(c("Model", "Predictors", 
                                   "ELPD Difference (SE)",
                                   "ELPD (SE)",
                                   "ELPD Difference Best (SE)",
                                   "LOO - ELPD Difference (SE)",
                                   "LOO - ELPD (SE)", 
                                   "LOO - Effective parameters (SE)",
                                   "Max. Pareto k",
                                   "F score", 
                                   "MCC",
                                   "% correctly predicted (train)",
                                   "% correctly predicted (test)",
                                   "% correctly predicted (positive, test)",
                                   "% correctly predicted (negative, test)"))


## Save table ------------------------------------------------------------------

write.csv(tblS2, file = "./outputs/tables/tblS2-sgRNA-log-model-selection-inf-formatted.csv",
          row.names = FALSE)



# Full hurdle model ------------------------------------------------------------

## Load models -----------------------------------------------------------------

model.full.xval <- cmdstan_model('./code/stan-model-sgRNA-full-xval.stan')
model.full.loo <- stan_model('./code/stan-model-sgRNA-full-loo-rstan.stan')
model.simple.xval <- cmdstan_model('./code/stan-model-sgRNA-simplest-xval.stan')
model.simple.loo <- cmdstan_model('./code/stan-model-sgRNA-simplest.stan')


## Assign folds ----------------------------------------------------------------

n_folds <- 10
set.seed(1235)
dat$fold <- kfold_split_stratified(K = n_folds, x = dat$article_idx)


## Set up results table --------------------------------------------------------

tbl_cols <- c("model", "predictors", "elpd", "elpd_se", "p_loo", "p_loo_se",
              "looic", "looic_se", "elpd_loo", "elpd_loo_se", "p_loo_loo",
              "p_loo_loo_se", "pareto_k_max", "mae_test", "mae_scaled_test",
              "mae_train", "mae_scaled_train",
              "n_within_50_test", "n_within_95_test", 
              "n_within_50_train", "n_within_95_train")
tblS3 <- as.data.frame(matrix(nrow = 0, ncol = length(tbl_cols)))
colnames(tblS3) <- tbl_cols
log_lik_log <- as.data.frame(matrix(nrow = n_chains * n_iter / 2, ncol = nrow(dat)))
log_lik_lin <- as.data.frame(matrix(nrow = n_chains * n_iter / 2, ncol = nrow(dat)))


## Set up df to store predictions ----------------------------------------------

pred.df <- list(true_val = c(),
                pred_val = c())


## Run simplest model ----------------------------------------------------------

### K-fold cross validation ----------------------------------------------------

n_within_95_test <- 0
n_within_50_test <- 0
n_within_95_train <- 0
n_within_50_train <- 0
n_corr_test <- 0
n_corr_train <- 0
mae_test <- c()
mae_scaled_test <- c()
mae_train <- c()
mae_scaled_train <- c()

for (k in 1:n_folds) {
  
  cat("Running fold ", k, "for the simplest model \n")
  
  # For ease of formatting data for stan
  dat.train <- subset(dat, fold != k)
  dat.test <- subset(dat, fold == k)
  
  # Formatting for stan
  dat.stan <- c(list(
    
    # Number of observations, articles, & levels 
    N_test = dim(dat.test)[1],
    N_train = dim(dat.train)[1],
    N_test_lin = nrow(subset(dat.test, cens_sg == 0)),
    
    # Response variables
    pos_sg_train = dat.train$pos_sg,
    pos_sg_test = dat.test$pos_sg,
    val_sg_train = dat.train$val_sg,
    val_sg_test = dat.test$val_sg,
    
    # Predictor variables
    t_train = dat.train$val_total,
    t_test = dat.test$val_total))
  
  # Run current model
  fit.current <- model.simple.xval$sample(data = dat.stan,
                                          chains = n_chains,
                                          parallel_chains = n_chains,
                                          iter_warmup = n_iter / 2,
                                          iter_sampling = n_iter / 2,
                                          refresh = 1000,
                                          seed = seed)
  
  # Store log likelihood for all test data
  log_lik_log_pars <- paste0("log_lik_log[", 1:nrow(dat.test), "]")
  log_lik_log[, dat$fold == k] <- as.data.frame(fit.current$draws(log_lik_log_pars, 
                                                              format = "matrix"))
  log_lik_lin_pars <- paste0("log_lik_lin[", 1:nrow(subset(dat.test, cens_sg == 0)), "]")
  log_lik_lin[, dat$fold == k & dat$cens_sg == 0] <- as.data.frame(fit.current$draws(log_lik_lin_pars, 
                                                                                 format = "matrix"))
  
  
  # TEST DATA
  dat.test.obs <- subset(dat.test, cens_sg == 0)
  
  # Determine number of samples within 50 and 95% prediction interval
  pred_val_test <- paste0("pred_val_test[", which(dat.test$cens_sg == 0), "]") 
  pred_test_df <- as.data.frame(fit.current$draws(pred_val_test, format = "df"))
  pred_test_df <- subset(pred_test_df, select = pred_val_test) #remove .chain, ...
  
  pred_0.025 <- sapply(pred_test_df, function(x) quantile(x, probs = 0.025))
  pred_0.25 <- sapply(pred_test_df, function(x) quantile(x, probs = 0.25))
  pred_0.75 <- sapply(pred_test_df, function(x) quantile(x, probs = 0.75))
  pred_0.975 <- sapply(pred_test_df, function(x) quantile(x, probs = 0.975))
  
  n_within_95 <- sum(pred_0.025 <= dat.test.obs$val_sg & dat.test.obs$val_sg <= pred_0.975)
  n_within_50 <- sum(pred_0.25 <= dat.test.obs$val_sg & dat.test.obs$val_sg <= pred_0.75)
  n_within_95_test <- n_within_95_test + n_within_95 
  n_within_50_test <- n_within_50_test + n_within_50 
  
  # Calculate MAE and MAE scaled
  pred_median <- sapply(pred_test_df, function(x) median(x))
  pred_sd <- sapply(pred_test_df, function(x) sd(x))
  mae_test <- c(mae_test, abs(dat.test.obs$val_sg - pred_median))
  mae_scaled_test <- c(mae_scaled_test, abs(dat.test.obs$val_sg - pred_median)/pred_sd)

  # TRAIN DATA
  dat.train.obs <- subset(dat.train, cens_sg == 0)
  
  # Determine number of samples within 50 and 95% prediction interval
  pred_val_train <- paste0("pred_val_train[", which(dat.train$cens_sg == 0), "]") 
  pred_train_df <- as.data.frame(fit.current$draws(pred_val_train, format = "df"))
  pred_train_df <- subset(pred_train_df, select = pred_val_train) #remove .chain, ...
  
  pred_0.025 <- sapply(pred_train_df, function(x) quantile(x, probs = 0.025))
  pred_0.25 <- sapply(pred_train_df, function(x) quantile(x, probs = 0.25))
  pred_0.75 <- sapply(pred_train_df, function(x) quantile(x, probs = 0.75))
  pred_0.975 <- sapply(pred_train_df, function(x) quantile(x, probs = 0.975))
  
  n_within_95 <- sum(pred_0.025 <= dat.train.obs$val_sg & dat.train.obs$val_sg <= pred_0.975)
  n_within_50 <- sum(pred_0.25 <= dat.train.obs$val_sg & dat.train.obs$val_sg <= pred_0.75)
  n_within_95_train <- n_within_95_train + n_within_95 
  n_within_50_train <- n_within_50_train + n_within_50 
  
  # Add MAE and MAE scaled terms 
  pred_median <- sapply(pred_train_df, function(x) median(x))
  pred_sd <- sapply(pred_train_df, function(x) sd(x))
  mae_train <- c(mae_train, abs(dat.train.obs$val_sg - pred_median))
  mae_scaled_train <- c(mae_scaled_train, abs(dat.train.obs$val_sg - pred_median)/pred_sd)
  
  
}

log_lik <- cbind(log_lik_log, subset(log_lik_lin, select = c(which(!is.na(log_lik_lin[1, ])))))

# Generate ELPD & similar statistics based on equations in Vehtari et al. 2017
elpd_current <- sum(log(colMeans(exp(as.matrix(log_lik)))))
elpd_se_current <- sqrt(ncol(log_lik)) * sd(log(colMeans(exp(as.matrix(log_lik)))))
looic_current <- -2 * elpd_current 
looic_se_current <- 2 * elpd_se_current
ploo_current <- sum(colVars(as.matrix(log_lik)))
ploo_se_current <- sqrt(ncol(log_lik)) * sd(colVars(as.matrix(log_lik)))

# Generate MAE statistics
mae_train_full <- median(mae_train)
mae_scaled_train_full <- median(mae_scaled_train)
mae_test_full <- median(mae_test)
mae_scaled_test_full <- median(mae_scaled_test)

# Store log_lik dataframe to use for model comparison
loglik_name <- paste0("log_lik_", "f1") # model name for future ref
assign(loglik_name, log_lik)



### PSIS-LOO method ------------------------------------------------------------

# Run model with no test data to generate ELPD & other statistics with LOO package

# Prep data for Stan
dat.stan <- c(list(
  
  # Number of observations, articles, & levels 
  N = nrow(dat),
  N_lin = nrow(subset(dat, cens_sg == 0)),
  
  # Response variables
  pos_sg = dat$pos_sg,
  val_sg = dat$val_sg,
  
  # Predictor variables
  t_val = dat$val_total))


## Run the model
fit.current <- model.simple.loo$sample(data = dat.stan,
                                       chains = n_chains,
                                       parallel_chains = n_chains,
                                       iter_warmup = n_iter / 2,
                                       iter_sampling = n_iter / 2,
                                       refresh = 1000)

# Extract log likelihood matrix
log_lik_log_pars <- paste0("log_lik_log[", 1:nrow(dat), "]")
log_lik_log <- as.data.frame(fit.current$draws(log_lik_log_pars, format = "matrix"))
log_lik_lin_pars <- paste0("log_lik_lin[", 1:nrow(subset(dat, cens_sg == 0)), "]")
log_lik_lin <- as.data.frame(fit.current$draws(log_lik_lin_pars, format = "matrix"))
log_lik <- cbind(log_lik_log, subset(log_lik_lin, select = c(which(!is.na(log_lik_lin[1, ])))))

# Obtain PSIS-based performance statistics
loo_current <- loo(as.matrix(log_lik), r_eff = NA)
elpd_loo_current <- loo_current$estimates[1, 1]
elpd_loo_se_current <- loo_current$estimates[1, 2]
ploo_loo_current <- loo_current$estimates[2, 1]
ploo_loo_se_current <- loo_current$estimates[2, 2]
pareto_k_max <- max(pareto_k_values(loo_current))

# Append statistics to performance table
next_row <- c(list(model = "f1",
                   predictors = "",
                   elpd = elpd_current,
                   elpd_se = elpd_se_current,
                   p_loo = ploo_current,
                   p_loo_se = ploo_se_current,
                   looic = looic_current,
                   looic_se = looic_se_current,
                   elpd_loo = elpd_loo_current,
                   elpd_loo_se = elpd_loo_se_current,
                   p_loo_loo = ploo_loo_current,
                   p_loo_loo_se = ploo_loo_se_current,
                   pareto_k_max = pareto_k_max,
                   mae_test = mae_test_full,
                   mae_scaled_test = mae_scaled_test_full,
                   mae_train = mae_train_full,
                   mae_scaled_train = mae_scaled_train_full,
                   n_within_50_test = n_within_50_test,
                   n_within_50_train = n_within_50_train,
                   n_within_95_test = n_within_95_test,
                   n_within_95_train = n_within_95_train))
tblS3 <- rbind(tblS3, next_row)
View(tblS3)

# Convert log likelihoods to loo object for future reference
loo_name <- paste0("loo_", "f1") # model name for future ref
assign(loo_name, loo_current)


## Run stepwise forward regression ---------------------------------------------

all_predictors <- c("log10_dose_pfu", "st_idx", "sp_idx", "age_idx", 
                    "sex_idx", "dpi_idx", "tg_idx")
req_predictors <- c()


for (n_predictors in 1:length(all_predictors)) {
  # Iterate over the number of predictors that can be in the model
  remaining_predictors <- all_predictors[all_predictors %notin% req_predictors]
  cat("Now running models with --", n_predictors, "-- predictors", "\n")
  
  for (ii_pred in 1:length(remaining_predictors)) {
    # Iterate over the possible combinations of remaining predictors
    log_lik_log <- as.data.frame(matrix(nrow = n_chains * n_iter / 2, ncol = nrow(dat)))
    log_lik_lin <- as.data.frame(matrix(nrow = n_chains * n_iter / 2, ncol = nrow(dat)))
    
    current_predictors <- c(req_predictors, remaining_predictors[ii_pred]) # list of predictors to use
    mdl_name <- paste0("f", length(current_predictors) + 1, ".", ii_pred) # model name for future ref
    
    # Toggle on/off parameters based on inclusion
    if ("log10_dose_pfu" %in% current_predictors) {dose_pred <- 1}  else {dose_pred <- 0}
    if ("dpi_idx" %in% current_predictors) {dpi_pred <- 1}  else {dpi_pred <- 0}
    if ("sp_idx" %in% current_predictors) {sp_pred <- 1}  else {sp_pred <- 0}
    if ("st_idx" %in% current_predictors) {st_pred <- 1}  else {st_pred <- 0}
    if ("tg_idx" %in% current_predictors) {tg_pred <- 1}  else {tg_pred <- 0}
    if ("age_idx" %in% current_predictors) {age_pred <- 1}  else {age_pred <- 0}
    if ("sex_idx" %in% current_predictors) {sex_pred <- 1}  else {sex_pred <- 0}
    
    # To calculate number of correctly classified values in train/test data
    n_within_95_test <- 0
    n_within_50_test <- 0
    n_within_95_train <- 0
    n_within_50_train <- 0
    n_corr_test <- 0
    n_corr_train <- 0
    mae_test <- c()
    mae_scaled_test <- c()
    mae_train <- c()
    mae_scaled_train <- c()
    
    ### K-fold cross validation ------------------------------------------------
    
    for (k in 1:n_folds) {
      
      cat("Running fold ", k, "for predictor", remaining_predictors[ii_pred],
          "and model", mdl_name, "\n")
      cat("Note: the chosen predictors so far are:", paste0(req_predictors, collapse = ", "), "\n")
      
      # For ease of formatting data for stan
      dat.train <- subset(dat, fold != k)
      dat.test <- subset(dat, fold == k)
      
      # Formatting for stan
      dat.stan <- c(list(
        
        # Number of observations, articles, & levels 
        N_test = dim(dat.test)[1],
        N_train = dim(dat.train)[1],
        N_test_lin = nrow(subset(dat.test, cens_sg == 0)),
        N_article = max(dat.train$article_idx),
        L_sp = max(dat.train$sp_idx),
        L_st = max(dat.train$st_idx),
        L_tg = max(dat.train$tg_idx),
        L_dpi = max(dat.train$dpi_idx),
        L_age = max(dat.train$age_idx),
        
        # Toggle predictors on/off
        dose_inc = dose_pred,
        dpi_inc = dpi_pred,
        sp_inc = sp_pred,
        tg_inc = tg_pred,
        st_inc = st_pred,
        age_inc = age_pred,
        sex_inc = sex_pred,
        
        # Response variables
        sg_pos_train = dat.train$pos_sg,
        sg_pos_test = dat.test$pos_sg,
        sg_train = dat.train$val_sg,
        sg_test = dat.test$val_sg,
        
        # Predictor variables
        t_train = dat.train$val_total,
        t_test = dat.test$val_total,
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
        
        # Housekeeping for tRNA censoring
        t_cens_train = dat.train$cens_total,
        t_cens_test = dat.test$cens_total,
        
        # Article hierarchy
        article_train = dat.train$article_idx, 
        article_test = dat.test$article_idx))
      
      # Run current model
      fit.current <- model.full.xval$sample(data = dat.stan,
                                           chains = n_chains,
                                           parallel_chains = n_chains,
                                           iter_warmup = n_iter / 2,
                                           iter_sampling = n_iter / 2,
                                           refresh = 1000,
                                           seed = seed)
      
      # Extract & save log likelihood for logistic model
      log_lik_log_pars <- paste0("log_lik_log[", 1:nrow(dat.test), "]")
      log_lik_log[, dat$fold == k] <- as.data.frame(fit.current$draws(log_lik_log_pars, 
                                                                      format = "matrix"))
      
      # Extract & save log likelihood for linear model
      log_lik_lin_pars <- paste0("log_lik_lin[", 1:nrow(subset(dat.test, cens_sg == 0)), "]")
      log_lik_lin[, dat$fold == k & dat$cens_sg == 0] <- as.data.frame(fit.current$draws(log_lik_lin_pars, 
                                                                                         format = "matrix"))
      
      # TEST DATA
      
      dat.test.obs <- subset(dat.test, cens_sg == 0)
      
      # Determine number of samples within 50 and 95% prediction interval
      pred_val_test <- paste0("pred_val_test[", which(dat.test$cens_sg == 0), "]") 
      pred_test_df <- as.data.frame(fit.current$draws(pred_val_test, format = "df"))
      pred_test_df <- subset(pred_test_df, select = pred_val_test) #remove .chain, ...
      
      pred_0.025 <- sapply(pred_test_df, function(x) quantile(x, probs = 0.025))
      pred_0.25 <- sapply(pred_test_df, function(x) quantile(x, probs = 0.25))
      pred_0.75 <- sapply(pred_test_df, function(x) quantile(x, probs = 0.75))
      pred_0.975 <- sapply(pred_test_df, function(x) quantile(x, probs = 0.975))
      
      n_within_95 <- sum(pred_0.025 <= dat.test.obs$val_sg & dat.test.obs$val_sg <= pred_0.975)
      n_within_50 <- sum(pred_0.25 <= dat.test.obs$val_sg & dat.test.obs$val_sg <= pred_0.75)
      n_within_95_test <- n_within_95_test + n_within_95 
      n_within_50_test <- n_within_50_test + n_within_50 
      
      # Calculate MAE and MAE scaled
      pred_median <- sapply(pred_test_df, function(x) median(x))
      pred_sd <- sapply(pred_test_df, function(x) sd(x))
      mae_test <- c(mae_test, abs(dat.test.obs$val_sg - pred_median))
      mae_scaled_test <- c(mae_scaled_test, abs(dat.test.obs$val_sg - pred_median)/pred_sd)
      
      
      # TRAIN DATA
      
      dat.train.obs <- subset(dat.train, cens_sg == 0)
      
      # Determine number of samples within 50 and 95% prediction interval
      pred_val_train <- paste0("pred_val_train[", which(dat.train$cens_sg == 0), "]") 
      pred_train_df <- as.data.frame(fit.current$draws(pred_val_train, format = "df"))
      pred_train_df <- subset(pred_train_df, select = pred_val_train) #remove .chain, ...
      
      pred_0.025 <- sapply(pred_train_df, function(x) quantile(x, probs = 0.025))
      pred_0.25 <- sapply(pred_train_df, function(x) quantile(x, probs = 0.25))
      pred_0.75 <- sapply(pred_train_df, function(x) quantile(x, probs = 0.75))
      pred_0.975 <- sapply(pred_train_df, function(x) quantile(x, probs = 0.975))
      
      n_within_95 <- sum(pred_0.025 <= dat.train.obs$val_sg & dat.train.obs$val_sg <= pred_0.975)
      n_within_50 <- sum(pred_0.25 <= dat.train.obs$val_sg & dat.train.obs$val_sg <= pred_0.75)
      n_within_95_train <- n_within_95_train + n_within_95 
      n_within_50_train <- n_within_50_train + n_within_50 
      
      # Add MAE and MAE scaled terms 
      pred_median <- sapply(pred_train_df, function(x) median(x))
      pred_sd <- sapply(pred_train_df, function(x) sd(x))
      mae_train <- c(mae_train, abs(dat.train.obs$val_sg - pred_median))
      mae_scaled_train <- c(mae_scaled_train, abs(dat.train.obs$val_sg - pred_median)/pred_sd)
    
    }
    
    # Combine linear & log likelihoods
    log_lik <- cbind(log_lik_log, subset(log_lik_lin, select = c(which(!is.na(log_lik_lin[1, ])))))
    
    # Generate ELPD & similar statistics based on equations in Vehtari et al. 2017
    elpd_current <- sum(log(colMeans(exp(as.matrix(log_lik)))))
    elpd_se_current <- sqrt(ncol(log_lik)) * sd(log(colMeans(exp(as.matrix(log_lik)))))
    looic_current <- -2 * elpd_current 
    looic_se_current <- 2 * elpd_se_current
    ploo_current <- sum(colVars(as.matrix(log_lik)))
    ploo_se_current <- sqrt(ncol(log_lik)) * sd(colVars(as.matrix(log_lik)))
    
    # Generate MAE statistics
    mae_train_full <- median(mae_train)
    mae_scaled_train_full <- median(mae_scaled_train)
    mae_test_full <- median(mae_test)
    mae_scaled_test_full <- median(mae_scaled_test)
    
    # Store log_lik dataframe for future reference
    loglik_name <- paste0("log_lik_", mdl_name) # x-val likelihood matrix
    assign(loglik_name, log_lik)
    
    
    ## PSIS-LOO method ---------------------------------------------------------
    
    cat("Now running PSIS-LOO method for predictor", remaining_predictors[ii_pred],
        "and model", mdl_name, "\n")
    
    # Formatting for stan
    dat.stan <- c(list(
      
      # Number of observations, articles, & levels 
      N = nrow(dat),
      N_lin = nrow(subset(dat, cens_sg == 0)),
      N_article = max(dat$article_idx),
      L_sp = max(dat$sp_idx),
      L_st = max(dat$st_idx),
      L_tg = max(dat$tg_idx),
      L_dpi = max(dat$dpi_idx),
      L_age = max(dat$age_idx),
      
      # Toggle predictors on/off
      dose_inc = dose_pred,
      dpi_inc = dpi_pred,
      sp_inc = sp_pred,
      tg_inc = tg_pred,
      st_inc = st_pred,
      age_inc = age_pred,
      sex_inc = sex_pred,
      
      # Response variables
      sg_pos = dat$pos_sg,
      sg = dat$val_sg,
      
      # Predictor variables
      t = dat$val_total,
      dose = dat$log10_dose_pfu,
      tg = dat$tg_idx,
      dpi = dat$dpi_idx,
      age = dat$age_idx,
      sex = dat$sex_idx,
      st = dat$st_idx,
      sp = dat$sp_idx,
      
      # Housekeeping for tRNA censoring
      t_cens = dat$cens_total,
      
      # Article hierarchy
      article = dat$article_idx))
    
    # Run the model with rstan
    fit.current <- sampling(model.full.loo, data = dat.stan, 
                            iter = n_iter, chains = n_chains,
                            seed = seed)
    loo_current <- loo(fit.current, moment_match = TRUE)
    
    # Generate ELPD & other statistics with LOO package, for comparison
    elpd_loo_current <- loo_current$estimates[1, 1]
    elpd_loo_se_current <- loo_current$estimates[1, 2]
    ploo_loo_current <- loo_current$estimates[2, 1]
    ploo_loo_se_current <- loo_current$estimates[2, 2]
    pareto_k_max <- max(pareto_k_values(loo_current))
    
    # Append statistics to performance table
    next_row <- c(list(model = mdl_name,
                       predictors = paste0(current_predictors, collapse = ", "),
                       elpd = elpd_current,
                       elpd_se = elpd_se_current,
                       p_loo = ploo_current,
                       p_loo_se = ploo_se_current,
                       looic = looic_current,
                       looic_se = looic_se_current,
                       elpd_loo = elpd_loo_current,
                       elpd_loo_se = elpd_loo_se_current,
                       p_loo_loo = ploo_loo_current,
                       p_loo_loo_se = ploo_loo_se_current,
                       pareto_k_max = pareto_k_max,
                       mae_test = mae_test_full,
                       mae_scaled_test = mae_scaled_test_full,
                       mae_train = mae_train_full,
                       mae_scaled_train = mae_scaled_train_full,
                       n_within_50_test = n_within_50_test,
                       n_within_50_train = n_within_50_train,
                       n_within_95_test = n_within_95_test,
                       n_within_95_train = n_within_95_train))
    tblS3 <- rbind(tblS3, next_row)
    write.csv(tblS3, file = "./outputs/tables/tblS3-sgRNA-full-model-selection-noninf-raw.csv",
              row.names = FALSE)
    View(tblS3)
    
    # Save PSIS-LOO loo object for future reference
    loo_name <- paste0("loo_", mdl_name) # model name for future ref
    assign(loo_name, loo_current)
    
  }
  
  # Find the best predictor for this step, and set as required moving forward
  npred.subset <- subset(tblS3, str_detect(tblS3$model, paste0("f", length(current_predictors) + 1, ".")))
  best_predictors <- npred.subset$predictors[which(npred.subset$elpd == max(npred.subset$elpd))] 
  req_predictors <- unique(c(req_predictors, unlist(strsplit(best_predictors,", "))))
  print(req_predictors)
  cat("The chosen predictors so far are:", paste0(req_predictors, collapse = ", "), "\n")
  
}


## Generate table S3 for diagnostics -------------------------------------------

tblS3.raw <- tblS3



### Add in elpd_diff -----------------------------------------------------------

best_model <- tblS3$model[which(tblS3$elpd == max(tblS3$elpd))]
best_elpd <- tblS3$elpd[tblS3$model == best_model]
best_log_lik <- get(paste0("log_lik_", best_model))

tblS3$elpd_diff <- 0
tblS3$elpd_diff_se <- 0

for (row_num in 1:nrow(tblS3)) {
  tblS3$elpd_diff[row_num] <- tblS3$elpd[row_num] - best_elpd
  current_log_lik <- get(paste0("log_lik_", tblS3$model[row_num]))
  tblS3$elpd_diff_se[row_num] <- sqrt(nrow(dat)) * sd(log(colMeans(exp(as.matrix(best_log_lik)))) -
                                                        log(colMeans(exp(as.matrix(current_log_lik)))))
}


### Add in loo-based elpd-diff -------------------------------------------------

comp_loo <- loo_compare(x = list("f1" = loo_f1,
                             "f2.1" = loo_f2.1, "f2.2" = loo_f2.2, "f2.3" = loo_f2.3,
                             "f2.4" = loo_f2.4, "f2.5" = loo_f2.5, "f2.6" = loo_f2.6,
                             "f2.7" = loo_f2.7, 
                             "f3.1" = loo_f3.1, "f3.2" = loo_f3.2, "f3.3" = loo_f3.3,
                             "f3.4" = loo_f3.4, "f3.5" = loo_f3.5, "f3.6" = loo_f3.6,
                             "f4.1" = loo_f4.1, "f4.2" = loo_f4.2, "f4.3" = loo_f4.3,
                             "f4.4" = loo_f4.4, "f4.5" = loo_f4.5, 
                             "f5.1" = loo_f5.1, "f5.2" = loo_f5.2, "f5.3" = loo_f5.3,
                             "f5.4" = loo_f5.4,
                             "f6.1" = loo_f6.1, "f6.2" = loo_f6.2, "f6.3" = loo_f6.3,
                             "f7.1" = loo_f7.1, "f7.2" = loo_f7.2,
                             "f8.1" = loo_f8.1
                             ))

tblS3$elpd_diff_loo <- 0
tblS3$elpd_diff_se_loo <- 0

for (row_num in 1:nrow(comp_loo)) {
  mdl <- rownames(comp_loo)[row_num]
  tblS3$elpd_diff_loo[tblS3$model == mdl] <- comp_loo[row_num, 1]
  tblS3$elpd_diff_se_loo[tblS3$model == mdl] <- comp_loo[row_num, 2]
}


### Compare to "best" model ----------------------------------------------------

# To add rows comparing ELPD between the "best" model
best_model <- "f5.1"
best_elpd <- tblS3$elpd[tblS3$model == best_model]
best_log_lik <- get(paste0("log_lik_", best_model))

tblS3$elpd_diff_best <- 0
tblS3$elpd_diff_best_se <- 0

for (row_num in 1:nrow(tblS3)) {
  tblS3$elpd_diff_best[row_num] <- tblS3$elpd[row_num] - best_elpd
  current_log_lik <- get(paste0("log_lik_", tblS3$model[row_num]))
  tblS3$elpd_diff_best_se[row_num] <- sqrt(nrow(dat)) * sd(log(colMeans(exp(as.matrix(best_log_lik)))) -
                                                             log(colMeans(exp(as.matrix(current_log_lik)))))
}

## Simple ELPD Diff statistic for LOO comparisons of key models
comp_loo_best <- loo_compare(x = list("f1" = loo_f1,
                                      "f5.1" = loo_f5.1, 
                                      "f8.1" = loo_f8.1))
write.csv(comp_loo_best, file = "./outputs/tables/EA-table-sgRNA-full-noninf-key-model-comparison-loo.csv",
          row.names = TRUE)



### Generate percentages for correctly predicted samples -----------------------

dat.obs <- subset(dat, cens_sg == 0)
tblS3$percent_50_test <- tblS3$n_within_50_test / nrow(dat.obs) * 100
tblS3$percent_50_train <- tblS3$n_within_50_train / (nrow(dat.obs) * (n_folds - 1)) * 100
tblS3$percent_95_test <- tblS3$n_within_95_test / nrow(dat.obs) * 100
tblS3$percent_95_train <- tblS3$n_within_95_train / (nrow(dat.obs) * (n_folds - 1)) * 100


### Update predictor names -----------------------------------------------------

tblS3$predictors <- str_remove_all(tblS3$predictors, "_idx")
tblS3$predictors <- str_remove_all(tblS3$predictors, "log10_|_pfu")
tblS3$predictors <- str_replace_all(tblS3$predictors, ",", " +" )
tblS3$predictors <- toupper(tblS3$predictors)


### Save raw version for fig2 --------------------------------------------------

write.csv(tblS3, file = "./outputs/tables/tblS3-sgRNA-full-model-selection-noninf-raw.csv",
          row.names = FALSE)


### Aesthetic changes ----------------------------------------------------------

# Rounding: prediction statistics
tblS3$percent_50_test <- round(tblS3$percent_50_test, digits = 1)
tblS3$percent_50_train <- round(tblS3$percent_50_train , digits = 1)
tblS3$percent_95_test <- round(tblS3$percent_95_test, digits = 1)
tblS3$percent_95_train <- round(tblS3$percent_95_train, digits = 1)
tblS3$mae_test <- round(tblS3$mae_test, digits = 2)
tblS3$mae_train <- round(tblS3$mae_train, digits = 2)
tblS3$mae_scaled_test <- round(tblS3$mae_scaled_test, digits = 2)
tblS3$mae_scaled_train <- round(tblS3$mae_scaled_train, digits = 2)

# Rounding: x-val statistics
tblS3$elpd <- round(tblS3$elpd, digits = 2)
tblS3$elpd_se <- round(tblS3$elpd_se, digits = 2)
tblS3$elpd_diff <- round(tblS3$elpd_diff, digits = 2)
tblS3$elpd_diff_se <- round(tblS3$elpd_diff_se, digits = 2)
tblS3$elpd_diff_best <- round(tblS3$elpd_diff_best, digits = 2)
tblS3$elpd_diff_best_se <- round(tblS3$elpd_diff_best_se, digits = 2)

# Rounding: LOO statistics
tblS3$elpd_diff_loo <- round(tblS3$elpd_diff_loo, digits = 2)
tblS3$elpd_diff_se_loo <- round(tblS3$elpd_diff_se_loo, digits = 2)
tblS3$elpd_loo <- round(tblS3$elpd_loo, digits = 2)
tblS3$elpd_loo_se <- round(tblS3$elpd_loo_se, digits = 2)
tblS3$p_loo_loo <- round(tblS3$p_loo_loo, digits = 2)
tblS3$p_loo_loo_se <- round(tblS3$p_loo_loo_se, digits = 2)


# Generate consolidated columns 
tblS3$"ELPD (SE)" <- paste0(round(tblS3$elpd, digits = 2),
                            " (", 
                            round(tblS3$elpd_se, digits = 2), 
                            ")")
tblS3$"ELPD Difference (SE)" <- paste0(round(tblS3$elpd_diff, digits = 2),
                                       " (", 
                                       round(tblS3$elpd_diff_se, digits = 2), 
                                       ")")
tblS3$"ELPD Difference Best (SE)" <- paste0(round(tblS3$elpd_diff_best, digits = 2),
                                            " (", 
                                            round(tblS3$elpd_diff_best_se, digits = 2), 
                                            ")")
tblS3$"LOOIC (SE)" <- paste0(round(tblS3$looic, digits = 2),
                             " (", 
                             round(tblS3$looic_se, digits = 2), 
                             ")")
tblS3$"Effective parameters (SE)" <- paste0(round(tblS3$p_loo, digits = 2),
                                            " (", 
                                            round(tblS3$p_loo_se, digits = 2), 
                                            ")")


tblS3$"LOO - ELPD (SE)" <- paste0(round(tblS3$elpd_loo, digits = 2),
                            " (", 
                            round(tblS3$elpd_loo_se, digits = 2), 
                            ")")
tblS3$"LOO - ELPD Difference (SE)" <- paste0(round(tblS3$elpd_diff_loo, digits = 2),
                                       " (", 
                                       round(tblS3$elpd_diff_se_loo, digits = 2), 
                                       ")")
tblS3$"LOO - Effective parameters (SE)" <- paste0(round(tblS3$p_loo_loo, digits = 2),
                                            " (", 
                                            round(tblS3$p_loo_loo_se, digits = 2), 
                                            ")")

tblS3$"MAE (Scaled) - Train" <- paste0(tblS3$mae_train,
                                       " (", 
                                       tblS3$mae_scaled_train, 
                                       ")")
tblS3$"MAE (Scaled) - Test" <- paste0(tblS3$mae_test,
                                       " (", 
                                       tblS3$mae_scaled_test, 
                                       ")")


# Rename columns to be more informative
colnames(tblS3)[colnames(tblS3) == "percent_50_test"] <- "% within 50% PI (test)"
colnames(tblS3)[colnames(tblS3) == "percent_50_train"] <- "% within 50% PI (train)"
colnames(tblS3)[colnames(tblS3) == "percent_95_test"] <- "% within 95% PI (test)"
colnames(tblS3)[colnames(tblS3) == "percent_95_train"] <- "% within 95% PI (train)"
colnames(tblS3)[colnames(tblS3) == "pareto_k_max"] <- "LOO - Max. Pareto k"
colnames(tblS3)[colnames(tblS3) == "model"] <- "Model"
colnames(tblS3)[colnames(tblS3) == "predictors"] <- "Predictors"

# Order by decreasing elpd value
tblS3 <- tblS3[order(tblS3$Model), ]

# Select only the columns of interest, and in the correct order
tblS3 <- tblS3 %>% select(c("Model", "Predictors", 
                            "ELPD Difference (SE)", "ELPD (SE)", 
                            "ELPD Difference Best (SE)",
                            "Effective parameters (SE)",
                            "LOO - ELPD Difference (SE)", 
                            "LOO - ELPD (SE)",
                            "LOO - Effective parameters (SE)",
                            "LOO - Max. Pareto k",
                            "MAE (Scaled) - Train", 
                            "MAE (Scaled) - Test", 
                            "% within 50% PI (train)", 
                            "% within 50% PI (test)", 
                            "% within 95% PI (train)",
                            "% within 95% PI (test)"))


## Save table ------------------------------------------------------------------

write.csv(tblS3, file = "./outputs/tables/tblS3-sgRNA-full-model-selection-noninf-formatted.csv",
          row.names = FALSE)

