# This runs cross-validation on the considered culture base models
#   including totRNA, sgRNA, or both as primary predictors


# Prep environment -------------------------------------------------------------

## Install CmdStanR, CmdStan, and RStan for model compilation ------------------

#  We recommend you review the installation instructions at:
#     https://mc-stan.org/cmdstanr/articles/cmdstanr.html 
#     https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started 

library("cmdstanr")
library("rstan")


## Install other packages ------------------------------------------------------

req_pkgs <- c("loo", "stringr", "tidyverse", "ggridges", "bayesplot",
              "gtools", "wesanderson", "cowplot", "ggpubr", "matrixStats",
              "ggExtra")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependences = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))

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

# Every model should generate predictions for all entries in this df
dat.inf <- subset(dat, !is.na(pos_inf)) 

# Values need to be numerics
dat.inf$val_sg <- as.numeric(dat.inf$val_sg)
dat.inf$val_total <- as.numeric(dat.inf$val_total)

# Update columns for Stan's sake, although not used for base models
dat.inf$article_idx <-  as.numeric(factor(dat.inf$article))


# Load model ------------------------------------------------------------------

model.culture.xval <- cmdstan_model('./code/stan-model-culture-xval.stan')


# Prep results table -----------------------------------------------------------

tbl_cols <- c("model", "predictors", "elpd", "elpd_se", "elpd_loo", "elpd_loo_se", 
              "p_loo_loo", "p_loo_loo_se", "pareto_k_max", 
              "F_score", "MCC",
              "n_correct_test", "n_correct_train", 
              "n_total_test", "n_total_train",
              "percent_pos", "percent_neg", "percent_total")

tblSX <- as.data.frame(matrix(nrow = 0, ncol = length(tbl_cols)))
colnames(tblSX) <- tbl_cols


# NO PREDICTIONS ---------------------------------------------------------------

# Every model should generate predictions for all entries in this df
dat.no.pred <- subset(dat.inf, pos_sg %in% c(0, 1) & pos_total %in% c(0, 1) & cens_sg != 3) 


## Assign folds ----------------------------------------------------------------

n_folds <- 10
set.seed(16669)
dat.no.pred$fold <- kfold_split_stratified(K = n_folds, x = dat.no.pred$article_idx)


## Evaluate RNA negative samples ------------------------------------------------

dat.neg <- subset(dat.no.pred, cens_total == 1 | pos_sg == 0)

# total RNA 
TN_T_neg <- nrow(subset(dat.neg, pos_inf == 0 & cens_total == 1))
FN_T_neg <- nrow(subset(dat.neg, pos_inf == 1 & cens_total == 1))

# sgRNA
TN_SG_neg <- nrow(subset(dat.neg, pos_inf == 0 & pos_sg == 0))
FN_SG_neg <- nrow(subset(dat.neg, pos_inf == 1 & pos_sg == 0))

# both together
TN_Both_neg <- nrow(subset(dat.neg, pos_inf == 0 & (pos_sg == 0 | cens_total == 1)))
FN_Both_neg <- nrow(subset(dat.neg, pos_inf == 1 & (pos_sg == 0 | cens_total == 1)))


## Total RNA ------------------------------------------------------------------

# Subset to total RNA positive data
dat.total <- subset(dat.no.pred, cens_total == 0)
log_lik <- as.data.frame(matrix(nrow = n_chains * n_iter / 2, ncol = nrow(dat.total)))


# Prep df to store predictions
pred.df <- list(true_val = c(),
                pred_val = c())

  
#### K-fold cross validation ---------------------------------------------------

n_corr_test <- 0
n_corr_train <- 0

for (k in 1:n_folds) {
  
  cat("Running fold ", k, "for the quantitative base model: T \n \n")
  
  # For ease of formatting data for stan
  dat.train <- subset(dat.total, fold != k)
  dat.test <- subset(dat.total, fold == k)
  
  # Formatting for stan
  dat.stan <- c(list(
    
    # Number of observations, articles, & levels 
    N_test = dim(dat.test)[1],
    N_train = dim(dat.train)[1],
    N_article = max(dat.train$article_idx),
    L_sp = 3,
    L_tg = 4,
    L_dpi = 3,
    L_age = 3,
    L_cell = 3,
    
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
    sg_train = dat.train$val_sg,
    sg_test = dat.test$val_sg,
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
  log_lik[, dat.total$fold == k] <- as.data.frame(fit.current$draws(log_lik_pars, 
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
loglik_name <- paste0("log_lik_", "q1.1") # model name for future ref
assign(loglik_name, log_lik)


#### PSIS-LOO method ------------------------------------------------------------

dat.train <- dat.total
dat.test <- dat.total

# Formatting for stan
dat.stan <- c(list(
  
  # Number of observations, articles, & levels 
  N_test = dim(dat.test)[1],
  N_train = dim(dat.train)[1],
  N_article = max(dat.train$article_idx),
  L_sp = 3,
  L_tg = 4,
  L_dpi = 3,
  L_age = 3,
  L_cell = 3,
  
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
  sg_train = dat.train$val_sg,
  sg_test = dat.test$val_sg,
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
log_lik_pars <- paste0("log_lik[", 1:nrow(dat.total), "]")
log_lik <- as.data.frame(fit.current$draws(log_lik_pars,format = "matrix"))

# Use loo package to extract elpd, p_loo, and pareto_k statistics
loo_current <- loo(as.matrix(log_lik), r_eff = NA)
elpd_loo_current <- loo_current$estimates[1, 1]
elpd_loo_se_current <- loo_current$estimates[1, 2]
ploo_current <- loo_current$estimates[2, 1]
ploo_se_current <- loo_current$estimates[2, 2]
pareto_k_max <- max(pareto_k_values(loo_current))

# Append statistics to performance table
next_row <- c(list(model = "q1.1",
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
                   n_correct_test = TP + TN,
                   n_correct_train = n_corr_train,
                   n_total_test = TP + TN + FP + FN,
                   n_total_train = nrow(dat.total) * (n_folds - 1),
                   percent_pos = percent_pos,
                   percent_neg = percent_neg,
                   percent_total = percent_total))
tblSX <- rbind(tblSX, next_row)
View(tblSX)

# Convert log likelihoods to loo object for future reference
loo_name <- paste0("loo_", "q1.1") # model name for future ref
assign(loo_name, loo_current)


## sgRNA ----------------------------------------------------------------------

# Subset to sgRNA positive data
dat.sg <- subset(dat.no.pred, pos_sg == 1)
log_lik <- as.data.frame(matrix(nrow = n_chains * n_iter / 2, ncol = nrow(dat.sg)))

# Prep df to store predictions
pred.df <- list(true_val = c(),
                pred_val = c())


#### K-fold cross validation ---------------------------------------------------

n_corr_test <- 0
n_corr_train <- 0

for (k in 1:n_folds) {
  
  cat("Running fold ", k, "for the quantitative base model: SG \n \n")
  
  # For ease of formatting data for stan
  dat.train <- subset(dat.sg, fold != k)
  dat.test <- subset(dat.sg, fold == k)
  
  # Formatting for stan
  dat.stan <- c(list(
    
    # Number of observations, articles, & levels 
    N_test = dim(dat.test)[1],
    N_train = dim(dat.train)[1],
    N_article = max(dat.train$article_idx),
    L_sp = 3,
    L_tg = 4,
    L_dpi = 3,
    L_age = 3,
    L_cell = 3,
    
    # Toggle predictors on/off
    t_inc = 0,
    sg_inc = 1,
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
    sg_train = dat.train$val_sg,
    sg_test = dat.test$val_sg,
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
  log_lik[, dat.sg$fold == k] <- as.data.frame(fit.current$draws(log_lik_pars, 
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
TN <- sum(pred.df$true_val == 0 & pred.df$pred_val == 0) + TN_SG_neg
FN <- sum(pred.df$true_val == 1 & pred.df$pred_val == 0) + FN_SG_neg
F_score <- (2 * TP) / (2 * TP + FP + FN)
MCC <- ((TP * TN) - (FP * FN))/(sqrt(TP + FP)*sqrt(TP + FN)*sqrt(TN + FP)*sqrt(TN + FN))

# Percent of samples correctly classified, stratified by Pos/Neg/All 
percent_pos <- TP / (TP + FN)  * 100
percent_neg <- TN / (TN + FP) * 100
percent_total <- (TN + TP) / (TP + FP + TN + FN) * 100

# Store log_lik data frame for future reference
loglik_name <- paste0("log_lik_", "c1.2") # model name for future ref
assign(loglik_name, log_lik)


#### PSIS-LOO method ------------------------------------------------------------

dat.train <- dat.sg
dat.test <- dat.sg

# Formatting for stan
dat.stan <- c(list(
  
  # Number of observations, articles, & levels 
  N_test = dim(dat.test)[1],
  N_train = dim(dat.train)[1],
  N_article = max(dat.train$article_idx),
  L_sp = 3,
  L_tg = 4,
  L_dpi = 3,
  L_age = 3,
  L_cell = 3,
  
  # Toggle predictors on/off
  t_inc = 0,
  sg_inc = 1,
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
  sg_train = dat.train$val_sg,
  sg_test = dat.test$val_sg,
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

fit.sg.quant <- fit.current

# Store log likelihood
log_lik_pars <- paste0("log_lik[", 1:nrow(dat.sg), "]")
log_lik <- as.data.frame(fit.current$draws(log_lik_pars,format = "matrix"))

# Use loo package to extract elpd, ic, p_loo, MC_se, and pareto_k statistics
loo_current <- loo(as.matrix(log_lik), r_eff = NA)
elpd_loo_current <- loo_current$estimates[1, 1]
elpd_loo_se_current <- loo_current$estimates[1, 2]
ploo_current <- loo_current$estimates[2, 1]
ploo_se_current <- loo_current$estimates[2, 2]
pareto_k_max <- max(pareto_k_values(loo_current))


# Append statistics to performance table
next_row <- c(list(model = "q1.2",
                   predictors = "SG",
                   elpd = elpd_current,
                   elpd_se = elpd_se_current,
                   elpd_loo = elpd_loo_current,
                   elpd_loo_se = elpd_loo_se_current,
                   p_loo_loo = ploo_current,
                   p_loo_loo_se = ploo_se_current,
                   pareto_k_max = pareto_k_max,
                   F_score = F_score,
                   MCC = MCC,
                   n_correct_test = TP + TN,
                   n_correct_train = n_corr_train,
                   n_total_test = TP + TN + FP + FN,
                   n_total_train = nrow(dat.sg) * (n_folds - 1),
                   percent_pos = percent_pos,
                   percent_neg = percent_neg,
                   percent_total = percent_total))
tblSX <- rbind(tblSX, next_row)

# Convert log likelihoods to loo object for future reference
loo_name <- paste0("loo_", "c1.2") # model name for future ref
assign(loo_name, loo_current)


## total RNA + sgRNA --------------------------------------------------------------

dat.both <- subset(dat.no.pred, pos_sg == 1 & cens_total == 0)
log_lik <- as.data.frame(matrix(nrow = n_chains * n_iter / 2, ncol = nrow(dat.both)))

# Prep df to store predictions
pred.df <- list(true_val = c(),
                pred_val = c())


#### K-fold cross validation ---------------------------------------------------

n_corr_test <- 0
n_corr_train <- 0

for (k in 1:n_folds) {
  
  cat("Running fold ", k, "for the quantitative base model: T + SG \n \n")
  
  # For ease of formatting data for stan
  dat.train <- subset(dat.both, fold != k)
  dat.test <- subset(dat.both, fold == k)
  
  # Formatting for stan
  dat.stan <- c(list(
    
    # Number of observations, articles, & levels 
    N_test = dim(dat.test)[1],
    N_train = dim(dat.train)[1],
    N_article = max(dat.train$article_idx),
    L_sp = 3,
    L_tg = 4,
    L_dpi = 3,
    L_age = 3,
    L_cell = 3,
    
    # Toggle predictors on/off
    t_inc = 1,
    sg_inc = 1,
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
    sg_train = dat.train$val_sg,
    sg_test = dat.test$val_sg,
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
  log_lik[, dat.both$fold == k] <- as.data.frame(fit.current$draws(log_lik_pars, 
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
TN <- sum(pred.df$true_val == 0 & pred.df$pred_val == 0) + TN_Both_neg
FN <- sum(pred.df$true_val == 1 & pred.df$pred_val == 0) + FN_Both_neg
F_score <- (2 * TP) / (2 * TP + FP + FN)
MCC <- ((TP * TN) - (FP * FN))/(sqrt(TP + FP)*sqrt(TP + FN)*sqrt(TN + FP)*sqrt(TN + FN))

# Percent of samples correctly classified, stratified by Pos/Neg/All 
percent_pos <- TP / (TP + FN)  * 100
percent_neg <- TN / (TN + FP) * 100
percent_total <- (TN + TP) / (TP + FP + TN + FN) * 100

# Store log_lik data frame for future reference
loglik_name <- paste0("log_lik_", "c1.3") # model name for future ref
assign(loglik_name, log_lik)


#### PSIS-LOO method ------------------------------------------------------------

# Will run slower, but prevents needing an additional model file
dat.train <- dat.both
dat.test <- dat.both

# Formatting for stan
dat.stan <- c(list(
  
  # Number of observations, articles, & levels 
  N_test = dim(dat.test)[1],
  N_train = dim(dat.train)[1],
  N_article = max(dat.train$article_idx),
  L_sp = 3,
  L_tg = 4,
  L_dpi = 3,
  L_age = 3,
  L_cell = 3,
  
  # Toggle predictors on/off
  t_inc = 1,
  sg_inc = 1,
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
  sg_train = dat.train$val_sg,
  sg_test = dat.test$val_sg,
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

fit.both.quant <- fit.current

# Store log likelihood
log_lik_pars <- paste0("log_lik[", 1:nrow(dat.both), "]")
log_lik <- as.data.frame(fit.current$draws(log_lik_pars,format = "matrix"))

# Use loo package to extract elpd, p_loo, and pareto_k statistics
loo_current <- loo(as.matrix(log_lik), r_eff = NA)
elpd_loo_current <- loo_current$estimates[1, 1]
elpd_loo_se_current <- loo_current$estimates[1, 2]
ploo_current <- loo_current$estimates[2, 1]
ploo_se_current <- loo_current$estimates[2, 2]
pareto_k_max <- max(pareto_k_values(loo_current))

# Append statistics to performance table
next_row <- c(list(model = "q1.3",
                   predictors = "T + SG",
                   elpd = elpd_current,
                   elpd_se = elpd_se_current,
                   elpd_loo = elpd_loo_current,
                   elpd_loo_se = elpd_loo_se_current,
                   p_loo_loo = ploo_current,
                   p_loo_loo_se = ploo_se_current,
                   pareto_k_max = pareto_k_max,
                   F_score = F_score,
                   MCC = MCC,
                   n_correct_test = TP + TN,
                   n_correct_train = n_corr_train,
                   n_total_test = TP + TN + FP + FN,
                   n_total_train = nrow(dat.both) * (n_folds - 1),
                   percent_pos = percent_pos,
                   percent_neg = percent_neg,
                   percent_total = percent_total))
tblSX <- rbind(tblSX, next_row)


# Convert log likelihoods to loo object for future reference
loo_name <- paste0("loo_", "c1.3") # model name for future ref
assign(loo_name, loo_current)


# WITH PREDICTIONS -------------------------------------------------------------

# Not including the samples without totRNA information
dat.pred <- subset(dat.inf, pos_total %in% c(0, 1) & pos_sg_combined %in% c(0, 1)) 

## Assign folds ----------------------------------------------------------------

n_folds <- 10
set.seed(16669)
dat.pred$fold <- kfold_split_stratified(K = n_folds, x = dat.pred$article_idx)


## Evaluate RNA negative samples ------------------------------------------------

dat.neg <- subset(dat.pred, cens_total == 1 | pos_sg_combined == 0)

# total RNA 
TN_T_neg <- nrow(subset(dat.neg, pos_inf == 0 & cens_total == 1))
FN_T_neg <- nrow(subset(dat.neg, pos_inf == 1 & cens_total == 1))

# sgRNA
TN_SG_neg <- nrow(subset(dat.neg, pos_inf == 0 & pos_sg_combined == 0))
FN_SG_neg <- nrow(subset(dat.neg, pos_inf == 1 & pos_sg_combined == 0))

# both together
TN_Both_neg <- nrow(subset(dat.neg, pos_inf == 0 & (pos_sg_combined == 0 | cens_total == 1)))
FN_Both_neg <- nrow(subset(dat.neg, pos_inf == 1 & (pos_sg_combined == 0 | cens_total == 1)))


## Total RNA ------------------------------------------------------------------

# Subset to total RNA positive data
dat.total <- subset(dat.pred, cens_total == 0)
log_lik <- as.data.frame(matrix(nrow = n_chains * n_iter / 2, ncol = nrow(dat.total)))


# Prep df to store predictions
pred.df <- list(true_val = c(),
                pred_val = c())


### K-fold cross validation ---------------------------------------------------

n_corr_test <- 0
n_corr_train <- 0

for (k in 1:n_folds) {
  
  cat("Running fold ", k, "for the quantitative base model: T \n \n")
  
  # For ease of formatting data for stan
  dat.train <- subset(dat.total, fold != k)
  dat.test <- subset(dat.total, fold == k)
  
  # Formatting for stan
  dat.stan <- c(list(
    
    # Number of observations, articles, & levels 
    N_test = dim(dat.test)[1],
    N_train = dim(dat.train)[1],
    N_article = max(dat.train$article_idx),
    L_sp = 3,
    L_tg = 4,
    L_dpi = 3,
    L_age = 3,
    L_cell = 3,
    
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
    sg_train = dat.train$val_sg_combined,
    sg_test = dat.test$val_sg_combined,
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
  log_lik[, dat.total$fold == k] <- as.data.frame(fit.current$draws(log_lik_pars, 
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
loglik_name <- paste0("log_lik_", "q1.1") # model name for future ref
assign(loglik_name, log_lik)


### PSIS-LOO method ------------------------------------------------------------

# Will run slightly slower, but prevents needing an additional model file
dat.train <- dat.total
dat.test <- dat.total

# Formatting for stan
dat.stan <- c(list(
  
  # Number of observations, articles, & levels 
  N_test = dim(dat.test)[1],
  N_train = dim(dat.train)[1],
  N_article = max(dat.train$article_idx),
  L_sp = 3,
  L_tg = 4,
  L_dpi = 3,
  L_age = 3,
  L_cell = 3,
  
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
  sg_train = dat.train$val_sg_combined,
  sg_test = dat.test$val_sg_combined,
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
log_lik_pars <- paste0("log_lik[", 1:nrow(dat.total), "]")
log_lik <- as.data.frame(fit.current$draws(log_lik_pars,format = "matrix"))

# Use loo package to extract elpd, ic, p_loo, MC_se, and pareto_k statistics
loo_current <- loo(as.matrix(log_lik), r_eff = NA)
elpd_loo_current <- loo_current$estimates[1, 1]
elpd_loo_se_current <- loo_current$estimates[1, 2]
ploo_current <- loo_current$estimates[2, 1]
ploo_se_current <- loo_current$estimates[2, 2]
pareto_k_max <- max(pareto_k_values(loo_current))

# Append statistics to performance table
next_row <- c(list(model = "q1.1",
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
                   n_correct_test = TP + TN,
                   n_correct_train = n_corr_train,
                   n_total_test = TP + TN + FP + FN,
                   n_total_train = nrow(dat.total) * (n_folds - 1),
                   percent_pos = percent_pos,
                   percent_neg = percent_neg,
                   percent_total = percent_total))
tblSX <- rbind(tblSX, next_row)
View(tblSX)

# Convert log likelihoods to loo object for future reference
loo_name <- paste0("loo_", "q1.1") # model name for future ref
assign(loo_name, loo_current)



## sgRNA ----------------------------------------------------------------------

# Subset to sgRNA positive data, including data with sgRNA & culture results
dat.sg <- subset(dat.inf,  pos_sg_combined == 1)
dat.sg$fold <- kfold_split_stratified(K = n_folds, x = dat.sg$article_idx)
dat.sg$val_total[is.na(dat.sg$val_total)] <- -999
dat.sg$tg_idx[is.na(dat.sg$tg_idx)] <- 3

TN_SG_neg <- nrow(subset(dat.inf, pos_inf == 0 & pos_sg_combined == 0))
FN_SG_neg <- nrow(subset(dat.inf, pos_inf == 1 & pos_sg_combined == 0))

log_lik <- as.data.frame(matrix(nrow = n_chains * n_iter / 2, ncol = nrow(dat.sg)))

# Prep df to store predictions
pred.df <- list(true_val = c(),
                pred_val = c())


#### K-fold cross validation ---------------------------------------------------

n_corr_test <- 0
n_corr_train <- 0

for (k in 1:n_folds) {
  
  cat("Running fold ", k, "for the quantitative base model: SG \n \n")
  
  # For ease of formatting data for stan
  dat.train <- subset(dat.sg, fold != k)
  dat.test <- subset(dat.sg, fold == k)
  
  # Formatting for stan
  dat.stan <- c(list(
    
    # Number of observations, articles, & levels 
    N_test = dim(dat.test)[1],
    N_train = dim(dat.train)[1],
    N_article = max(dat.train$article_idx),
    L_sp = 3,
    L_tg = 4,
    L_dpi = 3,
    L_age = 3,
    L_cell = 3,
    
    # Toggle predictors on/off
    t_inc = 0,
    sg_inc = 1,
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
    sg_train = dat.train$val_sg_combined,
    sg_test = dat.test$val_sg_combined,
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
  log_lik[, dat.sg$fold == k] <- as.data.frame(fit.current$draws(log_lik_pars, 
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
TN <- sum(pred.df$true_val == 0 & pred.df$pred_val == 0) + TN_SG_neg
FN <- sum(pred.df$true_val == 1 & pred.df$pred_val == 0) + FN_SG_neg
F_score <- (2 * TP) / (2 * TP + FP + FN)
MCC <- ((TP * TN) - (FP * FN))/(sqrt(TP + FP)*sqrt(TP + FN)*sqrt(TN + FP)*sqrt(TN + FN))

# Percent of samples correctly classified, stratified by Pos/Neg/All 
percent_pos <- TP / (TP + FN)  * 100
percent_neg <- TN / (TN + FP) * 100
percent_total <- (TN + TP) / (TP + FP + TN + FN) * 100

# Store log_lik data frame for future reference
loglik_name <- paste0("log_lik_", "c1.2") # model name for future ref
assign(loglik_name, log_lik)


#### PSIS-LOO method ------------------------------------------------------------

dat.train <- dat.sg
dat.test <- dat.sg

# Formatting for stan
dat.stan <- c(list(
  
  # Number of observations, articles, & levels 
  N_test = dim(dat.test)[1],
  N_train = dim(dat.train)[1],
  N_article = max(dat.train$article_idx),
  L_sp = 3,
  L_tg = 4,
  L_dpi = 3,
  L_age = 3,
  L_cell = 3,
  
  # Toggle predictors on/off
  t_inc = 0,
  sg_inc = 1,
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
  sg_train = dat.train$val_sg_combined,
  sg_test = dat.test$val_sg_combined,
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

fit.sg.quant <- fit.current

# Store log likelihood
log_lik_pars <- paste0("log_lik[", 1:nrow(dat.sg), "]")
log_lik <- as.data.frame(fit.current$draws(log_lik_pars,format = "matrix"))

# Use loo package to extract elpd, ic, p_loo, MC_se, and pareto_k statistics
loo_current <- loo(as.matrix(log_lik), r_eff = NA)
elpd_loo_current <- loo_current$estimates[1, 1]
elpd_loo_se_current <- loo_current$estimates[1, 2]
ploo_current <- loo_current$estimates[2, 1]
ploo_se_current <- loo_current$estimates[2, 2]
pareto_k_max <- max(pareto_k_values(loo_current))


# Append statistics to performance table
next_row <- c(list(model = "q1.2",
                   predictors = "SG*",
                   elpd = elpd_current,
                   elpd_se = elpd_se_current,
                   elpd_loo = elpd_loo_current,
                   elpd_loo_se = elpd_loo_se_current,
                   p_loo_loo = ploo_current,
                   p_loo_loo_se = ploo_se_current,
                   pareto_k_max = pareto_k_max,
                   F_score = F_score,
                   MCC = MCC,
                   n_correct_test = TP + TN,
                   n_correct_train = n_corr_train,
                   n_total_test = TP + TN + FP + FN,
                   n_total_train = nrow(dat.sg) * (n_folds - 1),
                   percent_pos = percent_pos,
                   percent_neg = percent_neg,
                   percent_total = percent_total))
tblSX <- rbind(tblSX, next_row)

# Convert log likelihoods to loo object for future reference
loo_name <- paste0("loo_", "c1.2") # model name for future ref
assign(loo_name, loo_current)


## total RNA + sgRNA --------------------------------------------------------------

dat.both <- subset(dat.pred, pos_sg_combined == 1 & cens_total == 0)
log_lik <- as.data.frame(matrix(nrow = n_chains * n_iter / 2, ncol = nrow(dat.both)))

# Prep df to store predictions
pred.df <- list(true_val = c(),
                pred_val = c())


#### K-fold cross validation ---------------------------------------------------

n_corr_test <- 0
n_corr_train <- 0

for (k in 1:n_folds) {
  
  cat("Running fold ", k, "for the quantitative base model: T + SG \n \n")
  
  # For ease of formatting data for stan
  dat.train <- subset(dat.both, fold != k)
  dat.test <- subset(dat.both, fold == k)
  
  # Formatting for stan
  dat.stan <- c(list(
    
    # Number of observations, articles, & levels 
    N_test = dim(dat.test)[1],
    N_train = dim(dat.train)[1],
    N_article = max(dat.train$article_idx),
    L_sp = 3,
    L_tg = 4,
    L_dpi = 3,
    L_age = 3,
    L_cell = 3,
    
    # Toggle predictors on/off
    t_inc = 1,
    sg_inc = 1,
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
    sg_train = dat.train$val_sg_combined,
    sg_test = dat.test$val_sg_combined,
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
  log_lik[, dat.both$fold == k] <- as.data.frame(fit.current$draws(log_lik_pars, 
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
TN <- sum(pred.df$true_val == 0 & pred.df$pred_val == 0) + TN_Both_neg
FN <- sum(pred.df$true_val == 1 & pred.df$pred_val == 0) + FN_Both_neg
F_score <- (2 * TP) / (2 * TP + FP + FN)
MCC <- ((TP * TN) - (FP * FN))/(sqrt(TP + FP)*sqrt(TP + FN)*sqrt(TN + FP)*sqrt(TN + FN))

# Percent of samples correctly classified, stratified by Pos/Neg/All 
percent_pos <- TP / (TP + FN)  * 100
percent_neg <- TN / (TN + FP) * 100
percent_total <- (TN + TP) / (TP + FP + TN + FN) * 100

# Store log_lik data frame for future reference
loglik_name <- paste0("log_lik_", "c1.3") # model name for future ref
assign(loglik_name, log_lik)


#### PSIS-LOO method ------------------------------------------------------------

# Will run slower, but prevents needing an additional model file
dat.train <- dat.both
dat.test <- dat.both

# Formatting for stan
dat.stan <- c(list(
  
  # Number of observations, articles, & levels 
  N_test = dim(dat.test)[1],
  N_train = dim(dat.train)[1],
  N_article = max(dat.train$article_idx),
  L_sp = 3,
  L_tg = 4,
  L_dpi = 3,
  L_age = 3,
  L_cell = 3,
  
  # Toggle predictors on/off
  t_inc = 1,
  sg_inc = 1,
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
  sg_train = dat.train$val_sg_combined,
  sg_test = dat.test$val_sg_combined,
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

fit.both.quant <- fit.current

# Store log likelihood
log_lik_pars <- paste0("log_lik[", 1:nrow(dat.both), "]")
log_lik <- as.data.frame(fit.current$draws(log_lik_pars,format = "matrix"))

# Use loo package to extract elpd, ic, p_loo, MC_se, and pareto_k statistics
loo_current <- loo(as.matrix(log_lik), r_eff = NA)
elpd_loo_current <- loo_current$estimates[1, 1]
elpd_loo_se_current <- loo_current$estimates[1, 2]
ploo_current <- loo_current$estimates[2, 1]
ploo_se_current <- loo_current$estimates[2, 2]
pareto_k_max <- max(pareto_k_values(loo_current))

# Append statistics to performance table
next_row <- c(list(model = "q1.3",
                   predictors = "T + SG*",
                   elpd = elpd_current,
                   elpd_se = elpd_se_current,
                   elpd_loo = elpd_loo_current,
                   elpd_loo_se = elpd_loo_se_current,
                   p_loo_loo = ploo_current,
                   p_loo_loo_se = ploo_se_current,
                   pareto_k_max = pareto_k_max,
                   F_score = F_score,
                   MCC = MCC,
                   n_correct_test = TP + TN,
                   n_correct_train = n_corr_train,
                   n_total_test = TP + TN + FP + FN,
                   n_total_train = nrow(dat.both) * (n_folds - 1),
                   percent_pos = percent_pos,
                   percent_neg = percent_neg,
                   percent_total = percent_total))
tblSX <- rbind(tblSX, next_row)


# Convert log likelihoods to loo object for future reference
loo_name <- paste0("loo_", "c1.3") # model name for future ref
assign(loo_name, loo_current)



# Update the table -----------------------------------------------------------

## Round for visualization 
tblSX$F_score <- round(tblSX$F_score, digits = 2)
tblSX$MCC <- round(tblSX$MCC, digits = 2)
tblSX$percent_total <- round(tblSX$percent_total, digits = 1)
tblSX$percent_pos <- round(tblSX$percent_pos, digits = 1)
tblSX$percent_neg <- round(tblSX$percent_neg, digits = 1)


## Change column names
colnames(tblSX)[colnames(tblSX) == "percent_total"] <- "Overall Accuracy (%)"
colnames(tblSX)[colnames(tblSX) == "model"] <- "Model"
colnames(tblSX)[colnames(tblSX) == "predictors"] <- "Predictors"
colnames(tblSX)[colnames(tblSX) == "n_correct_test"] <- "Num. correct (test)"
colnames(tblSX)[colnames(tblSX) == "n_total_test"] <- "Num. total"
colnames(tblSX)[colnames(tblSX) == "percent_pos"] <- "Positive Accuracy (%)"
colnames(tblSX)[colnames(tblSX) == "percent_neg"] <- "Negative Accuracy (%)"

## Get table in correct format
tblSX <- tblSX %>% dplyr::select(c("Model", "Predictors", 
                                   "Overall Accuracy (%)",
                                   "Positive Accuracy (%)", 
                                   "Negative Accuracy (%)",
                                   "MCC"))


# Save -------------------------------------------------------------------------

write.csv(tblSX, file = "./outputs/tables/tblS9-culture-base-models.csv")


