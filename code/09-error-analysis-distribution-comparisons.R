# This compares the distributions of values for various error types


# Prep environment -------------------------------------------------------------

## Install CmdStanR, CmdStan, and RStan for model compilation ------------------

#  We recommend you review the installation instructions at:
#     https://mc-stan.org/cmdstanr/articles/cmdstanr.html 
#     https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started 


library("cmdstanr")


## Install other packages ------------------------------------------------------

req_pkgs <- c("loo", "stringr", "tidyverse", "ggridges", "bayesplot",
              "gtools", "wesanderson", "cowplot", "ggpubr", "matrixStats")
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

seed <- 29


# Load model -------------------------------------------------------------------

stan.model <- cmdstan_model('./code/stan-model-distribution-comparison.stan')

# Prep model -------------------------------------------------------------------

dat <- read.csv("./data/pred-culture-data.csv")
dat <- subset(dat, !is.na(pos_inf) & !is.na(pos_total))

dat$pos_inf_pred_class <- NA
dat$pos_inf_pred_class[dat$pos_inf == 1 & dat$pos_inf_pred == 1] <- "True Positive"
dat$pos_inf_pred_class[dat$pos_inf == 0 & dat$pos_inf_pred == 0] <- "True Negative"
dat$pos_inf_pred_class[dat$pos_inf == 1 & dat$pos_inf_pred == 0] <- "False Negative"
dat$pos_inf_pred_class[dat$pos_inf == 0 & dat$pos_inf_pred == 1] <- "False Positive"

dat.test <- subset(dat, cens_total == 0)
dat.test$val_total <- as.numeric(dat.test$val_total)

dat.test$error_type_idx[dat.test$pos_inf_pred_class == "True Positive"] <- 1
dat.test$error_type_idx[dat.test$pos_inf_pred_class == "True Negative"] <- 2
dat.test$error_type_idx[dat.test$pos_inf_pred_class == "False Negative"] <- 3
dat.test$error_type_idx[dat.test$pos_inf_pred_class == "False Positive"] <- 4


# Run model -------------------------------------------------------------------

# Set chains and iterations
n_chains <- 6
n_iter <- 4000

dat.stan <- c(list(
  
  # Number of observations, articles, & levels 
  N = nrow(dat.test),
  N_types = max(dat.test$error_type_idx),
  y = dat.test$val_total, 
  y_type = dat.test$error_type_idx))

fit.inf <- stan.model$sample(data = dat.stan,
                             chains = n_chains,
                             parallel_chains = n_chains,
                             iter_warmup = n_iter / 2,
                             iter_sampling = n_iter / 2,
                             refresh = 1000,
                             seed = 78)
fit.inf


# Compare estimates of distribution means --------------------------------------

pars.inf <- c("mu[1]", "mu[2]", "mu[3]", "mu[4]")
pars.inf.df <- fit.inf$draws(variables = pars.inf,
                             inc_warmup = FALSE,
                             format = "df")
pars.inf.df <- as.data.frame(pars.inf.df)

# Histograms comparing distributions
hist(pars.inf.df$`mu[1]` - pars.inf.df$`mu[3]`, main = "TP vs. FN")
hist(pars.inf.df$`mu[2]` - pars.inf.df$`mu[3]`, main = "TN vs. FN")

hist(pars.inf.df$`mu[1]` - pars.inf.df$`mu[4]`, main = "TP vs. FP")
hist(pars.inf.df$`mu[2]` - pars.inf.df$`mu[4]`, main = "TN vs. FP")

# TP vs. FN
median(pars.inf.df$`mu[3]` - pars.inf.df$`mu[1]`); quantile(pars.inf.df$`mu[3]` - pars.inf.df$`mu[1]`, c(0.05, 0.95))
# TN vs. FN
median(pars.inf.df$`mu[3]` - pars.inf.df$`mu[2]`); quantile(pars.inf.df$`mu[3]` - pars.inf.df$`mu[2]`, c(0.05, 0.95))

# TN vs. FP
median(pars.inf.df$`mu[4]` - pars.inf.df$`mu[2]`); quantile(pars.inf.df$`mu[4]` - pars.inf.df$`mu[2]`, c(0.05, 0.95))

# TP vs. FP
median(pars.inf.df$`mu[4]` - pars.inf.df$`mu[1]`); quantile(pars.inf.df$`mu[4]` - pars.inf.df$`mu[1]`, c(0.05, 0.95))




