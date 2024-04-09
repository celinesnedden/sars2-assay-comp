# Generates and saves the fit for the best sgRNA hurdle model, with either
#     informative or non-informative priors
# Note that priors must be manually changed in the model file, and the model must 
#     then be recompiled


# Prep environment -------------------------------------------------------------

# Install & load necessary packages
req_pkgs <- c("cmdstanr", "tidyverse")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependences = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))

`%notin%` <- Negate(`%in%`)  # for convenience


# Load & prep data -------------------------------------------------------------

dat <- read.csv("./data/clean-data.csv")

# Removes NA & positively censored total RNA data for this analysis
dat <- subset(dat, pos_sg %in% c(0, 1) & cens_total == 0) 

# Set sgRNA positive but censored values to -9 (to screen out from linear component)
dat$val_sg[dat$cens_sg == 3] <- -9

# Convert columns to numerics
dat$val_sg <- as.numeric(dat$val_sg)
dat$val_total <- as.numeric(dat$val_total)

# Convert articles into factor
dat$article_idx <- as.numeric(factor(dat$article))

# Set a seed
seed <- 126402


# Informative priors -----------------------------------------------------------

## Load and run models ----------------------------------------------------------

# Load & compile models
model.simple <- cmdstan_model('./code/stan-model-sgRNA-simplest.stan')
model.best <- cmdstan_model('./code/stan-model-sgRNA-full-best.stan')

# Set chains and iterations
n_chains <- 6
n_iter <- 4000

# Prep data for Stan
dat.stan <- c(list(
  
  # Number of observations, articles, & levels 
  N = nrow(dat),
  N_lin = nrow(subset(dat, cens_sg == 0)),
  N_article = max(dat$article_idx),
  L_sp = max(dat$sp_idx),
  L_tg = max(dat$tg_idx),
  L_dpi = max(dat$dpi_idx),
  
  # Response variables
  pos_sg = dat$pos_sg,
  val_sg = dat$val_sg,
  
  # Predictor variables
  t_val = dat$val_total,
  dose = dat$log10_dose_pfu,
  tg = dat$tg_idx,
  dpi = dat$dpi_idx,
  age = dat$age_idx,
  st = dat$st_idx,
  sp = dat$sp_idx,
  
  # Article hierarchy
  article = dat$article_idx,
  
  # Predictions, dummy & not used here
  N_new = 1,
  t_val_new = 1,
  dose_new = 1,
  dpi_new = 1,
  st_new = 1,
  sp_new = 1,
  tg_new = 1,
  age_new = 1 ))

# Run simplest model
fit.simple <- model.simple$sample(data = dat.stan,
                                  chains = n_chains,
                                  parallel_chains = n_chains,
                                  iter_warmup = n_iter / 2,
                                  iter_sampling = n_iter / 2,
                                  refresh = 1000,
                                  seed = seed)

# Run the best model
fit.best <- model.best$sample(data = dat.stan,
                              chains = n_chains,
                              parallel_chains = n_chains,
                              iter_warmup = n_iter / 2,
                              iter_sampling = n_iter / 2,
                              refresh = 1000,
                              seed = seed)


## Save fits --------------------------------------------------------------------

fit.simple$save_object(file = "./outputs/fits/fit-sgRNA-simple-inf-priors.RDS")
fit.best$save_object(file = "./outputs/fits/fit-sgRNA-best-inf-priors.RDS")



# Non-informative priors -------------------------------------------------------

## Load and run models ----------------------------------------------------------

# Must manually change priors in model file before running!!
model.simple <- cmdstan_model('./code/stan-model-sgRNA-simplest.stan')
model.best <- cmdstan_model('./code/stan-model-sgRNA-full-best.stan')

# Set chains and iterations
n_chains <- 6
n_iter <- 4000

# Prep data for Stan
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
  
  # Response variables
  pos_sg = dat$pos_sg,
  val_sg = dat$val_sg,
  
  # Predictor variables
  t_val = dat$val_total,
  dose = dat$log10_dose_pfu,
  tg = dat$tg_idx,
  dpi = dat$dpi_idx,
  age = dat$age_idx,
  st = dat$st_idx,
  sp = dat$sp_idx,
  
  # Article hierarchy
  article = dat$article_idx,
  
  # Predictions, dummy & not used here
  N_new = 1,
  t_val_new = 1,
  dose_new = 1,
  dpi_new = 1,
  st_new = 1,
  sp_new = 1,
  tg_new = 1,
  age_new = 1 ))


# Run the best model
fit.simple.noninf <- model.simple$sample(data = dat.stan,
                                         chains = n_chains,
                                         parallel_chains = n_chains,
                                         iter_warmup = n_iter / 2,
                                         iter_sampling = n_iter / 2,
                                         refresh = 1000,
                                         seed = seed)

fit.best.noninf <- model.best$sample(data = dat.stan,
                                     chains = n_chains,
                                     parallel_chains = n_chains,
                                     iter_warmup = n_iter / 2,
                                     iter_sampling = n_iter / 2,
                                     refresh = 1000,
                                     seed = seed)


## Save fits --------------------------------------------------------------------

fit.simple.noninf$save_object(file = "./outputs/fits/fit-sgRNA-simple-noninf-priors.RDS")
fit.best.noninf$save_object(file = "./outputs/fits/fit-sgRNA-best-noninf-priors.RDS")

