# Generates and saves the fit for the best culture model, with either
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

dat <- read.csv("./data/pred-sg-data.csv")

dat.full <- subset(dat, !is.na(pos_inf) & cens_total %in% c(0, 1)) # for calculations below
dat.culture <- subset(dat, cens_total == 0 & !is.na(pos_inf)) # for fitting

# Convert columns to numerics
dat.culture$val_sg_pred <- as.numeric(dat.culture$val_sg_pred)
dat.culture$val_total <- as.numeric(dat.culture$val_total)

# Set articles as numeric factors
#dat.culture$article_idx <-  as.numeric(factor(dat.culture$article))
dat.culture$article_group <- dat.culture$article
dat.culture$article_group[dat.culture$article %in% c("Gabitzsch et al. 2021")] <- "Battelle Biomedical Research Center"
dat.culture$article_group[dat.culture$article %in% c("Aid et al. 2020",
                                     "Baum et al. 2020",
                                     "Chandrashekar et al. 2020",
                                     "Corbett et al. 2020",
                                     "Dagotto et al. 2021",
                                     "Francica et al. 2021",
                                     "Guebre Xabier et al. 2020",
                                     "He et al. 2021",
                                     "Jones et al. 2021",
                                     "Li Dapeng et al. 2021",
                                     "McMahan et al. 2020",
                                     "Mercado et al. 2020",
                                     "Patel et al. 2021",
                                     "Saunders et al. 2021",
                                     "Yu (Jingyou) et al. 2020",
                                     "Zost et al. 2020")] <- "Bioqual"
dat.culture$article_group[dat.culture$article %in% c("Lakshmanappa et al. 2021")] <- "California Primate Center"
dat.culture$article_group[dat.culture$article %in% c("Deng et al. 2020A",
                                     "Deng et al. 2020B",
                                     "Gao et al. 2020",
                                     "Yu (Pin) et al. 2020",
                                     "Guo et al. 2021",
                                     "Wang et al. 2020")] <- "CAMS (Chuan Qin)"
dat.culture$article_group[dat.culture$article %in% c("Rockx et al. 2020")] <- "Erasmus"
dat.culture$article_group[dat.culture$article %in% c("Yadav et al. 2021A",
                                     "Yadav et al. 2021B")] <- "ICMR-NIV, Pune"
dat.culture$article_group[dat.culture$article %in% c("Brouwer et al. 2021",
                                     "Maisonasse et al. 2021",
                                     "Maisonnasse et al. 2020",
                                     "Sokol et al. 2021")] <- "IDMIT"
dat.culture$article_group[dat.culture$article %in% c("Finch et al. 2020")] <- "IRF, Frederick"
dat.culture$article_group[dat.culture$article %in% c("Nagata et al. 2021")] <- "Japan NIID"
dat.culture$article_group[dat.culture$article %in% c("Kim et al. 2021",
                                     "Koo et al. 2020",
                                     "Seo et al. 2021")] <- "Korea Primate Center"
dat.culture$article_group[dat.culture$article %in% c("Liang et al. 2021",
                                     "Song et al. 2020",
                                     "Song et al. 2021",
                                     "An et al. 2021",
                                     "Jiao et al. 2021",
                                     "Nagata et al. 2021",
                                     "Lu et al. 2020A",
                                     "Lu et al. 2020B",
                                     "Yang et al. 2020")] <- "Kunming Primate Center / Kunming Institute of Zoology"
dat.culture$article_group[dat.culture$article %in% c("Bewley et al. 2020",
                                     "Lambe et al. 2021",
                                     "Rauch et al. 2020",
                                     "Salguero et al. 2021")] <- "Public Health England"
dat.culture$article_group[dat.culture$article %in% c("Boszormenyi et al. 2020",
                                     "Philippens et al. 2021",
                                     "Sanchez-Felipe et al. 2021")] <- "Rijswijk Primate Center"
dat.culture$article_group[dat.culture$article %in% c("Furuyama et al. 2022",
                                     "Hassan et al. 2021",
                                     "Munster et al. 2020",
                                     "Rosenke et al. 2020",
                                     "Speranza et al. 2020",
                                     "van Doremalen et al. 2020",
                                     "Williamson et al. 2020")] <- "RML"
dat.culture$article_group[dat.culture$article %in% c("Furuyama et al. 2022",
                                     "Hassan et al. 2021",
                                     "Munster et al. 2020",
                                     "Rosenke et al. 2020",
                                     "Speranza et al. 2020",
                                     "van Doremalen et al. 2020",
                                     "Williamson et al. 2020")] <- "RML"
dat.culture$article_group[dat.culture$article %in% c("Gorman et al. 2021",
                                     "Singh et al. 2020",
                                     "Vogel et al. 2021")] <- "Southwest Primate Center (San Antonio)"
dat.culture$article_group[dat.culture$article %in% c("Arunachalam et al. 2021",
                                     "Blair et al. 2021",
                                     "Fahlberg et al. 2020",
                                     "Huang et al. 2021")] <- "Tulane Primate Center"
dat.culture$article_group[dat.culture$article %in% c("Cross et al. 2020",
                                     "Cross et al. 2021",
                                     "Woolsey et al. 2020")] <- "UTMB, Galveston"
dat.culture$article_group[dat.culture$article %in% c("Hoang et al. 2021",
                                     "Routhu et al. 2021")] <- "Yerkes Primate Center"

dat.culture$article_idx <-  as.numeric(factor(dat.culture$article_group))


# Set TG indicator based on total RNA protocol: ordered by decreasing expected RNA quantity
dat.culture$tg_idx[dat.culture$tg_total %in% c("N")] <- 1
dat.culture$tg_idx[dat.culture$tg_total %in% c("E")] <- 2
dat.culture$tg_idx[dat.culture$tg_total %in% c("S")] <- 3

# Set a seed
options(mc.cores = parallel::detectCores())
seed <- 999

# Set chains & iteractions
n_iter <- 4000
n_chains <- 6


# Non-Informative priors -------------------------------------------------------

## Load and run models ---------------------------------------------------------

# Load & compile models
model.best <- cmdstan_model('./code/stan-model-culture-best.stan')
model.simple <- cmdstan_model('./code/stan-model-culture-simplest.stan')
model.article <- cmdstan_model('./code/stan-model-culture-best-with-article.stan')


# Formatting for stan
dat.stan <- c(list(
  
  # Number of observations, articles, & levels 
  N = dim(dat.culture)[1],
  N_article = max(dat.culture$article_idx),
  L_sp = max(dat.culture$sp_idx),
  L_tg = max(dat.culture$tg_idx),
  L_dpi = max(dat.culture$dpi_idx),
  L_age = max(dat.culture$age_idx),
  L_cell = max(dat.culture$cell_idx),
  
  # Response variables
  pos_inf = dat.culture$pos_inf,
  
  # Predictor variables
  t = dat.culture$val_total,
  sg = dat.culture$val_sg_pred,
  dose = dat.culture$log10_dose_pfu,
  tg = dat.culture$tg_idx,
  dpi = dat.culture$dpi_idx,
  age = dat.culture$age_idx,
  sp = dat.culture$sp_idx,
  cell = dat.culture$cell_idx,
  assay = dat.culture$assay_idx,
  article = dat.culture$article_idx))

# Run simplest model
fit.simple.noninf <- model.simple$sample(data = dat.stan,
                                         chains = n_chains,
                                         parallel_chains = n_chains,
                                         iter_warmup = n_iter / 2,
                                         iter_sampling = n_iter / 2,
                                         refresh = 1000,
                                         seed = seed)

# Run the best model
fit.best.noninf <- model.best$sample(data = dat.stan,
                                     chains = n_chains,
                                     parallel_chains = n_chains,
                                     iter_warmup = n_iter / 2,
                                     iter_sampling = n_iter / 2,
                                     refresh = 1000,
                                     seed = seed)


## Save fits -------------------------------------------------------------------

fit.simple.noninf$save_object(file = "./outputs/fits/fit-culture-simple-noninf-priors.RDS")
fit.best.noninf$save_object(file = "./outputs/fits/fit-culture-best-noninf-priors.RDS")



# Informative priors -----------------------------------------------------------

# Run models ----------------------------------------------------------

# Must manually change priors in model file before running!!
model.best <- cmdstan_model('./code/stan-model-culture-best.stan')
model.simple <- cmdstan_model('./code/stan-model-culture-simplest.stan')


# Run the best model
fit.best.inf <- model.best$sample(data = dat.stan,
                                  chains = n_chains,
                                  parallel_chains = n_chains,
                                  iter_warmup = n_iter / 2,
                                  iter_sampling = n_iter / 2,
                                  refresh = 1000,
                                  seed = seed)

# Run the simplest model
fit.simple.inf <- model.simple$sample(data = dat.stan,
                                      chains = n_chains,
                                      parallel_chains = n_chains,
                                      iter_warmup = n_iter / 2,
                                      iter_sampling = n_iter / 2,
                                      refresh = 1000,
                                      seed = seed)

fit.article.inf <- model.article$sample(data = dat.stan,
                                        chains = n_chains,
                                        parallel_chains = n_chains,
                                        iter_warmup = n_iter / 2,
                                        iter_sampling = n_iter / 2,
                                        refresh = 1000,
                                        seed = seed)


# Save fits --------------------------------------------------------------------

fit.best.inf$save_object(file = "./outputs/fits/fit-culture-best-inf-priors.RDS")
fit.simple.inf$save_object(file = "./outputs/fits/fit-culture-simple-inf-priors.RDS")
fit.article.inf$save_object(file = "./outputs/fits/fit-culture-article-inf-priors.RDS")

