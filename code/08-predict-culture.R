# Generates predictions for all culture data using the best model

# Load data --------------------------------------------------------------------

dat <- read.csv("./data/pred-sg-data.csv")

# Prep data --------------------------------------------------------------------

dat.full <- subset(dat, !is.na(pos_inf) & cens_total %in% c(0, 1))

# Set TG indicator based on total RNA protocol: ordered by decreasing RNA quantity
dat.full$tg_idx[dat.full$tg_total %in% c("N")] <- 1
dat.full$tg_idx[dat.full$tg_total %in% c("E")] <- 2
dat.full$tg_idx[dat.full$tg_total %in% c("S")] <- 3

# Subset to tRNA positive and negative samples
dat.pos <- subset(dat.full, cens_total == 0 & !is.na(pos_inf))
dat.neg <- subset(dat.full, cens_total == 1 & !is.na(pos_inf))

# Other data not included
dat.rest <- subset(dat, is.na(pos_inf) | cens_total %notin% c(0, 1))
nrow(dat.rest) + nrow(dat.full) == nrow(dat) # SHOULD BE TRUE. If false, problem!


# Load model fits --------------------------------------------------------------

fit.best <- readRDS(file = "./outputs/fits/fit-culture-best-inf-priors.RDS")
fit.simple <- readRDS(file = "./outputs/fits/fit-culture-simple-inf-priors.RDS")


# Extract culture predictions ---------------------------------------------------

# Determine which are predicted as < and > LOD
prob_pos_pred <- paste0("p_pos[", 1:nrow(dat.pos), "]") 
prob_pos_df <- fit.best$draws(prob_pos_pred, format = "df")
prob_pos_mean <- as.data.frame(colMeans(prob_pos_df))
prob_pos_mean <- prob_pos_mean[rownames(prob_pos_mean) %notin% c(".chain",
                                                                 ".iteration",
                                                                 ".draw"), ]

threshold <- 0.5 
pos_pred <- prob_pos_mean
pos_pred[prob_pos_mean >= threshold] <- 1
pos_pred[prob_pos_mean < threshold] <- 0


# Append predictions to data ---------------------------------------------------

# tRNA positive predictions from model fit
dat.pos$pos_inf_pred <- pos_pred
dat.pos$pos_inf_prob_pred <- prob_pos_mean

# tRNA negative predictions automatically classified as culture negative
dat.neg$pos_inf_pred <- 0
dat.neg$pos_inf_prob_pred <- 0

# Add empty columns for dat.rest
dat.rest$pos_inf_pred <- NA
dat.rest$pos_inf_prob_pred <- NA


# Combine them back into one df ------------------------------------------------

dat.combined <- rbind(dat.pos, dat.neg, dat.rest)

# Check for any problems
nrow(dat.combined) == nrow(dat) # TRUE
unique(dat.combined$X %in% dat$X) # TRUE: all dat.combined rows are in dat
unique(dat$X %in% dat.combined$X) # TRUE: all dat rows are in dat.combined
unique(sort(dat$X) == sort(dat.combined$X))


# Reorder columns & rows for visibility ----------------------------------------

dat.combined <- dat.combined %>%
  select(article, indiv, dpi, dpi_idx,
         sample_type, st_idx, sample_rep,
         val_total, pos_total, cens_total, 
         val_sg, pos_sg, cens_sg, 
         val_sg_pred, pos_sg_pred, pos_sg_prob_pred,
         val_sg_combined, pos_sg_combined,
         val_inf, pos_inf, pos_inf_pred, pos_inf_prob_pred,
         lodq_sg, lodq_total,
         tg_total, tg_sg, tg_idx, 
         cell_line, cell_idx, 
         culture_assay, assay_idx,
         animal_species, sp_idx,
         age_class_stand, age_years_rep, age_class_rep, age_idx,
         sex, sex_idx, 
         units_total, units_sg, units_inf,
         lod_total, loq_total, lod_sg, loq_sg,
         everything()) %>%
  arrange(X)


# Export -----------------------------------------------------------------------

write.csv(dat.combined, file = "./data/pred-culture-data.csv",
          row.names = FALSE)



# Get predictions also for SIMPLE MODEL ----------------------------------------

# Load data --------------------------------------------------------------------

dat <- read.csv("./data/pred-sg-data.csv")

# Prep data --------------------------------------------------------------------

dat.full <- subset(dat, !is.na(pos_inf) & cens_total %in% c(0, 1))

# Set TG indicator based on total RNA protocol: ordered by decreasing RNA quantity
dat.full$tg_idx[dat.full$tg_total %in% c("N")] <- 1
dat.full$tg_idx[dat.full$tg_total %in% c("E")] <- 2
dat.full$tg_idx[dat.full$tg_total %in% c("S")] <- 3

# Subset to tRNA positive and negative samples
dat.pos <- subset(dat.full, cens_total == 0 & !is.na(pos_inf))
dat.neg <- subset(dat.full, cens_total == 1 & !is.na(pos_inf))

# Other data not included
dat.rest <- subset(dat, is.na(pos_inf) | cens_total %notin% c(0, 1))
nrow(dat.rest) + nrow(dat.full) == nrow(dat) # SHOULD BE TRUE. If false, problem!


# Load model fits --------------------------------------------------------------

fit.best <- readRDS(file = "./outputs/fits/fit-culture-best-inf-priors.RDS")
fit.simple <- readRDS(file = "./outputs/fits/fit-culture-simple-inf-priors.RDS")


# Extract culture predictions ---------------------------------------------------

# Determine which are predicted as < and > LOD
prob_pos_pred <- paste0("p_pos[", 1:nrow(dat.pos), "]") 
prob_pos_df <- fit.simple$draws(prob_pos_pred, format = "df")
prob_pos_mean <- as.data.frame(colMeans(prob_pos_df))
prob_pos_mean <- prob_pos_mean[rownames(prob_pos_mean) %notin% c(".chain",
                                                                 ".iteration",
                                                                 ".draw"), ]

threshold <- 0.5 
pos_pred <- prob_pos_mean
pos_pred[prob_pos_mean >= threshold] <- 1
pos_pred[prob_pos_mean < threshold] <- 0


# Append predictions to data ---------------------------------------------------

# tRNA positive predictions from model fit
dat.pos$pos_inf_pred <- pos_pred
dat.pos$pos_inf_prob_pred <- prob_pos_mean

# tRNA negative predictions automatically classified as culture negative
dat.neg$pos_inf_pred <- 0
dat.neg$pos_inf_prob_pred <- 0

# Add empty columns for dat.rest
dat.rest$pos_inf_pred <- NA
dat.rest$pos_inf_prob_pred <- NA


# Combine them back into one df ------------------------------------------------

dat.combined <- rbind(dat.pos, dat.neg, dat.rest)

# Check for any problems
nrow(dat.combined) == nrow(dat) # TRUE
unique(dat.combined$X %in% dat$X) # TRUE: all dat.combined rows are in dat
unique(dat$X %in% dat.combined$X) # TRUE: all dat rows are in dat.combined
unique(sort(dat$X) == sort(dat.combined$X))


# Reorder columns & rows for visibility ----------------------------------------

dat.combined <- dat.combined %>%
  select(article, indiv, dpi, dpi_idx,
         sample_type, st_idx, sample_rep,
         val_total, pos_total, cens_total, 
         val_sg, pos_sg, cens_sg, 
         val_sg_pred, pos_sg_pred, pos_sg_prob_pred,
         val_sg_combined, pos_sg_combined,
         val_inf, pos_inf, pos_inf_pred, pos_inf_prob_pred,
         lodq_sg, lodq_total,
         tg_total, tg_sg, tg_idx, 
         cell_line, cell_idx, 
         culture_assay, assay_idx,
         animal_species, sp_idx,
         age_class_stand, age_years_rep, age_class_rep, age_idx,
         sex, sex_idx, 
         units_total, units_sg, units_inf,
         lod_total, loq_total, lod_sg, loq_sg,
         everything()) %>%
  arrange(X)


# Export -----------------------------------------------------------------------

write.csv(dat.combined, file = "./data/pred-culture-data-simple.csv",
          row.names = FALSE)





