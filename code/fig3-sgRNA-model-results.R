# Generates figure 3 from the manuscript
#     Includes fit lines drawn from the best and simplest models for all 
#        selected predictors, for the logistic & full PCR models
#     Includes some performance comparisons (CDFs of probabilities,
#        median error distributions, & obs. vs. predicted values)


# Prep environment -------------------------------------------------------------

# Install & load color palette & ggplot 
req_pkgs <- c("wesanderson", "ggplot2", "cmdstanr", "tidyverse", "ggridges",
              "ggpubr", "patchwork", "RColorBrewer", "scales")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependences = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))


# Extract desired colored palette
zissou_pal <- wes_palette("Zissou1", 12, type = "continuous")
brewer_pal <- brewer.pal(10, "RdBu")
purple_pal <- rev(c(brewer.pal(10, "RdBu")[1:4], "white"))


# Function for convenience
`%notin%` <- Negate(`%in%`)  

# Set number of lines to include in figures
n_lines <- 300


# Load & prep data -------------------------------------------------------------

dat <- read.csv("./data/clean-data.csv")
dat <- subset(dat, pos_sg %in% c(0, 1) & cens_total %in% c(0, 1)) 

# Count TN/FN for total RNA negative data 
dat.neg <- subset(dat, cens_total == 1)
TN_neg <- nrow(subset(dat.neg, pos_sg == 0))
FN_neg <- nrow(subset(dat.neg, pos_sg == 1))

# Set sgRNA positive but censored values to -9 (to screen out from linear component)
dat$val_sg[dat$cens_sg == 3] <- -9

# Convert columns to numerics
dat$val_sg <- as.numeric(dat$val_sg)
dat$val_total <- as.numeric(dat$val_total)


## Add jitter to < LOD values for visualization --------------------------------

dat$val_sg_jitter <- NA
dat$pos_sg_jitter <- NA
for (row_num in 1:nrow(dat)){
  if (dat$pos_sg[row_num] == 1){
    # sgRNA < LOD
    dat$val_sg_jitter[row_num] <- runif(1, 1.02, 1.18)
    dat$pos_sg_jitter[row_num] <- runif(1, 0.52, 0.98)
  }
  else if (dat$pos_sg[row_num] == 0) {
    # sgRNA > LOD
    dat$val_sg_jitter[row_num] <- runif(1, -0.18, -0.02)
    dat$pos_sg_jitter[row_num] <- runif(1, 0.02, 0.48)
  }
}


# Load model fits --------------------------------------------------------------

fit.simple <- readRDS(file = "./outputs/fits/fit-sgRNA-simple-inf-priors.RDS")
fit.best <- readRDS(file = "./outputs/fits/fit-sgRNA-best-inf-priors.RDS")



# LOGISTIC COMPONENT PANELS ----------------------------------------------------

## 3A: Raw Data ----------------------------------------------------------------

### Plot ------------------------------------------------------------------------

dat$val_total[dat$cens_total == 1] <- runif(length(dat$val_total[dat$cens_total == 1]),
                                                   -1.9, -0.1)
dat <- subset(dat, val_total != -1.1800000) # Because axes dont accomodate this point


fig.data <- ggplot() + 
  geom_rect(aes(xmin = -2, xmax = 0, ymin = 0, ymax = 1), fill = "grey56", color = "black") +
  geom_hline(yintercept = 0.5, linewidth = 0.6, color = "black") +
  geom_vline(xintercept = 0, linewidth = 1, color = "black") +
  geom_point(data = dat,
             aes(x = as.numeric(val_total), 
                 y = as.numeric(pos_sg_jitter),
                 fill = as.character(pos_sg)),
             alpha = 0.6, size = 1.3, shape = 21) +
  scale_fill_manual(values = c("0" = brewer_pal[9], "1" = brewer_pal[2])) +
  labs(y = "",
       x = "log10 total RNA copies / sample",
       color = "log10 \n dose \n (pfu)",
       tag = "A") +
  guides(fill = "none",
         color = guide_colorbar(frame.colour = "black",
                                title.vjust = 0)) +
  theme(legend.position = "none",
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1, size = 7.5),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(colour = c("transparent", rep("black", 7))),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        plot.margin = margin(0, 4, 0, 0)) + 
  scale_y_continuous(breaks = c(0.25, 0.75), limits = c(0, 1),
                     labels = c("sgRNA\nNegative", "sgRNA\nPositive"), 
                     expand = c(0, 0)) +
  scale_x_continuous(breaks = c(-1, 0, 2, 4, 6, 8, 10), limits = c(-2, 10),
                     labels = c("< LOD", 0, 2, 4, 6, 8, 10),
                     expand = c(0, 0)) +
  geom_hline(yintercept = 1, size = 0.6, color = "black") +
  geom_hline(yintercept = 0, size = 0.6, color = "black"); fig.data


## 3B: Model probabilities: intermediate sgRNA ---------------------------------

# Set threshold value for classifying as > or < the LOD
threshold <- 0.5

# Subset to data total RNA positive data
dat.sg <- subset(dat, cens_total == 0)

### Prep simple model results --------------------------------------------------

# Extract mean model predictions of classification
prob_pos <- paste0("prob_pos[", 1:nrow(dat.sg), "]") 
prob_df <- fit.simple$draws(prob_pos, format = "df")
prob_mean <- as.data.frame(colMeans(prob_df))
prob_mean <- prob_mean[rownames(prob_mean) %notin% c(".chain", ".iteration", ".draw"), ]

# Compare predictions with true observed classification
prediction <- prob_mean
prediction[prob_mean >= threshold] <- 1
prediction[prob_mean < threshold] <- 0

# Probabilities of true positives
pos_pred <- prob_mean[dat.sg$pos_sg == 1 & dat.sg$val_total <= 6 &
                        dat.sg$val_total >= 4]

# Probabilities of true negatives
neg_pred <- prob_mean[dat.sg$pos_sg == 0 & dat.sg$val_total <= 6 &
                        dat.sg$val_total >= 4]


# Place in a data frame
prediction.df <- as.data.frame(matrix(data = NA, 
                                      nrow = length(c(pos_pred, neg_pred)),
                                      ncol = 3))
prediction.df[, 1] <- "Simple"
prediction.df[, 2] <- c(pos_pred, neg_pred)
prediction.df[, 3] <- c(rep("Positive", length(pos_pred)),
                        rep("Negative", length(neg_pred)))



### Prep best model results -----------------------------------------------------

# Extract mean model predictions of classification
prob_pos <- paste0("prob_pos[", 1:nrow(dat.sg), "]") 
prob_df <- fit.best$draws(prob_pos, format = "df")
prob_mean <- as.data.frame(colMeans(prob_df))
prob_mean <- prob_mean[rownames(prob_mean) %notin% c(".chain", ".iteration", ".draw"), ]

# Compare predictions with true observed classification
prediction <- prob_mean
prediction[prob_mean >= threshold] <- 1
prediction[prob_mean < threshold] <- 0

# Probabilities of true positives
pos_pred <- prob_mean[dat.sg$pos_sg == 1 & dat.sg$val_total <= 6 &
                        dat.sg$val_total >= 4]

# Probabilities of true negatives
neg_pred <- prob_mean[dat.sg$pos_sg == 0 & dat.sg$val_total <= 6 &
                        dat.sg$val_total >= 4]


# Place in joint data frame
prediction.df <- rbind(prediction.df, 
                       cbind(rep("Best", length(c(pos_pred, neg_pred))),
                             c(pos_pred, neg_pred),
                             c(rep("Positive", length(pos_pred)), 
                               rep("Negative", length(neg_pred)))))
colnames(prediction.df) <- c("model", "prob_mean", "posneg")


### Plot -----------------------------------------------------------------------

# Plot all observed sgRNA negatives
fig3B.1 <- ggplot(subset(prediction.df, posneg == "Negative")) +
  geom_vline(xintercept = threshold, size = 0.7, color = zissou_pal[1], alpha = 0.8,
             linetype = "dashed") +
  geom_density(aes(x = (1-as.numeric(prob_mean)),
                   fill = factor(model, levels = c("Simple", "Best")),
                   color = factor(model, levels = c("Simple", "Best"))),
               alpha = 0.4) +
  facet_grid(.~ "Observed Negative") +
  scale_fill_manual(values = c("Best" = "#85B068",
                               "Simple" = "#0F4C5E")) +
  scale_color_manual(values = c("Best" = "#56783F",
                                "Simple" = "#0F4C5E")) +
  guides(color = "none",
         fill = guide_legend(override.aes = list(alpha = 0.7,
                                                 color = c("#0F4C5E", "#56783F")))) +
  labs(x = "Chance of Negative (%)", y = "Density",
       fill = "Model", tag = "B") + 
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0),
                     breaks = seq(0, 1, 0.25), 
                     labels = c("0","25", "50", "75", "")) +
  scale_y_continuous(limits = c(0, 7), breaks = seq(0, 6, 1),
                     expand = c(0, 0),
                     position = "left") +
  theme(legend.position = c(0.22, 0.7), 
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        legend.background = element_rect(fill = 'white'),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 5, 0, 0)); fig3B.1

# Plot all observed sgRNA positives
fig3B.2 <- ggplot(subset(prediction.df, posneg == "Positive")) +
  geom_vline(xintercept = threshold, size = 0.7, color = zissou_pal[1], alpha = 0.8,
             linetype = "dashed") +
  geom_density(aes(x = (as.numeric(prob_mean)),
                   fill = factor(model, levels = c("Simple", "Best")),
                   color = factor(model, levels = c("Simple", "Best"))),
               alpha = 0.4) +
  facet_grid(.~ "Observed Positive") +
  scale_fill_manual(values = c("Best" = "#85B068",
                               "Simple" = "#0F4C5E")) +
  scale_color_manual(values = c("Best" = "#56783F",
                                "Simple" = "#0F4C5E")) +
  guides(color = "none",
         fill = guide_legend(override.aes = list(alpha = 0.6,
                                                 color = c( "#0F4C5E", "#56783F")))) +
  labs(x = "Chance of Positive (%)", y = "Density",
       fill = "Model") + 
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0),
                     breaks = seq(0, 1, 0.25), 
                     labels = c("0","25", "50", "75", "100")) +
  scale_y_continuous(limits = c(0, 7), breaks = seq(0, 8, 1),
                     expand = c(0, 0),
                     position = "left") +
  theme(legend.position = "none", 
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        legend.background = element_rect(fill = 'white'),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 0, 0, 0)); fig3B.2

fig.perf <- fig3B.1 + fig3B.2 + plot_layout(ncol = 2); fig.perf





## 3B: Model probabilities -----------------------------------------------------

# Set threshold value for classifying as > or < the LOD
threshold <- 0.5

# Subset to data total RNA positive data
dat.sg <- subset(dat, cens_total == 0)

### Prep simple model results --------------------------------------------------

# Extract mean model predictions of classification
prob_pos <- paste0("prob_pos[", 1:nrow(dat.sg), "]") 
prob_df <- fit.simple$draws(prob_pos, format = "df")
prob_mean <- as.data.frame(colMeans(prob_df))
prob_mean <- prob_mean[rownames(prob_mean) %notin% c(".chain", ".iteration", ".draw"), ]

# Compare predictions with true observed classification
prediction <- prob_mean
prediction[prob_mean >= threshold] <- 1
prediction[prob_mean < threshold] <- 0

# Probabilities of true positives
pos_pred <- prob_mean[dat.sg$pos_sg == 1]

# Probabilities of true negatives
neg_pred <- prob_mean[dat.sg$pos_sg == 0]

# Place in a data frame
prediction.df <- as.data.frame(matrix(data = NA, 
                                      nrow = length(c(pos_pred, neg_pred)),
                                      ncol = 3))
prediction.df[, 1] <- "Simple"
prediction.df[, 2] <- c(pos_pred, neg_pred)
prediction.df[, 3] <- c(rep("Positive", length(pos_pred)),
                        rep("Negative", length(neg_pred)))



### Prep best model results -----------------------------------------------------

# Extract mean model predictions of classification
prob_pos <- paste0("prob_pos[", 1:nrow(dat.sg), "]") 
prob_df <- fit.best$draws(prob_pos, format = "df")
prob_mean <- as.data.frame(colMeans(prob_df))
prob_mean <- prob_mean[rownames(prob_mean) %notin% c(".chain", ".iteration", ".draw"), ]

# Compare predictions with true observed classification
prediction <- prob_mean
prediction[prob_mean >= threshold] <- 1
prediction[prob_mean < threshold] <- 0

# Probabilities of true positives
pos_pred <- prob_mean[dat.sg$pos_sg == 1]

# Probabilities of true negatives
neg_pred <- prob_mean[dat.sg$pos_sg == 0]


# Place in joint data frame
prediction.df <- rbind(prediction.df, 
                       cbind(rep("Best", length(c(pos_pred, neg_pred))),
                             c(pos_pred, neg_pred),
                             c(rep("Positive", length(pos_pred)), 
                               rep("Negative", length(neg_pred)))))
colnames(prediction.df) <- c("model", "prob_mean", "posneg")


### Plot -----------------------------------------------------------------------

# Plot all observed sgRNA negatives
fig3B.1 <- ggplot(subset(prediction.df, posneg == "Negative")) +
  geom_vline(xintercept = threshold, size = 0.7, color = zissou_pal[1], alpha = 0.8,
             linetype = "dashed") +
  geom_density(aes(x = (1-as.numeric(prob_mean)),
                   fill = factor(model, levels = c("Simple", "Best")),
                   color = factor(model, levels = c("Simple", "Best"))),
               alpha = 0.4) +
  facet_grid(.~ "Observed Negative") +
  scale_fill_manual(values = c("Best" = "#85B068",
                               "Simple" = "#0F4C5E")) +
  scale_color_manual(values = c("Best" = "#56783F",
                                "Simple" = "#0F4C5E")) +
  guides(color = "none",
         fill = guide_legend(override.aes = list(alpha = 0.7,
                                                 color = c("#0F4C5E", "#56783F")))) +
  labs(x = "Chance of Negative (%)", y = "Density",
       fill = "Model", tag = "B") + 
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0),
                     breaks = seq(0, 1, 0.25), 
                     labels = c("0","25", "50", "75", "")) +
  scale_y_continuous(limits = c(0, 7), breaks = seq(0, 6, 1),
                     expand = c(0, 0),
                     position = "left") +
  theme(legend.position = c(0.22, 0.7), 
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        legend.background = element_rect(fill = 'white'),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 5, 0, 0)); fig3B.1

# Plot all observed sgRNA positives
fig3B.2 <- ggplot(subset(prediction.df, posneg == "Positive")) +
  geom_vline(xintercept = threshold, size = 0.7, color = zissou_pal[1], alpha = 0.8,
             linetype = "dashed") +
  geom_density(aes(x = (as.numeric(prob_mean)),
                   fill = factor(model, levels = c("Simple", "Best")),
                   color = factor(model, levels = c("Simple", "Best"))),
               alpha = 0.4) +
  facet_grid(.~ "Observed Positive") +
  scale_fill_manual(values = c("Best" = "#85B068",
                               "Simple" = "#0F4C5E")) +
  scale_color_manual(values = c("Best" = "#56783F",
                                "Simple" = "#0F4C5E")) +
  guides(color = "none",
         fill = guide_legend(override.aes = list(alpha = 0.6,
                                                 color = c( "#0F4C5E", "#56783F")))) +
  labs(x = "Chance of Positive (%)", y = "Density",
       fill = "Model") + 
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0),
                     breaks = seq(0, 1, 0.25), 
                     labels = c("0","25", "50", "75", "100")) +
  scale_y_continuous(limits = c(0, 7), breaks = seq(0, 8, 1),
                     expand = c(0, 0),
                     position = "left") +
  theme(legend.position = "none", 
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        legend.background = element_rect(fill = 'white'),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 0, 0, 0)); fig3B.2

fig.perf <- fig3B.1 + fig3B.2 + plot_layout(ncol = 2); fig.perf





## 3AB: Top-left panels  -----------------------------------------------------

fig3AB <- fig.data + fig.perf + plot_layout(ncol = 2, widths = c(1, 0.8)); fig3AB


## 3D: DOSE ----------------------------------------------------------------------

x_vals <- seq(from = 0, to = 10, by = 0.05)
dose <- 5.5
set.seed(11)

### Prep simple model results ---------------------------------------------------

# Extract relevant parameter estimates
pars_simple <- c("gamma", "deltaT", "alpha", "betaT")
param.simple.df.wide <- fit.simple$draws(variables = pars_simple,
                                         inc_warmup = FALSE,
                                         format = "df")
param.simple.df <- param.simple.df.wide %>% pivot_longer(cols = everything(),
                                                         names_to = "param",
                                                         values_to = "value")
param.simple.df <- as.data.frame(param.simple.df)
param.simple.df <- subset(param.simple.df, param %notin% c(".chain", ".iteration", ".draw"))
param.simple.df$param <- factor(param.simple.df$param, levels = rev(pars_simple))

# Extract the mean for the reference line
delta_mean <- mean(param.simple.df.wide$deltaT)
gamma_mean <- mean(param.simple.df.wide$gamma)

# Generate values & a data frame for the reference line
y_vals <- exp(gamma_mean + delta_mean * x_vals)/(1 + exp(gamma_mean + delta_mean * x_vals))
lines.simple.df <- as.data.frame(cbind(x_vals, y_vals))


### Prep best model results -----------------------------------------------------

# Extract parameter estimates from best model
pars_best <- c("gamma", "deltaT", 
               "deltaDOSE",  
               "deltaTG[1]", "deltaTG[2]", "deltaTG[3]", "deltaTG[4]",
               "deltaSP[1]", "deltaSP[2]", "deltaSP[3]")
param.best.df.wide <- fit.best$draws(variables = pars_best,
                                   inc_warmup = FALSE,
                                   format = "df")
param.best.df <- param.best.df.wide %>% pivot_longer(cols = everything(),
                                                     names_to = "param",
                                                     values_to = "value")
param.best.df <- as.data.frame(param.best.df)
param.best.df <- subset(param.best.df, param %notin% c(".chain", ".iteration", ".draw"))
param.best.df$param <- factor(param.best.df$param, levels = rev(pars_best))


# Extract lines for different dose values, holding all else constant
dose.vary <- sample(runif(n_lines, min(dat$log10_dose_pfu), max(dat$log10_dose_pfu)))

lines.best.df <- as.data.frame(matrix(nrow = 0, ncol = 4))


for (ii in 1:length(dose.vary)){
  row <- sample(1:nrow(param.best.df.wide), size = 1)
  y <- param.best.df.wide$gamma[row] +
       param.best.df.wide$deltaT[row] * x_vals +
       param.best.df.wide$'deltaTG[3]'[row] +
       param.best.df.wide$deltaDOSE[row] * dose.vary[ii] +
       param.best.df.wide$'deltaSP[1]'[row]
  y_vals <- exp(y)/(1 + exp(y))
  line_num <- rep(ii, length(x_vals))
  log10_dose <- rep(dose.vary[ii], length(x_vals))
  
  lines.best.df <- rbind(lines.best.df,
                    cbind(line_num, log10_dose, x_vals, y_vals))
}


### Plot ------------------------------------------------------------------------

fig.dose.fits <- ggplot() + 
  geom_line(data = lines.best.df,
            aes(x = x_vals, y = y_vals, 
                group = line_num,
                color = log10_dose),
            size = 0.8, alpha = 0.5) +
  geom_line(data = as.data.frame(lines.simple.df),
            aes(x = x_vals, y = y_vals),
            color = "#0F4C5E",
            size = 1) +
  scale_color_gradientn(colours = zissou_pal) +
  scale_fill_gradientn(colours = zissou_pal) +
  labs(y = "Chance of Positive (%)", x = "",
       color = "log10\ndose", 
       tag = "D") +
  guides(fill = "none",
         color = guide_colorbar(frame.colour = "black",
                                #title.hjust = -1,
                                title.vjust = 0)) +
  theme(legend.position = c(0.85, 0.38), 
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, 'cm'),
    legend.key = element_blank(),
    axis.title.x = element_blank(),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    text = element_text(size = 12),
    panel.background = element_rect(fill = "white",
                                    colour = "black",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                    colour = "light grey"), 
    panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                    colour = "light grey"),
    plot.margin = margin(0, 0, 0, 0)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1),
                     labels = c(0, "20", "40", "60", "80", "100"), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10), limits = c(-0.01, 10),
                     labels = c(0, 2, 4, 6, 8, ""),
                     expand = c(0, 0)); fig.dose.fits


## 3D: SPECIES -----------------------------------------------------------------

# Generate lines for different target genes, holding all else constant
sp.vary <- sample(rep(1:3, each = n_lines / 3))
set.seed(1111)
lines.sp.df <- as.data.frame(matrix(nrow = 0, ncol = 4))


for (ii in 1:length(sp.vary)){
  row <- sample(1:nrow(param.best.df.wide), size = 1)
  
  if (sp.vary[ii] == 1){
    sp_sample <- param.best.df.wide$'deltaSP[1]'[row]
  }
  else if (sp.vary[ii] == 2){
    sp_sample <- param.best.df.wide$'deltaSP[2]'[row]
  }
  else if (sp.vary[ii] == 3){
    sp_sample <- param.best.df.wide$'deltaSP[3]'[row]
  }
  
  y <- param.best.df.wide$gamma[row] +
    param.best.df.wide$deltaT[row] * x_vals +
    param.best.df.wide$deltaDOSE[row] * dose +
    param.best.df.wide$'deltaTG[3]'[row] +
    sp_sample
  
  y_vals <- exp(y)/(1 + exp(y))
  line_num <- rep(ii, length(x_vals))
  sp_idx <- rep(sp.vary[ii], length(x_vals))
  
  lines.sp.df <- rbind(lines.sp.df,
                       cbind(line_num, sp_idx, x_vals, y_vals))
}

lines.sp.df$sp_idx[lines.sp.df$sp_idx == 1] <- "RM"
lines.sp.df$sp_idx[lines.sp.df$sp_idx == 2] <- "CM"
lines.sp.df$sp_idx[lines.sp.df$sp_idx == 3] <- "AGM"


fig.sp.fits <- ggplot() + 
  geom_line(data = lines.sp.df,
            aes(x = x_vals, y = y_vals, 
                group = line_num,
                color = as.character(sp_idx)),
            size = 0.8, alpha = 0.1) +
  geom_line(data = as.data.frame(lines.simple.df),
            aes(x = x_vals, y = y_vals),
            color = "#002C3D",
            size = 0.8) +
  scale_color_manual(values = c("RM" = "#0B6884",
                                "CM" = "#744468",
                                "AGM" = "#D56540")) +
  labs(y = element_blank(), x = element_blank(),
       color = "Species") +
  guides(color = guide_legend(override.aes = list(size = 2.7, alpha = 1,
                                                  linewidth = 2), 
                              keyheight = 0.6)) + 
  theme(legend.position = c(0.85, 0.252), 
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "white"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(0, 0, 0, 0)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1),
                     labels = c(0, "20", "40", "60", "80", "100"), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10), limits = c(-0.01, 10),
                     labels = c(0, 2, 4, 6, 8, ""),
                     expand = c(0, 0)); fig.sp.fits



## 3D: TARGET GENE -------------------------------------------------------------

# Generate lines for different target genes, holding all else constant
tg.vary <- sample(rep(1:4, each = n_lines / 4))
set.seed(11111)
lines.tg.df <- as.data.frame(matrix(nrow = 0, ncol = 4))

for (ii in 1:length(tg.vary)){
  row <- sample(1:nrow(param.best.df.wide), size = 1)
  
  if (tg.vary[ii] == 1){
    tg_sample <- param.best.df.wide$'deltaTG[1]'[row]
  }
  else if (tg.vary[ii] == 2){
    tg_sample <- param.best.df.wide$'deltaTG[2]'[row]
  }
  else if (tg.vary[ii] == 3){
    tg_sample <- param.best.df.wide$'deltaTG[3]'[row]
  }
  else if (tg.vary[ii] == 4){
    tg_sample <- param.best.df.wide$'deltaTG[4]'[row]
  }
  
  y <- param.best.df.wide$gamma[row] +
    param.best.df.wide$deltaT[row] * x_vals +
    param.best.df.wide$deltaDOSE[row] * dose +
    tg_sample + 
    param.best.df.wide$'deltaSP[1]'[row] 
  
  y_vals <- exp(y)/(1 + exp(y))
  line_num <- rep(ii, length(x_vals))
  tg_idx <- rep(tg.vary[ii], length(x_vals))
  
  lines.tg.df <- rbind(lines.tg.df,
                       cbind(line_num, tg_idx, x_vals, y_vals))
}

lines.tg.df$tg_idx[lines.tg.df$tg_idx == 1] <- "T↑ SG↑"
lines.tg.df$tg_idx[lines.tg.df$tg_idx == 2] <- "T↓ SG↑"
lines.tg.df$tg_idx[lines.tg.df$tg_idx == 3] <- "T↑ SG↓"
lines.tg.df$tg_idx[lines.tg.df$tg_idx == 4] <- "T↓ SG↓"


fig.tg.fits <- ggplot() + 
  geom_line(data = lines.tg.df,
            aes(x = x_vals, y = y_vals, 
                group = line_num,
                color = as.character(tg_idx)),
            size = 0.8, alpha = 0.1) +
  geom_line(data = as.data.frame(lines.simple.df),
            aes(x = x_vals, y = y_vals),
            color = "#002C3D",
            size = 0.8) +
  scale_color_manual(values = c("T↑ SG↑" = "#456187" , 
                                "T↓ SG↑" =  "#F39237", #"#E89005", 
                                "T↑ SG↓" = "#52A4BC" ,
                                "T↓ SG↓" = "#A1C084"),
                     breaks = c("T↑ SG↑", 
                                "T↓ SG↑",
                                "T↑ SG↓",
                                "T↓ SG↓")) +
  labs(y = element_blank(), x = element_blank(),
       color = "Target Gene") +
  guides(color = guide_legend(override.aes = list(size = 2.7, alpha = 1,
                                                  linewidth = 2), 
                              keyheight = 0.6)) + 
  theme(legend.position = c(0.8, 0.29), 
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(0, 0, 0, 0)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1),
                     labels = c(0, "20", "40", "60", "80", "100"), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10), limits = c(-0.01, 10),
                     labels = c(0, 2, 4, 6, 8, ""),
                     expand = c(0, 0)); fig.tg.fits






# 3D: Combine bottom left plots  ------------------------------------------------

fig3D <- fig.dose.fits +
  fig.sp.fits + 
  fig.tg.fits + 
  plot_layout(ncol = 3); fig3D

x_title <- ggplot() + 
  annotate(geom = "text", x = 1, y = 1.7, size = 2.8,
           label = "log10 total RNA copies / sample") +
  theme(plot.margin = c(0, -10, -10, -10)) +
  scale_y_continuous(limits = c(0, 2)) +
  theme_void()


fig3D <- ggarrange(fig3D, x_title, nrow = 2, heights = c(1, 0.1)); fig3D


# 3ABD: Combine left plots  ----------------------------------------------------

fig.left <- ggarrange(fig3AB, fig3D, nrow = 2, heights = c(0.9, 1)); fig.left


# 3C: Heat maps ----------------------------------------------------------------

x_vals <- seq(from = 0, to = 10, by = 1)

# Generate values & a data frame for the reference line
y_vals <- exp(gamma_mean + delta_mean * x_vals)/(1 + exp(gamma_mean + delta_mean * x_vals))
lines.simple.df <- as.data.frame(cbind(x_vals, y_vals))
lines.simple.df$predictor <- ""


## Determine range sgRNA has different outcome ---------------------------------

range.sg.df <- data.frame(dose = NA, sp = NA, tg = NA, 
                          x_val = NA, median_prob = NA, quantile_probs = NA,
                          simple_prob = NA)

for (dose.option in c(4, 5.5, 7)) {
  
  for (sp.option in c(1, 2, 3)){
    if (sp.option == 1){sp.sample <-  param.best.df.wide$'deltaSP[1]'}
    else if(sp.option == 2){sp.sample <- param.best.df.wide$'deltaSP[2]'}
    else if(sp.option == 3){sp.sample <- param.best.df.wide$'deltaSP[3]'}
    
    for (tg.option in c(1, 2, 3, 4)) {
      if (tg.option == 1){tg.sample <-  param.best.df.wide$'deltaTG[1]'}
      else if(tg.option == 2){tg.sample <- param.best.df.wide$'deltaTG[2]'}
      else if(tg.option == 3){tg.sample <- param.best.df.wide$'deltaTG[3]'}
      else if(tg.option == 4){tg.sample <- param.best.df.wide$'deltaTG[4]'}
      
     
       for (x_val in x_vals) {
        y <- param.best.df.wide$gamma +
          param.best.df.wide$deltaT * x_val +
          sp.sample +
          tg.sample +
          param.best.df.wide$deltaDOSE * dose.option 
        y_prob <- exp(y)/(1+exp(y))
        y_med <- median(y_prob) 
        y_quant <- paste0(round(quantile(y_prob, c(0.05, 0.95)), 2), collapse = ", ")
        
        range.sg.df <- rbind(range.sg.df,
                             data.frame(
                               dose = dose.option,
                               sp = sp.option,
                               tg = tg.option,
                               x_val = x_val, 
                               median_prob = y_med, 
                               quantile_probs = y_quant,
                               simple_prob = lines.simple.df$y_vals[lines.simple.df$x_vals == x_val]))
        
      }
    }
  }
}

# Get the relevant ranges
sort(unique(range.sg.df$x_val[(range.sg.df$median_prob > 0.5 & range.sg.df$simple_prob <= 0.5) |
                               range.sg.df$median_prob  <= 0.5 & range.sg.df$simple_prob > 0.5]))

range.sg.df <- range.sg.df[-1, ]


## Make the Supplemental Table for 90% CrI -------------------------------------

tblS7 <- subset(range.sg.df, select = -c(median_prob, simple_prob))

tblS7 <- tblS7 %>%
  pivot_wider(names_from = x_val, values_from = c(quantile_probs))

tblS7.doses <- subset(tblS7, dose %in% c(4, 5.5, 7) & sp == 1 & tg == 3) %>%
  rename(predictor = dose) %>%
  select(-c(sp, tg)); tblS7.doses

tblS7.species <- subset(tblS7, dose == 5.5 & tg == 3) %>%
  rename(predictor = sp) %>%
  select(-c(dose, tg)); tblS7.species

tblS7.tg <- subset(tblS7, dose == 5.5 & sp == 1 & tg %in% c(1, 2, 3, 4)) %>%
  rename(predictor = tg) %>%
  select(-c(dose, sp)); tblS7.tg

tblS7.combined <- rbind(tblS7.doses, tblS7.species, tblS7.tg); tblS7.combined
tblS7.combined$group <- c(rep("Dose", 3),
                          rep("Species", 3),
                          rep("Target Gene", 4))
tblS7.combined <- tblS7.combined %>%
  select(c(group, predictor, everything())); tblS7.combined


write.csv(tblS7.combined, file = "./outputs/tables/tblS6A-sgRNA-logistic-90CrI.csv",
          row.names = FALSE)





## 3C: Simple predictions ---------------------------------------------------------

fig.simple <- ggplot(as.data.frame(lines.simple.df), 
                     aes(x = as.numeric(x_vals), 
                         y = predictor, 
                         fill = as.numeric(y_vals))) +
  geom_tile() + 
  geom_text(aes(label = (round(as.numeric(y_vals), 2)) * 100), 
            color = "black", size = 2.5) +
  geom_rect(aes(xmin = 1.5, xmax = 6.5, ymin = 0.5, ymax = 1.5), 
            fill = "transparent", color = "grey33", size = 1.5) +
  scale_fill_distiller(palette = "RdBu") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0.5, 1.5), clip = "on") +
  labs(x = "log10 total RNA copies / sample",
       y = "Simple") +
  theme(legend.position = "none",
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 7),
        axis.text.y = element_text(angle = 0, 
                                   vjust = 1, hjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(0, 0, 0.75, 0)) + 
  #scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12), expand = c(0, 0))
  scale_y_discrete(expand = c(0, 0)); fig.simple


## 3C: DOSE --------------------------------------------------------------------

# Extract parameter estimates from best model
pars_best <- c("gamma", "deltaT", 
               "deltaDOSE", 
               "deltaTG[1]", "deltaTG[2]", "deltaTG[3]", "deltaTG[4]",
               "deltaSP[1]", "deltaSP[2]", "deltaSP[3]")
param.best.df.wide <- fit.best$draws(variables = pars_best,
                                     inc_warmup = FALSE,
                                     format = "df")
param.best.df <- param.best.df.wide %>% pivot_longer(cols = everything(),
                                                     names_to = "param",
                                                     values_to = "value")
param.best.df <- as.data.frame(param.best.df)
param.best.df <- subset(param.best.df, param %notin% c(".chain", ".iteration", ".draw"))
param.best.df$param <- factor(param.best.df$param, levels = rev(pars_best))


dose.df <- data.frame(tRNA = NA,
                      dose = NA,
                      probs = NA)

for (x_val in x_vals) {
  y <- param.best.df.wide$gamma +
    param.best.df.wide$deltaT * x_val +
    param.best.df.wide$'deltaTG[3]' +
    param.best.df.wide$'deltaSP[1]' +
    param.best.df.wide$deltaDOSE * 4 
  
  y_small <- median(exp(y)/(1+exp(y)))
  
  y <- param.best.df.wide$gamma +
    param.best.df.wide$deltaT * x_val +
    param.best.df.wide$'deltaTG[3]' +
    param.best.df.wide$'deltaSP[1]' +
    param.best.df.wide$deltaDOSE * 5.5 
  
  y_mid <- median(exp(y)/(1+exp(y)))
  
  y <- param.best.df.wide$gamma +
    param.best.df.wide$deltaT * x_val +
    param.best.df.wide$'deltaTG[3]' +
    param.best.df.wide$'deltaSP[1]' +
    param.best.df.wide$deltaDOSE * 7 
  
  y_big <- median(exp(y)/(1+exp(y)))
  
  dose.df <- rbind(dose.df, 
                   data.frame(tRNA = x_val, 
                              dose = "10^4", 
                              probs = y_small),
                   data.frame(tRNA = x_val, 
                              dose = "10^5.5",
                              probs = y_mid),
                   data.frame(tRNA = x_val, 
                              dose = "10^7",
                              probs = y_big))
}

dose.df <- dose.df[-1, ]

fig.dose <- ggplot(dose.df, aes(x = as.numeric(tRNA), 
                                y = factor(dose, levels = c("10^7", "10^5.5", "10^4")), 
                                fill = as.numeric(probs))) +
  geom_tile() + 
  geom_text(aes(label = (round(as.numeric(probs), 2)) * 100), color = "black", size = 2.5) +
  geom_rect(aes(xmin = 3.5, xmax = 5.5, ymin = 0.5, ymax = 3.5), 
            fill = "transparent", color = "grey33", size = 1.5) +
  scale_fill_distiller(palette = "RdBu") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0.5, 3.5), clip = "on") +
  labs(x = "log10 total RNA copies / sample", 
       y = "log10 Dose") + 
  scale_y_discrete(expand = c(0, 0), labels = c(7, 5.5, 4)) +
  theme(legend.position = "none",
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 7),
        axis.text.y = element_text(angle = 45, size = 6.5, face = c("plain", "bold.italic", "plain"), 
                                   vjust = 1, hjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(0, 0, 0.75, 0)) ; fig.dose



## 3C: SPECIES -----------------------------------------------------------------

sp.df <- data.frame(tRNA = NA,
                    sp = NA,
                    probs = NA)

for (x_val in x_vals) {
  y <- param.best.df.wide$gamma +
    param.best.df.wide$deltaT * x_val +
    param.best.df.wide$'deltaTG[3]' +
    param.best.df.wide$'deltaSP[1]' +
    param.best.df.wide$deltaDOSE * 5.5 
  
  y_rm <- median(exp(y)/(1+exp(y)))
  
  y <- param.best.df.wide$gamma +
    param.best.df.wide$deltaT * x_val +
    param.best.df.wide$'deltaTG[3]' +
    param.best.df.wide$'deltaSP[2]' +
    param.best.df.wide$deltaDOSE * 5.5 
  
  y_cm <- median(exp(y)/(1+exp(y)))
  
  y <- param.best.df.wide$gamma +
    param.best.df.wide$deltaT * x_val +
    param.best.df.wide$'deltaTG[3]' +
    param.best.df.wide$'deltaSP[3]' +
    param.best.df.wide$deltaDOSE * 5.5 
  
  y_agm <- median(exp(y)/(1+exp(y)))
  
  sp.df <- rbind(sp.df, 
                  data.frame(tRNA = x_val, 
                             sp = "RM", 
                             probs = y_rm),
                  data.frame(tRNA = x_val, 
                             sp = "CM",
                             probs = y_cm),
                  data.frame(tRNA = x_val, 
                             sp = "AGM",
                             probs = y_agm))
}

sp.df <- sp.df[-1, ]

fig.sp <- ggplot(sp.df, aes(x = as.numeric(tRNA), 
                            y = factor(sp, levels = rev(c("RM", "CM", "AGM"))),
                            fill = as.numeric(probs))) +
  geom_tile() + 
  geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = 0.5, ymax = 3.5), 
            fill = "transparent", color = "grey33", size = 1.5) +
  geom_text(aes(label = (round(as.numeric(probs), 2)) * 100), color = "black", size = 2.5) +
  scale_fill_distiller(palette = "RdBu", breaks = c(0, 0.5, 1), 
                       labels = c(0, "50", "100"), limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0.5, 3.5), clip = "on") +
  labs(x = "log10 total RNA copies / sample", fill = "Chance of Positive (%)",
       y = "Species") +
  guides(fill = guide_colorbar(frame.colour = "black")) +
  theme(legend.position = "none",
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 7),
        axis.text.y = element_text(angle = 45, size = 6.5,
                                   face = c("plain", "plain", "bold"), 
                                   vjust = 1, hjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(0, 0, 0.75, 0)) + 
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(-0.5, 10, 2), labels = seq(0, 10, 2)); fig.sp


## 3C: TARGET GENE -------------------------------------------------------------

tg.df <- data.frame(tRNA = NA,
                    tg = NA,
                    probs = NA)

for (x_val in x_vals) {
  y <- param.best.df.wide$gamma +
    param.best.df.wide$deltaT * x_val +
    param.best.df.wide$'deltaSP[1]' +
    param.best.df.wide$'deltaTG[1]' +
    param.best.df.wide$deltaDOSE * 5.5 
  
  y_tg1 <- median(exp(y)/(1+exp(y)))
  quant_tg1 <- quantile(exp(y)/(1+exp(y)), 0.05, 0.95)
  
  y <- param.best.df.wide$gamma +
    param.best.df.wide$deltaT * x_val +
    param.best.df.wide$'deltaSP[1]' +
    param.best.df.wide$'deltaTG[2]' +
    param.best.df.wide$deltaDOSE * 5.5 
  
  y_tg2 <- median(exp(y)/(1+exp(y)))
  quant_tg2 <- quantile(exp(y)/(1+exp(y)), 0.05, 0.95)
  
  y <- param.best.df.wide$gamma +
    param.best.df.wide$deltaT * x_val +
    param.best.df.wide$'deltaSP[1]' +
    param.best.df.wide$'deltaTG[3]' +
    param.best.df.wide$deltaDOSE * 5.5 
  
  y_tg3 <- median(exp(y)/(1+exp(y)))
  quant_tg3 <- quantile(exp(y)/(1+exp(y)), 0.05, 0.95)
  
  y <- param.best.df.wide$gamma +
    param.best.df.wide$deltaT * x_val +
    param.best.df.wide$'deltaSP[1]' +
    param.best.df.wide$'deltaTG[4]' +
    param.best.df.wide$deltaDOSE * 5.5 
  
  y_tg4 <- median(exp(y)/(1+exp(y)))
  quant_tg4 <- quantile(exp(y)/(1+exp(y)), 0.05, 0.95)
  
  tg.df <- rbind(tg.df, 
                 data.frame(tRNA = x_val, 
                            tg = "T↑ SG↑", 
                            probs = y_tg1),
                 data.frame(tRNA = x_val, 
                            tg = "T↓ SG↑",
                            probs = y_tg2),
                 data.frame(tRNA = x_val, 
                            tg = "T↑ SG↓",
                            probs = y_tg3),
                 data.frame(tRNA = x_val, 
                            tg = "T↓ SG↓",
                            probs = y_tg4))
}

tg.df <- tg.df[-1, ]

fig.tg <- ggplot(tg.df, aes(x = as.numeric(tRNA), 
                            y = factor(tg, levels = rev(c("T↑ SG↑", 
                                                          "T↓ SG↑",
                                                          "T↑ SG↓",
                                                          "T↓ SG↓"))),
                            fill = as.numeric(probs))) +
  geom_tile() + 
  geom_rect(aes(xmin = 2.5, xmax = 4.5, ymin = 0.5, ymax = 4.5), 
            fill = "transparent", color = "grey33", size = 1.5) +
  geom_text(aes(label = (round(as.numeric(probs), 2)) * 100), color = "black", size = 2.5) +
  scale_fill_distiller(palette = "RdBu", breaks = c(0, 0.5, 1), 
                       labels = c(0, "50", "100"), limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0.5, 4.5), clip = "on") +
  labs(x = "log10 total RNA copies / sample", fill = "Chance of Positive (%)",
       y = "Target Gene") +
  guides(fill = guide_colorbar(frame.colour = "black")) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 7.5, vjust = 0.7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        legend.margin=margin(t=-1),
        axis.title = element_text(size = 7),
        axis.text.y = element_text(angle = 45, size = 6.5,
                                   face = c("plain", "bold", "plain",  "plain"), 
                                   vjust = 1, hjust = 0.5),
        axis.text.x = element_text(size = 6.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(0, 0, 0.75, 0)) + 
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(-0.5, 10, 2), labels = seq(0, 10, 2)); fig.tg



## 3C: Combine -----------------------------------------------------------------

fig.right <- fig.simple + labs(tag = "C") + 
  fig.dose + 
  fig.sp + 
  fig.tg + 
  plot_layout(ncol = 1, heights = c(1, 3, 3, 4)); fig.right



# 3A-D: All Logistic Panels ----------------------------------------------------

fig.log <- ggarrange(fig.left, fig.right, ncol = 2, widths = c(1, 0.4))

# Generate text-only figure for joint labels
log_label <- ggplot() + 
  annotate(geom = "text", x = 1, y = 1, size = 3.2, 
           label = "Logistic Component", angle = -90,
           fontface = 2) +
  geom_segment(aes(x = 0.6, xend = 0.6, y = -0.5, yend = 2.5), color = "grey42") +
  theme(plot.margin = c(-10, -10, -10, -10)) +
  scale_x_continuous(limits = c(0, 2)) +
  theme_void(); log_label


fig.log <- ggarrange(fig.log, log_label, ncol = 2, widths = c(1, 0.04)); fig.log


# LINEAR COMPONENT PANELS ------------------------------------------------------

## Add jitter to < LOD values for visualization --------------------------------

dat$val_sg_jitter <- dat$val_sg
dat$val_total_jitter <- dat$val_total
for (row_num in 1:nrow(dat)){
  if (dat$pos_sg[row_num] == 0){
    dat$val_sg_jitter[row_num] <- runif(1, -1.8, -0.2)
  }
  if (dat$pos_total[row_num] == 0){
    dat$val_total_jitter[row_num] <- runif(1, -1.8, -0.2)
  }
}


## 3E: Raw data ----------------------------------------------------------------

dat.quant <- subset(dat, cens_sg != 3 & cens_total != 3)
dat.quant$combo_idx[dat.quant$pos_sg == 1 & dat.quant$pos_total == 1] <- "Both Positive"
dat.quant$combo_idx[dat.quant$pos_sg == 0 & dat.quant$pos_total == 1] <- "Only totRNA Positive"
dat.quant$combo_idx[dat.quant$pos_sg == 1 & dat.quant$pos_total == 0] <- "Only sgRNA Positive"
dat.quant$combo_idx[dat.quant$pos_sg == 0 & dat.quant$pos_total == 0] <- "Both Negative"


fig.data.quant <- ggplot(subset(dat.quant, cens_total == 0)) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "grey55", alpha = 0.2) +
  annotate("rect", ymin = -Inf, ymax = Inf, xmin = -Inf, xmax = 0,
           fill = "grey55", alpha = 0.2) +
  annotate("text", label = "< LOD", x = 8.8, y = -0.5, size = 2.7) +
  geom_vline(xintercept = 0, linewidth = 0.5, color = "black")+
  geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +

  geom_point(aes(x = val_total_jitter, y = val_sg_jitter, 
                 fill = combo_idx),
             shape = 21, size = 1.3, alpha = 0.6) +
  labs(y = "log10 sgRNA copies / sample", 
       x = "log10 total RNA copies / sample",
       fill = element_blank(),
       tag = "E") +
  scale_fill_manual(values = c("Both Positive" = brewer_pal[2],
                               "Only totRNA Positive" = brewer_pal[9],
                               "Only sgRNA Positive" = brewer_pal[2],
                               "Both Negative" = brewer_pal[9])) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), limits = c(-2, 10),
                     labels = c(0, 2, 4, 6, 8, 10), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10), limits = c(0, 10),
                     expand = c(0, 0)) +
  theme(legend.position = "none", 
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.key=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        text = element_text(size=12),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        #panel.grid.minor = element_blank(),
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(4, 4, 2, 2)); fig.data.quant
  

## 3F: Median error distributions ----------------------------------------------

dat <- read.csv("./data/clean-data.csv")

# Removes positively censored tRNA samples 
dat <- subset(dat, pos_sg %in% c(0, 1) & cens_total == 0) 

# Set sgRNA positive but censored values to -9 (to screen out from linear component)
dat$val_sg[dat$cens_sg == 3] <- -9

# Convert columns to numerics
dat$val_sg <- as.numeric(dat$val_sg)
dat$val_total <- as.numeric(dat$val_total)

dat <- subset(dat, cens_total == 0)

### Prep simple model results --------------------------------------------------

pred_val <- paste0("pred_val[", which(dat$cens_sg == 0), "]") 
pred_df <- as.data.frame(fit.simple$draws(pred_val, format = "df"))
pred_df <- subset(pred_df, select = pred_val) #remove .chain, ...

# Add ME/MAE and ME scaled terms 
dat.obs <- subset(dat, cens_sg == 0)
pred_median <- sapply(pred_df, function(x) median(x)) # CHANGE BACK TO MEDIAN ONCE UPDATED
pred_sd <- sapply(pred_df, function(x) sd(x))
me <- dat.obs$val_sg - pred_median
mae <- abs(dat.obs$val_sg - pred_median)
mae_scaled <- abs(dat.obs$val_sg - pred_median)/pred_sd
pred.simple.df <- as.data.frame(matrix(data = NA, nrow = length(me),
                                       ncol = 2))
pred.simple.df[, 1] <- "Simple"
pred.simple.df[, 2] <- me



### Prep best model results -----------------------------------------------------

pred_val <- paste0("pred_val[", which(dat$cens_sg == 0), "]") 
pred_df <- as.data.frame(fit.best$draws(pred_val, format = "df"))
pred_df <- subset(pred_df, select = pred_val) #remove .chain, ...

# Add ME/MAE and ME scaled terms 
dat.obs <- subset(dat, cens_sg == 0)
pred_median <- sapply(pred_df, function(x) median(x)) # CHANGE BACK TO MEDIAN ONCE UPDATED
pred_sd <- sapply(pred_df, function(x) sd(x))
me <- dat.obs$val_sg - pred_median
mae <- abs(dat.obs$val_sg - pred_median)
mae_scaled <- abs(dat.obs$val_sg - pred_median)/pred_sd
pred.df <- as.data.frame(matrix(data = NA, nrow = length(me),
                                ncol = 2))
pred.df[, 1] <- "Best"
pred.df[, 2] <- me

pred.df <- rbind(pred.df, pred.simple.df)


### Plot -----------------------------------------------------------------------

mae_best <- median(abs(subset(pred.df, V1 == "Best")$V2))
mae_simple <- median(abs(subset(pred.df, V1 == "Simple")$V2))

fig3F <- ggplot(pred.df) +
  geom_density(aes(x = abs(V2), fill = factor(V1, levels = c("Simple", "Best")),
                   color = factor(V1, levels = c("Simple", "Best"))),
               alpha = 0.4) +
  geom_density(aes(x = abs(V2),
                   color = factor(V1, levels = c("Simple", "Best"))),
               alpha = 1) +
  geom_vline(xintercept = mae_best, size = 0.7, 
             color = "#56783F", alpha = 1,
             linetype = "dotted") +
  geom_vline(xintercept = mae_simple, size = 0.7, 
             color = "#0F4C5E", alpha = 1,
             linetype = "dotted") +
  #facet_grid(V1 ~.) +
  scale_fill_manual(values = c("Best" = "#85B068",
                               "Simple" = "#0F4C5E")) +
  scale_color_manual(values = c("Best" = "#56783F",
                                "Simple" = "#0F4C5E")) +
  guides(color = "none",
         fill = guide_legend(override.aes = list(alpha = 0.6,
                                                 color = c( "#0F4C5E", "#56783F")))) +
  labs(x = "Median absolute prediction error (log10)", y = "Density",
       fill = "Model", tag = "F") +
  theme(legend.position = c(0.87, 0.60), 
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        legend.background = element_rect(fill = 'white'),
        text = element_text(size = 11),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor = element_blank()) + 
  scale_x_continuous(limits = c(0, 3), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.25), breaks = seq(0, 1.25, 0.5),
                     expand = c(0, 0),
                     position = "left"); fig3F


##3EF: Combine -----------------------------------------------------------------

fig3EF <- fig.data.quant + fig3F + plot_layout(nrow = 2, heights = c(1, 0.3)); fig3EF


## 3G: Fit lines ---------------------------------------------------------------

### Prep simple model results --------------------------------------------------

pars_simple <- c("alpha", "betaT", "sigma")
param.simple.df <- fit.simple$draws(variables = pars_simple,
                                    inc_warmup = FALSE,
                                    format = "df")
param.simple.df <- as.data.frame(param.simple.df)

# Generate values for the line
beta_mean <- mean(param.simple.df$betaT)
alpha_mean <- mean(param.simple.df$alpha)

x_vals <- seq(from = 0, to = 10, by = 0.05)
y_vals <- alpha_mean + beta_mean * x_vals

# Combine into df for plotting
lines.simple.df <- as.data.frame(cbind(x_vals, y_vals))


### Prep best model results -----------------------------------------------------

# Extract parameter estimates
pars_best <- c("alpha", "betaT", 
               "betaTG[1]", "betaTG[2]", "betaTG[3]", "betaTG[4]",
               "betaSP[1]", "betaSP[2]", "betaSP[3]",
               "betaDPI[1]", "betaDPI[2]", "betaDPI[3]",
               "betaDOSE",
               "sigma_bar",
               "sigma_sd")
param.best.df <- fit.best$draws(variables = pars_best,
                                inc_warmup = FALSE,
                                format = "df")
param.best.df <- as.data.frame(param.best.df)


### 3G: TARGET GENE -----------------------------------------------------

# Generate lines for different target genes, holding all else constant
# There are no samples with quantitative sgRNA data that classify as TG1 (Both high)
#    So we do not plot them

tg <- sample(rep(2:4, each = n_lines / 3))
set.seed(1)
lines.df <- as.data.frame(matrix(nrow = 0, ncol = 4))

for (ii in 1:length(tg)) {
  row <- sample(1:nrow(param.best.df), size = 1)
  
  if (tg[ii] == 1){
    tg_sample <- param.best.df$'betaTG[1]'[row]
  }
  else if (tg[ii] == 2){
    tg_sample <- param.best.df$'betaTG[2]'[row]
  }
  else if (tg[ii] == 3){
    tg_sample <- param.best.df$'betaTG[3]'[row]
  }
  else if (tg[ii] == 4){
    tg_sample <- param.best.df$'betaTG[4]'[row]
  }
  
  y_vals <- param.best.df$alpha[row] +
    param.best.df$betaT[row] * x_vals +
    tg_sample +
    param.best.df$'betaSP[1]'[row] + 
    param.best.df$'betaDPI[2]'[row] +
    param.best.df$betaDOSE[row] * 5.5
  
  line_num <- rep(ii, length(x_vals))
  tg_idx <- rep(tg[ii], length(x_vals))
  
  lines.df <- rbind(lines.df,
                    cbind(line_num, tg_idx, x_vals, y_vals))
}

lines.df$tg_idx[lines.df$tg_idx == 1] <- "T↑ SG↑"
lines.df$tg_idx[lines.df$tg_idx == 2] <- "T↓ SG↑"
lines.df$tg_idx[lines.df$tg_idx == 3] <- "T↑ SG↓"
lines.df$tg_idx[lines.df$tg_idx == 4] <- "T↓ SG↓"



fig.quant.tg.fits <- ggplot() + 
  geom_line(data = lines.df,
            aes(x = x_vals, y = y_vals, 
                group = line_num,
                color = tg_idx),
            size = 0.8, alpha = 0.1) +
  geom_line(data = lines.simple.df,
            aes(x = x_vals, y = y_vals),
            color = "grey25", size = 1) +
  scale_color_manual(values = c("T↑ SG↑" = "#456187" , 
                                "T↓ SG↑" =  "#F39237", #"#E89005", 
                                "T↑ SG↓" = "#52A4BC" ,
                                "T↓ SG↓" = "#A1C084"),
                                #"T↑ SG↑" = "#E3D081" , 
                                #"T↓ SG↑" =  "#8ABFCD", #"#E89005", 
                                #"T↑ SG↓" = "#52A4BC" ,
                                #"T↓ SG↓" = "#C9546C"),
                     breaks = c("T↑ SG↑", 
                                "T↓ SG↑",
                                "T↑ SG↓",
                                "T↓ SG↓")) +
  labs(y = "log10 sgRNA copies / sample", 
       x = "log10 total RNA copies / sample",
       color = "Target Gene") +
  guides(color=guide_legend(keywidth=0.35,
                           override.aes = list(size = 2.7,
                                               linewidth = 2, alpha = 1))) +
  theme(legend.position = c(0.2, 0.81), 
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        axis.title = element_text(size = 8),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 7),
        axis.title.x = element_blank(),
        legend.key=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        text = element_text(size=12),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(0, 0, 0, 0)) + 
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), limits = c(0, 10),
                     labels = c(0, 2, 4, 6, 8, ""), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10), limits = c(0.01, 10),
                     expand = c(0, 0)); fig.quant.tg.fits


### 3G: SPECIES ----------------------------------------------------------------

# Generate lines for different species, holding all else constant
sp <- sample(rep(1:3, each = n_lines / 3))
set.seed(11)
lines.sp.df <- as.data.frame(matrix(nrow = 0, ncol = 4))


for (ii in 1:length(sp)) {
  row <- sample(1:nrow(param.best.df), size = 1)
  
  if (sp[ii] == 1){
    sp_sample <- param.best.df$'betaSP[1]'[row]
  }
  else if (sp[ii] == 2){
    sp_sample <- param.best.df$'betaSP[2]'[row]
  }
  else if (sp[ii] == 3){
    sp_sample <- param.best.df$'betaSP[3]'[row]
  }
  
  y_vals <- param.best.df$alpha[row] +
    param.best.df$betaT[row] * x_vals +
    param.best.df$'betaTG[3]'[row] + 
    sp_sample + 
    param.best.df$'betaDPI[2]'[row] +
    param.best.df$betaDOSE[row] * 5.5
  
  line_num <- rep(ii, length(x_vals))
  sp_idx <- rep(sp[ii], length(x_vals))
  dpi_idx <- NA
  dose_val <- NA
  
  lines.sp.df <- rbind(lines.sp.df,
                        cbind(line_num, sp_idx, dpi_idx, dose_val, x_vals, y_vals))
}

lines.sp.df$sp_idx[lines.sp.df$sp_idx == 1] <- "RM"
lines.sp.df$sp_idx[lines.sp.df$sp_idx == 2] <- "CM"
lines.sp.df$sp_idx[lines.sp.df$sp_idx == 3] <- "AGM"


fig.quant.sp.fits <- ggplot(lines.sp.df) + 
  geom_line(aes(x = x_vals, y = y_vals, 
                group = line_num,
                color = sp_idx),
            size = 0.8, alpha = 0.1) +
  geom_line(data = as.data.frame(lines.simple.df),
            aes(x = x_vals, y = y_vals),
            color = "#002C3D",
            size = 1) +
  scale_color_manual(values = c("RM" = "#0B6884",
                                "CM" = "#744468",
                                "AGM" = "#D56540")) +
  labs(y = "", color = "Species") +
  guides(color = guide_legend(override.aes = list(size = 2.7, alpha = 1,
                                                  linewidth = 2), 
                              keyheight = 0.6)) + #,
  theme(legend.position = c(0.14, 0.81), 
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 7),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(0, 0, 0, 0)) + 
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), limits = c(-0.01, 10),
                     labels = c("", 2, 4, 6, 8, 10),
                     expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10), limits = c(-0.01, 10),
                     expand = c(0, 0)); fig.quant.sp.fits


### 3G: DOSE -------------------------------------------------------------------

# Generate lines for different target genes, holding all else constant
dose <- runif(n_lines, min(dat$log10_dose_pfu), max(dat$log10_dose_pfu))
set.seed(111)
lines.dose.df <- as.data.frame(matrix(nrow = 0, ncol = 4))

for (ii in 1:length(dose)) {
  row <- sample(1:nrow(param.best.df), size = 1)
  
  y_vals <- param.best.df$alpha[row] +
    param.best.df$betaT[row] * x_vals +
    param.best.df$'betaTG[3]'[row] + 
    param.best.df$'betaSP[1]'[row] +
    param.best.df$'betaDPI[2]'[row] +
    param.best.df$betaDOSE[row] * dose[ii]
  
  line_num <- rep(ii, length(x_vals))
  dpi_idx <- NA
  sp_idx <- NA
  dose_val <- dose[ii]
  
  lines.dose.df <- rbind(lines.dose.df,
                         cbind(line_num, sp_idx, dpi_idx, dose_val, x_vals, y_vals))
  
}

fig.quant.dose.fits <- ggplot(data = lines.dose.df) + 
  geom_line(aes(x = x_vals, y = y_vals, 
                group = line_num,
                color = as.numeric(dose_val)),
            size = 0.8, alpha = 0.3) +
  geom_line(data = as.data.frame(lines.simple.df),
            aes(x = x_vals, y = y_vals),
            color = "#002C3D",
            size = 1) +
  scale_color_gradientn(colours = zissou_pal) +
  labs(y = "", x = "log10 total RNA copies / sample",
       color = "log10\ndose") +
  guides(color = guide_colorbar(frame.colour = "black",
                                title.vjust = 0)) +
  theme(legend.position = c(0.12, 0.72), 
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 7),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(0, 0, 0, 0)) + 
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), limits = c(-0.01, 10),
                     labels = c(0, 2, 4, 6, 8, 10),
                     expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10), limits = c(-0.01, 10),
                     labels = c(0, 2, 4, 6, 8, ""),
                     expand = c(0, 0)); fig.quant.dose.fits


## Calculate tRNA copies neeeded for sgRNA to be 4
#x_small <- (4 - (param.best.df.wide$alpha +
#                   param.best.df.wide$'betaTG[1]' +
#                   param.best.df.wide$'betaSP[1]' +
#                   param.best.df.wide$'betaDPI[2]' +
#                   param.best.df.wide$betaDOSE * 4))/param.best.df.wide$betaT
#
#
#median(x_small)
#quantile(x_small, c(0.05, 0.95))


#x_big <- (4 - (param.best.df.wide$alpha +
#                 param.best.df.wide$'betaTG[1]' +
#                 param.best.df.wide$'betaSP[1]' +
#                 param.best.df.wide$'betaDPI[2]' +
#                 param.best.df.wide$betaDOSE * 7))/param.best.df.wide$betaT


#median(x_big)
#quantile(x_big, c(0.05, 0.95))

### 3G: DPI --------------------------------------------------------------------

# Generate lines for different target genes, holding all else constant
dpi <- sample(rep(1:3, each = n_lines / 3))
set.seed(1111)
lines.dpi.df <- as.data.frame(matrix(nrow = 0, ncol = 4))

for (ii in 1:length(dpi)) {
  row <- sample(1:nrow(param.best.df), size = 1)
  
  if (dpi[ii] == 1){
    dpi_sample <- param.best.df$'betaDPI[1]'[row]
  }
  else if (dpi[ii] == 2){
    dpi_sample <- param.best.df$'betaDPI[2]'[row]
  }
  else if (dpi[ii] == 3){
    dpi_sample <- param.best.df$'betaDPI[3]'[row]
  }
  
  y_vals <- param.best.df$alpha[row] +
    param.best.df$betaT[row] * x_vals +
    param.best.df$'betaTG[3]'[row] + 
    param.best.df$'betaSP[1]'[row] +
    dpi_sample +
    param.best.df$betaDOSE[row] * 5.5
  
  line_num <- rep(ii, length(x_vals))
  dpi_idx <- rep(dpi[ii], length(x_vals))
  sp_idx <- NA
  dose_val <- NA
  
  lines.dpi.df <- rbind(lines.dpi.df,
                       cbind(line_num, sp_idx, dpi_idx, dose_val, x_vals, y_vals))
}

lines.dpi.df$dpi_idx[lines.dpi.df$dpi_idx == 1] <- "I, 1"
lines.dpi.df$dpi_idx[lines.dpi.df$dpi_idx == 2] <- "I, 2+"
lines.dpi.df$dpi_idx[lines.dpi.df$dpi_idx == 3] <- "NI, 1+"


fig.quant.dpi.fits <- ggplot(lines.dpi.df) + 
  geom_line(aes(x = x_vals, y = y_vals, 
                group = line_num,
                color = dpi_idx),
            size = 0.8, alpha = 0.1) +
  geom_line(data = as.data.frame(lines.simple.df),
            aes(x = x_vals, y = y_vals),
            color = "#002C3D",
            size = 1) +
  scale_color_manual(values = c("I, 1" = "#37123C",
                                "I, 2+" = "#A83869",
                                "NI, 1+" = "#D986B0")) +
  labs(y = "", color = "DPI") +
  guides(color = guide_legend(override.aes = list(size = 2.7, alpha = 1,
                                                  linewidth = 2), 
                              keyheight = 0.6)) + #,
  theme(legend.position = c(0.15, 0.81), 
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(0, 0, 0, 0)) + 
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), limits = c(-0.01, 10),
                     expand = c(0, 0),
                     labels = c("", 2, 4, 6, 8, 10)) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10), limits = c(-0.01, 10),
                     labels = c(0, 2, 4, 6, 8, 10),
                     expand = c(0, 0)); fig.quant.dpi.fits

# Calculate tRNA copies for certain sgRNA quantity
#x_dpi1 <- (4 - (param.best.df.wide$alpha +
#                  param.best.df.wide$'betaTG[1]' +
#                  param.best.df.wide$'betaSP[1]' +
#                  param.best.df.wide$'betaDPI[1]' +
#                  param.best.df.wide$betaDOSE * 5.5))/param.best.df.wide$betaT
#x_dpi2 <- (4 - (param.best.df.wide$alpha +
#                  param.best.df.wide$'betaTG[1]' +
#                  param.best.df.wide$'betaSP[1]' +
#                  param.best.df.wide$'betaDPI[2]' +
#                  param.best.df.wide$betaDOSE * 5.5))/param.best.df.wide$betaT
#x_dpi3 <- (4 - (param.best.df.wide$alpha +
#                  param.best.df.wide$'betaTG[1]' +
#                  param.best.df.wide$'betaSP[1]' +
#                  param.best.df.wide$'betaDPI[3]' +
#                  param.best.df.wide$betaDOSE * 5.5))/param.best.df.wide$betaT
#
#
#median(x_dpi1); median(x_dpi2); median(x_dpi3)
#quantile(x_dpi1, c(0.05, 0.95)); quantile(x_dpi2, c(0.05, 0.95)); quantile(x_dpi3, c(0.05, 0.95))



### Combine into one plot ------------------------------------------------------

fig3G <- fig.quant.dose.fits + labs(tag = "G") + fig.quant.dpi.fits + fig.quant.tg.fits +
  fig.quant.sp.fits + plot_layout(ncol = 2); fig3G

x_title <- ggplot() + 
  annotate(geom = "text", x = 1, y = 1.7, size = 2.8,
           label = "log10 total RNA copies / sample") +
  theme(plot.margin = c(-10, -10, -10, -10)) +
  scale_y_continuous(limits = c(0, 2)) +
  theme_void()


fig3G <- ggarrange(fig3G, x_title, nrow = 2, heights = c(1, 0.05)); fig3G


## Combine 3EF and 3G  ---------------------------------------------------------

fig.lin.left <- ggarrange(fig3EF, fig3G, ncol = 2, widths = c(1.5, 2))


## 3H: Heat maps ---------------------------------------------------------------

x_vals <- seq(from = 0, to = 10, by = 1)


### 3F: Simple predictions -----------------------------------------------------

pars_simple <- c("alpha", "betaT", "sigma")
param.simple.df <- fit.simple$draws(variables = pars_simple,
                                    inc_warmup = FALSE,
                                    format = "df")
param.simple.df <- as.data.frame(param.simple.df)

# Generate values for the line
beta_mean <- median(param.simple.df$betaT)
alpha_mean <- median(param.simple.df$alpha)

# Generate values & a data frame for the reference line
y_vals <- alpha_mean + beta_mean * x_vals

# Combine into df for plotting
lines.simple.df <- as.data.frame(cbind(x_vals, y_vals))
lines.simple.df$predictor <- ""


### Get best model params -------------------------------------------------------

pars_best <- c("alpha", "betaT", 
               "betaTG[1]", "betaTG[2]", "betaTG[3]", "betaTG[4]",
               "betaSP[1]", "betaSP[2]", "betaSP[3]",
               "betaDPI[1]", "betaDPI[2]", "betaDPI[3]",
               "betaDOSE",
               "sigma_bar",
               "sigma_sd")
param.best.df.wide <- fit.best$draws(variables = pars_best,
                                     inc_warmup = FALSE,
                                     format = "df")
param.best.df.wide <- as.data.frame(param.best.df.wide)


#### Determine in what range simple model has different outcome ----------------

val.sg.df <- data.frame(dose = NA, sp = NA, dpi = NA, tg = NA, 
                        x_val = NA, median_val = NA, quantile_val = NA,
                        simple_val = NA)

for (dose.option in c(4, 5.5, 7)) {
  for (dpi.option in c(1, 2, 3)){
    if (dpi.option == 1){dpi.sample <-  param.best.df.wide$'betaDPI[1]'}
    else if(dpi.option == 2){dpi.sample <- param.best.df.wide$'betaDPI[2]'}
    else if(dpi.option == 3){dpi.sample <- param.best.df.wide$'betaDPI[3]'}
    
    for (sp.option in c(1, 2, 3)){
      if (sp.option == 1){sp.sample <-  param.best.df.wide$'betaSP[1]'}
      else if(sp.option == 2){sp.sample <- param.best.df.wide$'betaSP[2]'}
      else if(sp.option == 3){sp.sample <- param.best.df.wide$'betaSP[3]'}
      
      for (tg.option in c( 2, 3, 4)){
        if (tg.option == 2){tg.sample <-  param.best.df.wide$'betaTG[2]'}
        else if(tg.option == 3){tg.sample <- param.best.df.wide$'betaTG[3]'}
        else if(tg.option == 4){tg.sample <- param.best.df.wide$'betaTG[4]'}
        
        for (x_val in x_vals) {
         y <- param.best.df.wide$alpha +
              param.best.df.wide$betaT * x_val +
              tg.sample +
              sp.sample +
              dpi.sample +
              param.best.df.wide$betaDOSE * dose.option

          y_med <- median(y) 
          y_quant <- paste0(round(quantile(y, c(0.05, 0.95)), 2), collapse = ", ")
          
          val.sg.df <- rbind(val.sg.df,
                               data.frame(
                                 dose = dose.option,
                                 sp = sp.option,
                                 dpi = dpi.option,
                                 tg = tg.option,
                                 x_val = x_val, 
                                 median_val = y_med, 
                                 quantile_val = y_quant,
                                 simple_val = lines.simple.df$y_vals[lines.simple.df$x_vals == x_val]))
        }
      }
    }
  }
}

sort(unique(val.sg.df$x_val[(val.sg.df$median_val > 1.69 & val.sg.df$simple_val <= 1.69) |
                                val.sg.df$median_val  <= 1.69 & val.sg.df$simple_val > 1.69]))

val.sg.df <- val.sg.df[-1, ]

#### Make the Supplemental Table for 90% CrI -------------------------------------

tblS7 <- subset(val.sg.df, select = -c(median_val, simple_val))

tblS7 <- tblS7 %>%
  pivot_wider(names_from = x_val, values_from = c(quantile_val))

tblS7.doses <- subset(tblS7, dose %in% c(4, 5.5, 7) & sp == 1 & tg == 3 & dpi == 2) %>%
  rename(predictor = dose) %>%
  select(-c(sp, tg, dpi)); tblS7.doses

tblS7.dpi <- subset(tblS7, dose == 5.5 & sp == 1 & tg == 3) %>%
  rename(predictor = dpi) %>%
  select(-c(sp, tg, dose)); tblS7.dpi

tblS7.species <- subset(tblS7, dose == 5.5 & tg == 3 & dpi == 2) %>%
  rename(predictor = sp) %>%
  select(-c(dose, tg, dpi)); tblS7.species

tblS7.tg <- subset(tblS7, dose == 5.5 & sp == 1 & tg %in% c(2, 3, 4) & dpi == 2) %>%
  rename(predictor = tg) %>%
  select(-c(dose, sp, dpi)); tblS7.tg

tblS7.combined <- rbind(tblS7.doses, tblS7.dpi, tblS7.species, tblS7.tg); tblS7.combined
tblS7.combined$group <- c(rep("Dose", 3),
                          rep("DPI", 3),
                          rep("Species", 3),
                          rep("Target Gene", 3))
tblS7.combined <- tblS7.combined %>%
  select(c(group, predictor, everything())); tblS7.combined


write.csv(tblS7.combined, file = "./outputs/tables/tblS6B-sgRNA-linear-90CrI.csv",
          row.names = FALSE)




### 3F: SIMPLE MODEL  ----------------------------------------------------------

brewer.pal <- rev(brewer.pal(10, "RdBu"))[-c(1, 10)]
brewer.pal[4] <- "white"

fig.quant.simple <- ggplot(as.data.frame(lines.simple.df), 
                     aes(x = as.numeric(x_vals), 
                         y = predictor, 
                         fill = as.numeric(y_vals))) +
  geom_tile() + 
  geom_rect(aes(xmin = 0.5, xmax = 5.5, ymin = 0.5, ymax = 1.5), 
            fill = "transparent", color = "grey33", size = 1.5) +
  geom_text(aes(label = (round(as.numeric(y_vals), 1))), color = "black", size = 2.5) +
  scale_fill_gradientn(colours = brewer.pal,
                       limits = c(-2.2, 10),
                       values = rescale(c(seq(-2.2, 1.7, length.out = 4),
                                          seq(1.7 + (10-1.7)/4, 10, length.out = 4))),
                       breaks = c(0, 3, 6, 9)) +
  labs(x = "log10 total RNA copies / sample",
       y = "Simple") +
  theme(legend.position = "none",
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 7),
        axis.text.y = element_text(angle = 0, 
                                   vjust = 1, hjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(0, 0, 0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12), expand = c(0, 0)); fig.quant.simple





### 3F: TARGET GENE ------------------------------------------------------------

tg.df <- data.frame(tRNA = NA,
                    tg = NA,
                    probs = NA)

for (x_val in x_vals) {
  #y <- param.best.df.wide$alpha +
  #  param.best.df.wide$betaT * x_val +
  #  param.best.df.wide$'betaTG[1]' +
  #  param.best.df.wide$'betaSP[1]' +
  #  param.best.df.wide$'betaDPI[2]' +
  #  param.best.df.wide$betaDOSE * 5.5

  #
  #y_tg1 <- median(y)
  
  y <- param.best.df.wide$alpha +
    param.best.df.wide$betaT * x_val +
    param.best.df.wide$'betaTG[2]' +
    param.best.df.wide$'betaSP[1]' +
    param.best.df.wide$'betaDPI[2]' +
    param.best.df.wide$betaDOSE * 5.5
  
  y_tg2 <- median(y)
  
  y <- param.best.df.wide$alpha +
    param.best.df.wide$betaT * x_val +
    param.best.df.wide$'betaTG[3]' +
    param.best.df.wide$'betaSP[1]' +
    param.best.df.wide$'betaDPI[2]' +
    param.best.df.wide$betaDOSE * 5.5
  
  y_tg3<- median(y)
  
  y <- param.best.df.wide$alpha +
    param.best.df.wide$betaT * x_val +
    param.best.df.wide$'betaTG[4]' +
    param.best.df.wide$'betaSP[1]' +
    param.best.df.wide$'betaDPI[2]' +
    param.best.df.wide$betaDOSE * 5.5
  
  y_tg4 <- median(y)
  
  tg.df <- rbind(tg.df, 
                 #data.frame(tRNA = x_val, 
                 #           tg = "T↑ SG↑", 
                 #           probs = y_tg1),
                 data.frame(tRNA = x_val, 
                            tg = "T↓ SG↑",
                            probs = y_tg2),
                 data.frame(tRNA = x_val, 
                            tg = "T↑ SG↓",
                            probs = y_tg3),
                 data.frame(tRNA = x_val, 
                            tg = "T↓ SG↓",
                            probs = y_tg4))

}

tg.df <- tg.df[-1, ]

fig.quant.tg <- ggplot(tg.df, aes(x = as.numeric(tRNA), 
                                  y = factor(tg, levels = rev(c("T↑ SG↑", 
                                                                "T↓ SG↑",
                                                                "T↑ SG↓",
                                                                "T↓ SG↓"))),
                                  fill = as.numeric(probs))) +
  geom_tile() + 
  geom_text(aes(label = (round(as.numeric(probs), 1))), color = "black", size = 2.5) +
  geom_rect(aes(xmin = 1.5, xmax = 3.5, ymin = 0.5, ymax = 3.5), 
            fill = "transparent", color = "grey33", size = 1.5) +
  scale_fill_gradientn(colours = brewer.pal,
                       limits = c(-2.2, 10),
                       values = rescale(c(seq(-2.2, 1.7, length.out = 4),
                                          seq(1.7 + (10-1.7)/4, 10, length.out = 4))),
                       breaks = c(0, 3, 6, 9)) +
  labs(x = "log10 total RNA copies / sample", 
       y = "Target Gene") + 
  theme(legend.position = "none",
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 7),
        axis.text.y = element_text(angle = 45, size = 6.5, face = c("plain","bold","plain"), 
                                   vjust = 1, hjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(0, 0, 0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12), expand = c(0, 0)); fig.quant.tg


### 3F: SPECIES ----------------------------------------------------------------

sp.df <- data.frame(tRNA = NA,
                    sp = NA,
                    probs = NA)

for (x_val in x_vals) {
  y <- param.best.df.wide$alpha +
    param.best.df.wide$betaT * x_val +
    param.best.df.wide$'betaTG[3]' +
    param.best.df.wide$'betaSP[1]' +
    param.best.df.wide$'betaDPI[2]' +
    param.best.df.wide$betaDOSE * 5.5
  
  y_rm <- median(y)
  
  y <- param.best.df.wide$alpha +
    param.best.df.wide$betaT * x_val +
    param.best.df.wide$'betaTG[3]' +
    param.best.df.wide$'betaSP[2]' +
    param.best.df.wide$'betaDPI[2]' +
    param.best.df.wide$betaDOSE * 5.5
  
  y_cm <- median(y)
  
  y <- param.best.df.wide$alpha +
    param.best.df.wide$betaT * x_val +
    param.best.df.wide$'betaTG[3]' +
    param.best.df.wide$'betaSP[3]' +
    param.best.df.wide$'betaDPI[2]' +
    param.best.df.wide$betaDOSE * 5.5
  
  y_agm <- median(y)
  
  sp.df <- rbind(sp.df, 
                 data.frame(tRNA = x_val, 
                            sp = "RM", 
                            probs = y_rm),
                 data.frame(tRNA = x_val, 
                            sp = "CM",
                            probs = y_cm),
                 data.frame(tRNA = x_val, 
                            sp = "AGM",
                            probs = y_agm))
}

sp.df <- sp.df[-1, ]

fig.quant.sp <- ggplot(sp.df, aes(x = as.numeric(tRNA), 
                                  y = factor(sp, levels = c("AGM", "CM", "RM")), 
                                  fill = as.numeric(probs))) +
  geom_tile() + 
  geom_text(aes(label = (round(as.numeric(probs), 1))), color = "black", size = 2.5) +
  geom_rect(aes(xmin = 2.5, xmax = 4.5, ymin = 0.5, ymax = 3.5), 
            fill = "transparent", color = "grey33", size = 1.5) +
  scale_fill_gradientn(colours = brewer.pal,
                       limits = c(-2.2, 10),
                       values = rescale(c(seq(-2.2, 1.7, length.out = 4),
                                          seq(1.7 + (10-1.7)/4, 10, length.out = 4))),
                       breaks = c(0, 3, 6, 9)) +
  guides(fill = guide_colorbar(frame.colour = "black")) +
  labs(x = "log10 total RNA copies / sample", 
       y = "Species",
       fill = "log10 sgRNA copy numbers") + 
  theme(legend.position = "bottom",
        legend.title = element_text(size = 7.5, vjust = 0.7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.margin = margin(t = -1),
        axis.title = element_text(size = 7),
        axis.text.y = element_text(angle = 45, size = 6.5, 
                                   face = c("plain","plain", "bold"), 
                                   vjust = 1, hjust = 0.5),
        axis.text.x = element_text(size = 7),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(0, 0, 0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(-0.5, 10, 2), labels = seq(0, 10, 2),
                     limits = c(-0.5, 10.5), expand = c(0, 0)); fig.quant.sp


### 3F: DOSE -------------------------------------------------------------------

dose.df <- data.frame(tRNA = NA,
                      dose = NA,
                      probs = NA)

for (x_val in x_vals) {
  y <- param.best.df.wide$alpha +
       param.best.df.wide$betaT * x_val +
       param.best.df.wide$'betaTG[3]' +
       param.best.df.wide$'betaSP[1]' +
       param.best.df.wide$'betaDPI[2]' +
       param.best.df.wide$betaDOSE * 4
  
  y_small <- median(y)
  
  y <- param.best.df.wide$alpha +
    param.best.df.wide$betaT * x_val +
    param.best.df.wide$'betaTG[3]' +
    param.best.df.wide$'betaSP[1]' +
    param.best.df.wide$'betaDPI[2]' +
    param.best.df.wide$betaDOSE * 5.5
  
  y_mid <- median(y)
  
  y <- param.best.df.wide$alpha +
    param.best.df.wide$betaT * x_val +
    param.best.df.wide$'betaTG[3]' +
    param.best.df.wide$'betaSP[1]' +
    param.best.df.wide$'betaDPI[2]' +
    param.best.df.wide$betaDOSE * 7
  
  y_big <- median(y)
  
  dose.df <- rbind(dose.df, 
                   data.frame(tRNA = x_val, 
                              dose = "10^4", 
                              probs = y_small),
                   data.frame(tRNA = x_val, 
                              dose = "10^5.5",
                              probs = y_mid),
                   data.frame(tRNA = x_val, 
                              dose = "10^7",
                              probs = y_big))
}

dose.df <- dose.df[-1, ]





fig.quant.dose <- ggplot(dose.df, aes(x = as.numeric(tRNA), 
                                      y = factor(dose, levels = c("10^7", "10^5.5", "10^4")), 
                                      fill = as.numeric(probs))) +
  geom_tile() + 
  geom_text(aes(label = (round(as.numeric(probs), 1))), color = "black", size = 2.5) +
  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = 0.5, ymax = 3.5), 
            fill = "transparent", color = "grey33", size = 1.5) +
  scale_fill_gradientn(colours = brewer.pal,
                       limits = c(-2.2, 10),
                       values = rescale(c(seq(-2.2, 1.7, length.out = 4),
                                          seq(1.7 + (10-1.7)/4, 10, length.out = 4))),
                       breaks = c(0, 3, 6, 9)) +
  labs(x = "log10 total RNA copies / sample", 
       y = "log10 Dose") + 
  theme(legend.position = "none",
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 7),
        axis.text.y = element_text(angle = 45, size = 6.5, face = c("plain", "bold", "plain"), 
                                   vjust = 1, hjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(0, 0, 0, 0)) + 
  scale_y_discrete(expand = c(0, 0), labels = c("7", "5.5", "4")) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12), expand = c(0, 0)); fig.quant.dose


### Predictions of totRNA quantities per sgRNA quant

x <- (4 - (param.best.df.wide$alpha +
          param.best.df.wide$'betaTG[3]' +
          param.best.df.wide$'betaSP[1]' +
          param.best.df.wide$'betaDPI[2]' +
          param.best.df.wide$betaDOSE * 4))/ param.best.df.wide$betaT

median(x); quantile(x, c(0.05, 0.95))

x <- (4 - (param.best.df.wide$alpha +
             param.best.df.wide$'betaTG[3]' +
             param.best.df.wide$'betaSP[1]' +
             param.best.df.wide$'betaDPI[2]' +
             param.best.df.wide$betaDOSE * 7))/ param.best.df.wide$betaT

median(x); quantile(x, c(0.05, 0.95))


### 3F: DPI --------------------------------------------------------------------

dpi.df <- data.frame(tRNA = NA,
                     dpi = NA,
                     probs = NA)

for (x_val in x_vals) {
  y <- param.best.df.wide$alpha +
    param.best.df.wide$betaT * x_val +
    param.best.df.wide$'betaTG[3]' +
    param.best.df.wide$'betaSP[1]' +
    param.best.df.wide$'betaDPI[1]' +
    param.best.df.wide$betaDOSE * 5.5
  
  y_dpi1 <- median(y)
  
  y <- param.best.df.wide$alpha +
    param.best.df.wide$betaT * x_val +
    param.best.df.wide$'betaTG[3]' +
    param.best.df.wide$'betaSP[1]' +
    param.best.df.wide$'betaDPI[2]' +
    param.best.df.wide$betaDOSE * 5.5
  
  y_dpi2 <- median(y)
  
  y <- param.best.df.wide$alpha +
    param.best.df.wide$betaT * x_val +
    param.best.df.wide$'betaTG[3]' +
    param.best.df.wide$'betaSP[1]' +
    param.best.df.wide$'betaDPI[3]' +
    param.best.df.wide$betaDOSE * 5.5
  
  y_dpi3 <- median(y)
  
  dpi.df <- rbind(dpi.df, 
                 data.frame(tRNA = x_val, 
                            dpi = "I, 1", 
                            probs = y_dpi1),
                 data.frame(tRNA = x_val, 
                            dpi = "I, 2+",
                            probs = y_dpi2),
                 data.frame(tRNA = x_val, 
                            dpi = "NI, 1+",
                            probs = y_dpi3))
}

dpi.df <- dpi.df[-1, ]

fig.quant.dpi <- ggplot(dpi.df, aes(x = as.numeric(tRNA), 
                                   y = factor(dpi, levels = c("NI, 1+", "I, 2+", "I, 1")), 
                                   fill = as.numeric(probs))) +
  geom_tile() + 
  geom_text(aes(label = (round(as.numeric(probs), 1))), color = "black", size = 2.5) +
  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = 0.5, ymax = 3.5), 
            fill = "transparent", color = "grey33", size = 1.5) +
  scale_fill_gradientn(colours = brewer.pal,
                       limits = c(-2.2, 10),
                       values = rescale(c(seq(-2.2, 1.7, length.out = 4),
                                          seq(1.7 + (10-1.7)/4, 10, length.out = 4))),
                       breaks = c(0, 3, 6, 9)) +
  labs(x = "log10 total RNA copies / sample", 
       y = "DPI",
       fill = "log10 sgRNA copy numbers") + 
  guides(fill = guide_colorbar(frame.colour = "black")) +
  theme(legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 7),
        axis.text.y = element_text(angle = 45, size = 6.5, face = c("plain","bold","plain"), 
                                   vjust = 1, hjust = 0.5),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(0, 0, 0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(-0.5, 10.5, 2), 
                     labels = seq(0, 10, 2),
                     expand = c(0, 0)); fig.quant.dpi



### Predictions of totRNA quantities per sgRNA quant

x <- (4 - (param.best.df.wide$alpha +
             param.best.df.wide$'betaTG[3]' +
             param.best.df.wide$'betaSP[1]' +
             param.best.df.wide$'betaDPI[1]' +
             param.best.df.wide$betaDOSE * 5.5))/ param.best.df.wide$betaT

median(x); quantile(x, c(0.05, 0.95))

x <- (4 - (param.best.df.wide$alpha +
             param.best.df.wide$'betaTG[3]' +
             param.best.df.wide$'betaSP[1]' +
             param.best.df.wide$'betaDPI[2]' +
             param.best.df.wide$betaDOSE * 5.5))/ param.best.df.wide$betaT

median(x); quantile(x, c(0.05, 0.95))

x <- (4 - (param.best.df.wide$alpha +
             param.best.df.wide$'betaTG[3]' +
             param.best.df.wide$'betaSP[1]' +
             param.best.df.wide$'betaDPI[3]' +
             param.best.df.wide$betaDOSE * 5.5))/ param.best.df.wide$betaT

median(x); quantile(x, c(0.05, 0.95))


### 3H: Combine -----------------------------------------------------------------

fig.lin.right <- fig.quant.simple + labs(tag = "H") + fig.quant.dose + 
  fig.quant.dpi + fig.quant.tg + fig.quant.sp +
  plot_layout(ncol = 1, heights = c(1, 3, 3, 3, 3)); fig.lin.right



## 3E-H: All Linear Panels -----------------------------------------------------

fig.lin <- ggarrange(fig.lin.left, fig.lin.right, ncol = 2, widths = c(1, 0.4))

# Generate text-only figure for joint labels
linear_label <- ggplot() + 
  annotate(geom = "text", x = 1, y = 1, size = 3.2, 
           label = "Linear Component", angle = -90,
           fontface = 2) +
  geom_segment(aes(x = 0.6, xend = 0.6, y = -0.5, yend = 2.5), color = "grey42") +
  theme(plot.margin = c(-10, -10, -10, -10)) +
  scale_x_continuous(limits = c(0, 2)) +
  theme_void(); linear_label


fig.lin <- ggarrange(fig.lin, linear_label, ncol = 2, widths = c(1, 0.04))


# Combine into one figure ------------------------------------------------------

fig3 <- ggarrange(fig.log, fig.lin, nrow = 2, heights = c(1, 1.05))


# Save -------------------------------------------------------------------------

ggsave('./outputs/figures/fig3-sgRNA-model-results.tiff',
       plot = fig3,
       device = 'tiff',
       height = 8.5,
       width = 11,
       units = 'in',
       bg = 'white')

