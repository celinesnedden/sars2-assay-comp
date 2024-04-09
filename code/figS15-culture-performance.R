# Generates figure SX from the manuscript, including raw data & model fits,
#     and it generates many of the stats/predictions described in the text
# No dependencies, this file will run without running any other scripts


# Prep environment -------------------------------------------------------------

req_pkgs <- c("loo", "stringr", "tidyverse", "ggridges", "bayesplot", "patchwork",
              "gtools", "wesanderson", "cowplot", "ggpubr", "matrixStats", "pROC",
              "RColorBrewer")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependences = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))


# Small convenience function
`%notin%` <- Negate(`%in%`)  # for convenience

# Extract desired colored palette
zissou_pal <- wes_palette("Zissou1", 12, type = "continuous")
brewer_pal <- brewer.pal(11, "RdBu")[2:10]

# Load data --------------------------------------------------------------------

dat <- read.csv("./data/pred-sg-data.csv")

# Convert columns to numerics
dat$val_total <- as.numeric(dat$val_total) # total RNA pos censored data will trigger NA conversion warning


# Prep data --------------------------------------------------------------------

dat.culture <- subset(dat, !is.na(pos_inf) & pos_total %in% c(0, 1))

# Convert columns to numerics
dat.culture$val_total <- as.numeric(dat.culture$val_total)

# Set TG indicator based on total RNA protocol: ordered by decreasing RNA quantity
dat.culture$tg_idx[dat.culture$tg_total %in% c("N")] <- 1
dat.culture$tg_idx[dat.culture$tg_total %in% c("E")] <- 2
dat.culture$tg_idx[dat.culture$tg_total %in% c("S")] <- 3


# Load model fits --------------------------------------------------------------

fit.best <- readRDS(file = "./outputs/fits/fit-culture-best-inf-priors.RDS")
fit.simple <- readRDS(file = "./outputs/fits/fit-culture-simple-inf-priors.RDS")


# Extract parameters -----------------------------------------------------------

## Best model ------------------------------------------------------------------

# Extract parameter estimates from best model
pars_best <- c("gamma", "psiT", 
               "psiDPI[1]",  "psiDPI[2]", "psiDPI[3]",
               "psiAGE[1]", "psiAGE[2]", "psiAGE[3]",
               "psiSP[1]", "psiSP[2]", "psiSP[3]",
               "psiASSAY",
               "psiTG[1]", "psiTG[2]", "psiTG[3]",
               "psiDOSE", 
               "psiCELL[1]", "psiCELL[2]", "psiCELL[3]")

param.best.df.wide <- fit.best$draws(variables = pars_best,
                                     inc_warmup = FALSE,
                                     format = "df")
param.best.df <- param.best.df.wide %>% pivot_longer(cols = everything(),
                                                     names_to = "param",
                                                     values_to = "value")
param.best.df <- as.data.frame(param.best.df)
param.best.df <- subset(param.best.df, param %notin% c(".chain", ".iteration", ".draw"))
param.best.df$param <- factor(param.best.df$param, levels = rev(pars_best))


## Simple model ----------------------------------------------------------------

# Extract relevant parameter estimates
pars_simple <- c("gamma", "psiT")
param.simple.df.wide <- fit.simple$draws(variables = pars_simple,
                                         inc_warmup = FALSE,
                                         format = "df")
param.simple.df <- param.simple.df.wide %>% pivot_longer(cols = everything(),
                                                         names_to = "param",
                                                         values_to = "value")
param.simple.df <- as.data.frame(param.simple.df)
param.simple.df <- subset(param.simple.df, param %notin% c(".chain", ".iteration", ".draw"))
param.simple.df$param <- factor(param.simple.df$param, levels = rev(pars_simple))

# Extract the median reference line
psi_median <- median(param.simple.df.wide$psiT)
gamma_median <- median(param.simple.df.wide$gamma)

# Generate values & a data frame for the reference line
x_vals <- seq(from = 0.5, to = 12.5, by = 1)
y_vals <- exp(gamma_median + psi_median * x_vals)/(1 + exp(gamma_median + psi_median * x_vals))
lines.simple.df <- as.data.frame(cbind(x_vals, y_vals))


# S15A: Intermediate totRNA quantities ------------------------------------------

# Set threshold value for classifying as > or < the LOD
threshold <- 0.5

# Subset to total RNA positive data
dat.culture <- subset(dat, !is.na(pos_inf) & cens_total == 0)

### Prep simple model results --------------------------------------------------

# Extract mean model predictions of classification
prob_pos <- paste0("p_pos[", 1:nrow(dat.culture), "]") 
prob_df <- fit.simple$draws(prob_pos, format = "df")
prob_mean <- as.data.frame(colMeans(prob_df))
prob_mean <- prob_mean[rownames(prob_mean) %notin% c(".chain", ".iteration", ".draw"), ]

# Compare predictions with true observed classification
prediction <- prob_mean
prediction[prob_mean >= threshold] <- 1
prediction[prob_mean < threshold] <- 0
n_correct <- sum(dat.culture$pos_inf == prediction)

# Calculate number of samples correctly predicted 
n_correct

# Probabilities of true positives
pos_pred_simple <- prob_mean[dat.culture$pos_inf == 1 & dat.culture$val_total <= 8 &
                        dat.culture$val_total >= 6]

# Probabilities of true negatives
neg_pred_simple <- prob_mean[dat.culture$pos_inf == 0 & dat.culture$val_total <= 8 &
                        dat.culture$val_total >= 6]

median(pos_pred_simple)
median(neg_pred_simple)

# Percentages
pos_pred_correct_percent_simple <- length(pos_pred_simple[pos_pred >= 0.5]) / length(pos_pred_simple) * 100
neg_pred_correct_percent_simple <- length(neg_pred_simple[neg_pred < 0.5]) / length(neg_pred_simple) * 100


# Place in a data frame
prediction.df <- as.data.frame(matrix(data = NA, 
                                      nrow = length(c(pos_pred_simple, neg_pred_simple)),
                                      ncol = 3))
prediction.df[, 1] <- "Simple"
prediction.df[, 2] <- c(pos_pred_simple, neg_pred_simple)
prediction.df[, 3] <- c(rep("Positive", length(pos_pred_simple)),
                        rep("Negative", length(neg_pred_simple)))



### Prep full model results -----------------------------------------------------

# Extract mean model predictions of classification
prob_pos <- paste0("p_pos[", 1:nrow(dat.culture), "]") 
prob_df <- fit.best$draws(prob_pos, format = "df")
prob_mean <- as.data.frame(colMeans(prob_df))
prob_mean <- prob_mean[rownames(prob_mean) %notin% c(".chain", ".iteration", ".draw"), ]

# Compare predictions with true observed classification
prediction <- prob_mean
prediction[prob_mean >= threshold] <- 1
prediction[prob_mean < threshold] <- 0
n_correct <- sum(dat.culture$pos_inf == prediction)

# Probabilities of true positives
pos_pred_best <- prob_mean[dat.culture$pos_inf == 1  & dat.culture$val_total <= 8 &
                        dat.culture$val_total >= 6]

# Probabilities of true negatives
neg_pred_best <- prob_mean[dat.culture$pos_inf == 0  & dat.culture$val_total <= 8 &
                        dat.culture$val_total >= 6]

# Percentages
pos_pred_correct_percent_full <- length(pos_pred_best[pos_pred_best >= 0.5]) / length(pos_pred_best) * 100
neg_pred_correct_percent_full <- length(neg_pred_best[neg_pred_best < 0.5]) / length(neg_pred_best) * 100

# Median probabilities
median(pos_pred_best)
median(neg_pred_best)


# Place in joint data frame
prediction.df <- rbind(prediction.df, 
                       cbind(rep("Best", length(c(pos_pred_best, neg_pred_best))),
                             c(pos_pred_best, neg_pred_best),
                             c(rep("Positive", length(pos_pred_best)), 
                               rep("Negative", length(neg_pred_best)))))
colnames(prediction.df) <- c("model", "prob_mean", "posneg")


### Make the plots -------------------------------------------------------------

# Plot all observed culture negatives
fig3B.1 <- ggplot(subset(prediction.df, posneg == "Negative")) +
  geom_vline(xintercept = threshold, size = 0.7, color = zissou_pal[1], alpha = 0.8,
             linetype = "dashed") +
  geom_density(aes(x = (1-as.numeric(prob_mean)),
                   fill = factor(model, levels = c("Simple", "Best")),
                   color = factor(model, levels = c("Simple", "Best"))),
               alpha = 0.4) +
  annotate("text", x = 0.9, y = 4.5, label = paste0(round(neg_pred_correct_percent_full), "%"),
           color = "#85B068", size = 3.5) +
  annotate("text", x = 0.9, y = 5, label = paste0(round(neg_pred_correct_percent_simple), "%"),
           color = "#0F4C5E", size = 3.5) +
  facet_grid(.~ "Observed Culture Negative") +
  scale_fill_manual(values = c("Best" = "#85B068",
                               "Simple" = "#0F4C5E")) +
  scale_color_manual(values = c("Best" = "#56783F",
                                "Simple" = "#0F4C5E")) +
  guides(color = "none",
         fill = guide_legend(override.aes = list(alpha = 0.6,
                                                 color = c( "#0F4C5E", "#56783F")))) +
  labs(x = "Chance of Culture Neg. (%)", y = "Density",
       fill = "Model", tag = "A") + 
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0),
                     breaks = seq(0, 1, 0.25), 
                     labels = c("0","25", "50", "75", "100")) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 2),
                     expand = c(0, 0),
                     position = "left") +
  theme(legend.position = c(0.225, 0.8), 
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        legend.background = element_rect(fill = 'white'),
        axis.title = element_text(size = 8.5),
        axis.text = element_text(size = 7),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(2, 6, 2, 2)); fig3B.1


# Plot all observed culture positives
fig3B.2 <- ggplot(subset(prediction.df, posneg == "Positive")) +
  geom_vline(xintercept = threshold, size = 0.7, color = zissou_pal[1], alpha = 0.8,
             linetype = "dashed") +
  geom_density(aes(x = (as.numeric(prob_mean)),
                   fill = factor(model, levels = c("Simple", "Best")),
                   color = factor(model, levels = c("Simple", "Best"))),
               alpha = 0.4) +
  annotate("text", x = 0.9, y = 4.5, label = paste0(round(pos_pred_correct_percent_full), "%"),
           color = "#85B068", size = 3.5) +
  annotate("text", x = 0.9, y = 5, label = paste0(round(pos_pred_correct_percent_simple), "%"),
           color = "#0F4C5E", size = 3.5) +
  facet_grid(.~ "Observed Culture Positive") +
  scale_fill_manual(values = c("Best" = "#85B068",
                               "Simple" = "#0F4C5E")) +
  scale_color_manual(values = c("Best" = "#56783F",
                                "Simple" = "#0F4C5E")) +
  guides(color = "none") +
  labs(x = "Chance of Culture Pos. (%)", y = "Density",
       fill = "Model") + 
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0),
                     breaks = seq(0, 1, 0.25), 
                     labels = c("", "25", "50", "75", "100")) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 2),
                     expand = c(0, 0),
                     position = "left") +
  theme(legend.position = "none", 
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        legend.background = element_rect(fill = 'white'),
        text = element_text(size = 11),
        axis.title = element_text(size = 8.5),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(2, 2, 2, 0)); fig3B.2

fig.int <- fig3B.1 + labs(tag = "B") + fig3B.2 + plot_layout(ncol = 2); fig.int



# S15B: Per-sample intermediate quantities ----------------------------------------

pos_diffs <- pos_pred_best - pos_pred_simple
neg_diffs <- neg_pred_simple - neg_pred_best

df.prob <- data.frame(culture = c(rep("Culture Pos", length(pos_diffs)),
                                  rep("Culture Neg", length(neg_diffs))),
                      diffs = c(pos_diffs, neg_diffs))


# Make the plot
fig.int.dist <- ggplot(df.prob) +
  geom_vline(xintercept = 0, linetype = "dashed", color = zissou_pal[1]) +
  geom_density(aes(x =  diffs, color = culture, fill = culture),
               alpha = 0.4) +
  annotate("text", label = "Best Does\nBetter", x = 0.3, y = 3, size = 3) +
  annotate("text", label = "Simple Does\nBetter", x = -0.27, y = 3, size = 3) +
  scale_fill_manual(values = c("Culture Pos" = "#E4BC11",
                               "Culture Neg" = "grey62")) +
  scale_color_manual(values = c("Culture Pos" = "#E4BC11",
                                "Culture Neg" = "grey62")) +
  labs(x = "Difference in the probability\nof correct classification",
       y = "Density", fill = "", color = "",
       tag = "C") +
  scale_x_continuous(limits = c(-0.5, 0.5)) +
  scale_y_continuous(limits = c(0, 3.4), expand = c(0, 0)) +
  facet_wrap(.~ "Intermediate totRNA quantities") +
  theme(legend.position = c(0.25, 0.6),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.margin= margin(c(-4,3,3,3)),
        legend.key = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        text = element_text(size = 11),
        axis.title = element_text(size = 8.5),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(2, 5, 2, 0)); fig.int.dist






# S15B: Compare per-sample probabilities ----------------------------------------

dat.best <- read.csv("./data/pred-culture-data.csv")
dat.best <- subset(dat.best, !is.na(pos_total) & !is.na(pos_inf))

dat.simple <- read.csv("./data/pred-culture-data-simple.csv")
dat.simple <- subset(dat.simple, !is.na(pos_total) & !is.na(pos_inf))

# ALL OF THESE NEED TO BE TRUE TO BE OK TO COMPARE
unique(dat.simple$indiv == dat.best$indiv)
unique(dat.simple$sample_rep == dat.best$sample_rep)
unique(dat.simple$dpi == dat.best$dpi)
nrow(dat.simple) == nrow(dat.best)

# Calculate per sample differences
dat.best$prob_diff <- dat.best$pos_inf_prob_pred - dat.simple$pos_inf_prob_pred
dat.best$prob_diff1m <- (1-dat.best$pos_inf_prob_pred) - (1-dat.simple$pos_inf_prob_pred)


# Make the plot
fig.dist <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = zissou_pal[1]) +
  geom_density(data = subset(dat.best, pos_inf == 1 & pos_total == 1),
               aes(x =  prob_diff, color = "Culture Pos", fill = "Culture Pos"),
               alpha = 0.4) +
  geom_density(data = subset(dat.best, pos_inf == 0 & pos_total == 1),
               aes(x =  prob_diff1m, color = "Culture Neg", fill = "Culture Neg"),
               alpha = 0.4) +
  annotate("text", label = "Best Does\nBetter", x = 0.3, y = 6, size = 3) +
  annotate("text", label = "Simple Does\nBetter", x = -0.27, y = 6, size = 3) +
  scale_fill_manual(values = c("Culture Pos" = "#E4BC11",
                               "Culture Neg" = "grey62")) +
  scale_color_manual(values = c("Culture Pos" = "#E4BC11",
                                "Culture Neg" = "grey62")) +
  labs(x = "Difference in the probability\nof correct classification",
       y = "Density", fill = "", color = "",
       tag = "A") +
  scale_x_continuous(limits = c(-0.5, 0.5)) +
  scale_y_continuous(limits = c(0, 6.8), expand = c(0, 0)) +
  facet_wrap(.~ "All totRNA quantities") +
  theme(legend.position = c(0.25, 0.6),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.margin= margin(c(-4,3,3,3)),
        legend.key = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        text = element_text(size = 11),
        axis.title = element_text(size = 8.5),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(2, 2, 2, 0)); fig.dist


# Save -------------------------------------------------------------------------

fig.comb <- fig.dist + fig.int + fig.int.dist  + plot_layout(ncol = 3,
                                                             widths = c(1, 2, 1)); fig.comb


ggsave('./outputs/figures/figS15-culture-performance.tiff',
       plot = fig.comb,
       device = 'tiff',
       height = 2.8,
       width = 8.5,
       units = 'in',
       bg = 'white')

