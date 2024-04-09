# Generates Figure S13, the error analysis for the sgRNA model 
#     also calculates key statistics reported in reconstruction paragraph in the main text

# Prep environment -------------------------------------------------------------

# Install & load color palette & ggplot 
req_pkgs <- c("wesanderson", "ggplot2", "cmdstanr", "tidyverse", "ggridges",
              "ggpubr", "patchwork", "ggExtra")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependences = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))


`%notin%` <- Negate(`%in%`)  # for convenience


# Load data with predictions ---------------------------------------------------

dat <- read.csv("./data/pred-culture-data.csv")
dat <- subset(dat, pos_sg %in% c(0, 1) & cens_total %in% c(0, 1))

## Combine indiv names with sample rep to get the right connecting lines
dat$indiv_sample <- paste0(dat$indiv, "_", dat$sample_rep, 
                           "_", dat$tg_idx)


# Error classifier columns -----------------------------------------------------

# Add TP/TN/FP/FN classifier
dat$pos_sg_pred_class <- NA
dat$pos_sg_pred_class[dat$pos_sg == 1 & dat$pos_sg_pred == 1] <- "True Positive"
dat$pos_sg_pred_class[dat$pos_sg == 0 & dat$pos_sg_pred == 0] <- "True Negative"
dat$pos_sg_pred_class[dat$pos_sg == 1 & dat$pos_sg_pred == 0] <- "False Negative"
dat$pos_sg_pred_class[dat$pos_sg == 0 & dat$pos_sg_pred == 1] <- "False Positive"

# Add correct/incorrect classifier 
dat$pos_sg_pred_correct <- NA
dat$pos_sg_pred_correct[dat$pos_sg == 1 & dat$pos_sg_pred == 1] <- "Correct"
dat$pos_sg_pred_correct[dat$pos_sg == 0 & dat$pos_sg_pred == 0] <- "Correct"
dat$pos_sg_pred_correct[dat$pos_sg == 1 & dat$pos_sg_pred == 0] <- "Incorrect"
dat$pos_sg_pred_correct[dat$pos_sg == 0 & dat$pos_sg_pred == 1] <- "Incorrect"


# Plot individual trajectories by errors ---------------------------------------

dat.ni <- subset(dat, sample_type == "Non-invasive")

indiv_list <- c()
indiv_correct_list <- c()
for (indiv.ii in unique(dat.ni$indiv_sample)) {
  indiv.sub <- subset(dat.ni, indiv_sample == indiv.ii)
  if (nrow(indiv.sub) >= 2) {
    if (str_detect(paste0(indiv.sub$pos_sg_pred_correct, collapse = ""), "Incorrect")){
      indiv_list <- c(indiv_list, indiv.ii)
    }
    else {
      indiv_correct_list <- c(indiv_correct_list, indiv.ii)
    }
  }
}

dat.ni$error_idx[dat.ni$indiv_sample %in% indiv_list] <- "Some errors"
dat.ni$error_idx[dat.ni$indiv_sample %in% indiv_correct_list] <- "No errors"

dat.traj <- subset(dat.ni, !is.na(error_idx))


figSXA <- ggplot(dat.traj) + 
  
  geom_line(aes(x = dpi, y = indiv_sample),
            color = "#E4BC11", 
            alpha = 0.3) +
  geom_point(aes(x = dpi, y = indiv_sample, fill = as.character(pos_sg),
                 color = pos_sg_pred_correct),
             shape = 21, size = 1.5) +
  
  # Set colors
  scale_fill_manual(values = c("0" = subset(traj_palette, assay == "culture" & pos_neg == "neg")$hex,
                               "1" = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex)) +  
  scale_color_manual(values = c("Correct" = "black",
                                "Incorrect" = "transparent")) +  
  # Other plot features
  facet_wrap(. ~ factor(error_idx, levels = c("Some errors", "No errors")),
                        ncol = 2, scales = "free_y") +
  labs(y = "individual", x = "day post infection", 
       fill = "Pos/Neg", tag = "A") +
  theme(legend.position = "none", 
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill='white'),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        axis.text = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 10),
        axis.ticks.y= element_blank(),
        text = element_text(size = 11),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        plot.margin = margin(2, 0, 0, 2)) + 
  guides(color = "none") +
  scale_x_continuous(breaks = seq(0, 28, by = 4), limits = c(-0.01, 29),
                     labels = seq(0, 28, by = 4), 
                     expand = c(0, 0)); figSXA



# Plot tRNA vals for diff error types ------------------------------------------

dat$val_total[dat$cens_total == 1] <- runif(length(dat$val_total[dat$cens_total == 1]), -1.95, -0.3)

figSXB <- ggplot(dat) + 
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey90",
            alpha = 0.1, color = "black") +
  geom_point(aes(y = as.numeric(val_total), x = dpi, 
                 color = pos_sg_pred_class, fill = pos_sg_pred_class),
             shape = 21, alpha = 0.8, size = 2) +
  scale_color_manual(values = c("True Positive" = "black",
                                "True Negative" = "black",
                                "False Positive" = "transparent",
                                "False Negative" = "transparent")) +
  scale_fill_manual(values = c("True Positive" = "#E4BC11",
                                "True Negative" = "grey62",
                                "False Positive" = "grey62",
                                "False Negative" = "#E4BC11")) +
  facet_wrap(.~ factor(pos_sg_pred_class, 
                       levels = c("True Positive", "False Negative", 
                                  "False Positive", "True Negative"))) +
  labs(y = "log10 total RNA copies / sample", x = "day post infection",
       tag = "B") +
  scale_x_continuous(breaks = seq(0, 28, 4), limits = c(0, 28)) +
  scale_y_continuous(breaks = seq(-1, 12, by = 2),
                     labels = c("< LOD", seq(1, 12, 2)),
                     limits = c(-2, 10)) +
  theme(legend.position = "none",
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill='white'),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        text = element_text(size = 11),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        strip.background  = element_rect(colour = "black"),
        strip.text.x = element_text(face = "bold"),
        plot.margin = margin(2, 2, 2, 2)); figSXB


dat$val_total[dat$cens_total == 1] <- runif(length(dat$val_total[dat$cens_total == 1]), -0.95, -0.3)
figSXC <- ggplot(dat) + 
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey90",
            alpha = 0.1, color = "black") +
  geom_histogram(aes(y = as.numeric(val_total),
                 color = pos_sg_pred_class, fill = pos_sg_pred_class),
                 breaks = seq(-1, 10, 1)) +
  scale_color_manual(values = c("True Positive" = "black",
                                "True Negative" = "black",
                                "False Positive" = "transparent",
                                "False Negative" = "transparent")) +
  scale_fill_manual(values = c("True Positive" = "#E4BC11",
                               "True Negative" = "grey62",
                               "False Positive" = "grey62",
                               "False Negative" = "#E4BC11")) +
  facet_wrap(.~ factor(pos_sg_pred_class, 
                       levels = c("True Positive", "False Negative", 
                                  "False Positive", "True Negative"))) +
  labs(y = "log10 total RNA copies / sample",
       tag = "C") +
  scale_y_continuous(breaks = seq(-1, 12, by = 2),
                     labels = c("< LOD", seq(1, 12, 2))) +
  theme(legend.position = "none",
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill='white'),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        text = element_text(size = 11),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        strip.background  = element_rect(colour = "black"),
        strip.text.x = element_text(face = "bold"),
        plot.margin = margin(20, 2, 2, 2)); figSXC


# Combine & save ---------------------------------------------------------------

figSX.right <- figSXB + figSXC + plot_layout(nrow = 2, heights = c(1, 1)); figSXB

figSX <- figSXA + figSX.right + plot_layout(ncol = 2, widths = c(1.1, 1)); figSX

ggsave("./outputs/figures/figS14-sgRNA-error-analysis.tiff", 
       figSX, width = 10, height = 10, units = "in")


# Calculate statistics for results ------------------------------------------------

## Percentages -----------------------------------------------------------------

# Percent of individual samples with/out errors
n_traj <- length(indiv_correct_list) + length(indiv_list)
n_traj_no_error <- length(indiv_correct_list)

percent_traj_no_error <- n_traj_no_error / n_traj * 100
cat(round(percent_traj_no_error, 2), "% of trajectories were predicted without error: ",
    n_traj_no_error, "/", n_traj)


# Calculate timing of trajectories' first and last positives
indiv.df <- data.frame(indiv_sample = NA,
                       first_true_pos_dpi = NA,
                       first_pred_pos_dpi = NA,
                       last_true_pos_dpi = NA,
                       last_pred_pos_dpi = NA)

for (indiv.ii in unique(dat.ni$indiv_sample)){
  indiv.sub <- subset(dat.ni, indiv_sample == indiv.ii)
  
  if (nrow(indiv.sub) >= 2) {
    first_true_pos_dpi <- min(indiv.sub$dpi[indiv.sub$pos_sg == 1])
    first_pred_pos_dpi <- min(indiv.sub$dpi[indiv.sub$pos_sg_pred == 1])
    last_true_pos_dpi <- max(indiv.sub$dpi[indiv.sub$pos_sg == 1])
    last_pred_pos_dpi <- max(indiv.sub$dpi[indiv.sub$pos_sg_pred == 1])
    
    indiv.df <- rbind(indiv.df,
                      c(indiv.ii, first_true_pos_dpi, first_pred_pos_dpi,
                        last_true_pos_dpi, last_pred_pos_dpi))
  }
  
}

indiv.df <- indiv.df[-1, ]

n_traj <- nrow(indiv.df) # should be the same as above

# Correct prediction of first positive
n_firstpos <- nrow(subset(indiv.df, first_true_pos_dpi == first_pred_pos_dpi))
percent_firstpos <- n_firstpos / n_traj * 100
cat(round(percent_firstpos, 1), "% of trajectories had their first positive correctly predicted: ",
    n_firstpos, "/", n_traj)


# Correct prediction of last positive
n_lastpos <- nrow(subset(indiv.df, last_true_pos_dpi == last_pred_pos_dpi))
percent_lastpos <- n_lastpos / n_traj * 100
cat(round(percent_lastpos, 1), "% of trajectories had their last positive correctly predicted: ",
    n_lastpos, "/", n_traj)

# Incorrect prediction of last positive
nrow(subset(indiv.df, last_true_pos_dpi != last_pred_pos_dpi))
nrow(subset(indiv.df, last_true_pos_dpi != last_pred_pos_dpi)) / nrow(indiv.df) * 100


## Distribution comparisons ----------------------------------------------------

# Compare sgRNA predicted distribution with totRNA and observed sgRNA distributions

# Load the data in again
dat <- read.csv("./convert-units/manuscript/data/pred-culture-data.csv")
dat <- subset(dat, pos_sg %in% c(0, 1) & cens_total %in% c(0, 1))
dat$val_sg <- as.numeric(dat$val_sg)
dat$val_sg_pred <- as.numeric(dat$val_sg_pred)

# Load the model
stan.model <- cmdstan_model('./code/stan-model-distribution-comparison.stan')

# Prep the data
predicted_vals <- subset(dat, val_sg_pred != -9)$val_sg_pred
observed_vals <- subset(dat, val_sg != -9)$val_sg
totRNA_vals <- subset(dat, pos_total == 1)$val_total

# Convert to format for stan
dat.stan <- c(list(
  
  # Number of observations, articles, & levels 
  N = length(c(predicted_vals, observed_vals, totRNA_vals)),
  N_types = 3,
  y = c(predicted_vals, observed_vals, totRNA_vals), 
  y_type = c(rep(1, length(predicted_vals)),
             rep(2, length(observed_vals)),
             rep(3, length(totRNA_vals)))
))

fit.dist <- stan.model$sample(data = dat.stan,
                              chains = n_chains,
                              parallel_chains = n_chains,
                              iter_warmup = n_iter / 2,
                              iter_sampling = n_iter / 2,
                              refresh = 1000,
                              seed = seed)

pars <- c("mu[1]", "mu[2]", "mu[3]")
pars.df <- fit.dist$draws(variables = pars,
                          inc_warmup = FALSE,
                          format = "df")
pars.df <- as.data.frame(pars.df)

# Histograms comparing distributions
hist(pars.df$`mu[1]` - pars.df$`mu[2]`, main = "Predicted sgRNA vs. Observed sgRNA")
hist(pars.df$`mu[1]` - pars.df$`mu[3]`, main = "Predicted sgRNA vs. totRNA")

# Median difference
median(pars.df$`mu[1]` - pars.df$`mu[2]`)
median(pars.df$`mu[1]` - pars.df$`mu[3]`)

quantile(pars.df$`mu[1]` - pars.df$`mu[2]`, probs = c(0.05, 0.95))
quantile(pars.df$`mu[1]` - pars.df$`mu[3]`, probs = c(0.05, 0.95))


