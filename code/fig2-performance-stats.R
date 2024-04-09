# Creates Figure 2, including the performance statistics of various
#   models, including the simplest, best, and full models, for the sgRNA
#   and culture predictions

# Prep environment -------------------------------------------------------------

# Install & load various packages
req_pkgs <- c("wesanderson", "ggplot2", "cmdstanr", "tidyverse", "ggridges",
              "ggpubr", "patchwork", "stringr")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependences = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))

# Other prep
zissou_pal <- wes_palette("Zissou1", 12, type = "continuous") # color palette
`%notin%` <- Negate(`%in%`)  # for convenience


# Prep results -----------------------------------------------------------------

## sgRNA: Logistic Component -----------------------------------------------------

# Read in the raw selection output table
tbl_log <- read.csv(file = './outputs/tables/tblS2-sgRNA-log-model-selection-inf-raw.csv', 
                    stringsAsFactors = FALSE)
tbl_log$n_predictors <- as.numeric(substring(tbl_log$model, 2, 2)) # add # predictors

# Get best performing models for each # of predictors based on ELPD
top_log <- c("l1", "l8.1")
for (n_pred in 2:7){
  pred.sub <- subset(tbl_log, n_predictors == n_pred)
  top_model <-tbl_log$model[tbl_log$elpd == max(pred.sub$elpd)]
  top_log <- c(top_log, top_model)
}

# Subset to best performing models per predictor number
tbl_log <- subset(tbl_log, model %in% top_log)

# Add "T" to indicate totRNA included in all models
tbl_log$predictors[tbl_log$model != "l1"] <- paste0("T + ", tbl_log$predictors[tbl_log$model != "l1"])
tbl_log$predictors[tbl_log$model == "l1"] <- "T"


## sgRNA: Linear Component -----------------------------------------------------

tbl_lin <- read.csv(file = './outputs/tables/tblS3-sgRNA-full-model-selection-inf-raw.csv', 
                    stringsAsFactors = FALSE)
tbl_lin$n_predictors <- as.numeric(substring(tbl_lin$model, 2, 2))

# Get best performing models for each # of predictors, based on ELPD
top_lin <- c("f1", "f8.1")
for (n_pred in 2:7){
  pred.sub <- subset(tbl_lin, n_predictors == n_pred)
  top_model <-tbl_lin$model[tbl_lin$elpd == max(pred.sub$elpd)]
  top_lin <- c(top_lin, top_model)
}

tbl_lin <- subset(tbl_lin, model %in% top_lin)
tbl_lin$predictors[tbl_lin$model != "f1"] <- paste0("T + ", tbl_lin$predictors[tbl_lin$model != "f1"])
tbl_lin$predictors[tbl_lin$model == "f1"] <- "T"

tbl_lin$mae <- as.numeric(tbl_lin$mae_test)
tbl_lin$Predictors <- tbl_lin$predictors
tbl_lin$Model <- tbl_lin$model


## Culture: Logistic Component -------------------------------------------------

tbl_cul <- read.csv(file = './outputs/tables/tblS4-culture-totRNA-model-selection-inf-raw.csv', 
                    stringsAsFactors = FALSE)
tbl_cul$n_predictors <- as.numeric(substring(tbl_cul$model, 2, 2))
tbl_cul$n_predictors[tbl_cul$model == "c10.1"] <- 10

# Get best performing models for each # of predictors, based on ELPD
top_log <- c("c1", "c10.1")
for (n_pred in 2:9){
  pred.sub <- subset(tbl_cul, n_predictors == n_pred)
  top_model <-tbl_cul$model[tbl_cul$elpd == max(pred.sub$elpd)]
  top_log <- c(top_log, top_model)
}

tbl_cul <- subset(tbl_cul, model %in% top_log)
tbl_cul$predictors <- str_remove_all(tbl_cul$predictors, "_idx")
tbl_cul$predictors <- str_remove_all(tbl_cul$predictors, "log10_|_pfu")
tbl_cul$predictors <- str_replace_all(tbl_cul$predictors, ",", " +" )
tbl_cul$predictors <- toupper(tbl_cul$predictors)
tbl_cul$predictors[tbl_cul$model != "c1"] <- paste0("T + ", tbl_cul$predictors[tbl_cul$model != "c1"])
tbl_cul$predictors[tbl_cul$model == "c1"] <- "T"


# Get cross-model statistics ---------------------------------------------------

# Need to make sure these aren't generated for the subset of models

range(abs(tbl_log$percent_test - tbl_log$percent_train))
median(abs(tbl_log$percent_test - tbl_log$percent_train))

range(abs(tbl_lin$mae_test - tbl_lin$mae_train))
median(abs(tbl_lin$mae_test - tbl_lin$mae_train))

range(abs(tbl_cul$percent_test_adj - tbl_cul$percent_train_adj))
median(abs(tbl_cul$percent_test_adj - tbl_cul$percent_train_adj))


# Create a legend --------------------------------------------------------------

color.df <- data.frame(model = c("Simple", "Best", "Full"),
                       color = c("#0F4C5E", "#85B068", "#5e468e"))
color.df$model <- factor(color.df$model, levels = color.df$model)

fig_legend <- ggplot(color.df) +
  geom_point(aes(x = 1, y = color, fill = model), shape = 21) +
  guides(fill = guide_legend(override.aes = list(size = 2.2))) +
  labs(fill = "Model") +
  scale_fill_manual(values = color.df$color) +
  theme(legend.margin = margin(0.1,0.1,0.1,0.1, unit="cm"),
        legend.background = element_rect(colour = "grey35"),
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank()); fig_legend


# 2A: Logistic -----------------------------------------------------------------

## ELPD ------------------------------------------------------------------------

tbl_log$stat_label <- "ELPD"
fig_log_elpd <- ggplot(data = tbl_log) +
  
  # Reference lines
  geom_hline(yintercept = tbl_log$elpd[tbl_log$n_predictors == 8], color = "#5e468e",
             linetype = "dashed") +
  geom_vline(xintercept = 4, color = "#85B068", alpha = 0.7) +
  
  # Results
  geom_line(aes(x = n_predictors, y = elpd), color = "grey") +
  geom_point(data = subset(tbl_log, n_predictors %notin% c(1, 4, 8)),
             aes(x = n_predictors, y = elpd), 
             size = 1.5, shape = 21, fill = "grey") +
  geom_point(data = subset(tbl_log, n_predictors == 1),
             aes(x = n_predictors, y = elpd), 
             size = 2.3, shape = 21, fill = "#0F4C5E") +
  geom_point(data = subset(tbl_log, n_predictors == 4),
             aes(x = n_predictors, y = elpd), 
             size = 2.3, shape = 21, fill = "#85B068") +
  geom_point(data = subset(tbl_log, n_predictors == 8),
             aes(x = n_predictors, y = elpd), 
             size = 2.3, shape = 21, fill = "#5e468e") +
  
  # Axes
  labs(y = element_blank(), x = element_blank()) +
  ggtitle("sgRNA Logistic") +
  scale_y_continuous(limits = c(-315, -246), breaks = seq(-310, -250, 20)) +
  
  # Other aesthetics
  facet_wrap(.~ stat_label) +
  theme(legend.position = "none", 
        axis.text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        text = element_text(size = 11),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        plot.margin = margin(2, 2, 2, 2)) +
    annotation_custom(get_legend(fig_legend), 
                    xmin  = 13, xmax = 1); fig_log_elpd  



## Prediction Accuracy ---------------------------------------------------------

tbl_log$stat_label <- "Prediction Accuracy (%)"
fig_log_pred <- ggplot(data = tbl_log) +
  
  # Reference lines
  geom_hline(yintercept = tbl_log$percent_test[tbl_log$n_predictors == 8], 
             color = "#5e468e", linetype = "dashed") +
  geom_vline(xintercept = 4, color = "#85B068", alpha = 0.7) +
  
  # Results
  geom_line(aes(x = n_predictors, y = percent_test), color = "grey") +
  geom_point(data = subset(tbl_log, n_predictors %notin% c(1, 4, 8)),
             aes(x = n_predictors, y = percent_test), 
             size = 1.5, shape = 21, fill = "grey") +
  geom_point(data = subset(tbl_log, n_predictors == 1),
             aes(x = n_predictors, y = percent_test), 
             size = 2.3, shape = 21, fill = "#0F4C5E") +
  geom_point(data = subset(tbl_log, n_predictors == 4),
             aes(x = n_predictors, y = percent_test), 
             size = 2.3, shape = 21, fill = "#85B068") +
  geom_point(data = subset(tbl_log, n_predictors == 8),
             aes(x = n_predictors, y = percent_test), 
             size = 2.3, shape = 21, fill = "#5e468e") +
  
  # Axes
  labs(y = element_blank(), x = element_blank()) +
  scale_y_continuous(limits = c(87, 91.3), breaks = seq(87, 92, 2)) +
  
  # Other aesthetics
  facet_wrap(.~ stat_label) +
  theme(legend.position = "none", 
        axis.text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_blank(),
        text = element_text(size = 11),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        plot.margin = margin(2, 2, 2, 2)); fig_log_pred


## MCC -------------------------------------------------------------------------

tbl_log$stat_label <- "MCC"
fig_log_mcc <- ggplot(tbl_log) + 
  
  # Reference lines
  geom_hline(yintercept = tbl_log$MCC[tbl_log$n_predictors == 8], 
             color = "#5e468e", linetype = "dashed") +
  geom_vline(xintercept = 4, color = "#85B068", alpha = 0.7) +
  
  # Results
  geom_line(aes(x = n_predictors, y = MCC), color = "grey") +
  geom_point(data = subset(tbl_log, n_predictors %notin% c(1, 4, 8)),
             aes(x = n_predictors, y = MCC),
             size = 1.5, shape = 21, fill = "grey") +
  geom_point(data = subset(tbl_log, n_predictors == 1),
             aes(x = n_predictors, y = MCC),
             size = 2.3, shape = 21, fill = "#0F4C5E") +
  geom_point(data = subset(tbl_log, n_predictors == 4),
             aes(x = n_predictors, y = MCC),
             size = 2.3, shape = 21, fill = "#85B068") +
  geom_point(data = subset(tbl_log, n_predictors == 8),
             aes(x = n_predictors, y = MCC),
             size = 2.3, shape = 21, fill = "#5e468e") +
  
  # Axes
  scale_x_continuous(breaks = 1:8, 
                     labels = tbl_log[order(tbl_log$n_predictors), ]$predictors) +
  scale_y_continuous(limits = c(0.695, 0.84), breaks = seq(0.7, 0.84, 0.1),
                     expand = c(0, 0)) +
  labs(y = element_blank(), x = "") +
  
  # Other aesthetics
  facet_wrap(.~ stat_label) +
  theme(legend.position = "none", 
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = c("#0F4C5E", "black",
                                                                                "black", "#85B068", 
                                                                                "black", "black",
                                                                                "black", "#5e468e"),
                                   face = c("bold", "plain",
                                            "plain", "bold", 
                                            "plain", "plain",
                                            "plain", "bold")),
        axis.title = element_text(size = 10),
        text = element_text(size = 11),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor.x = element_blank(),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        plot.margin = margin(2, 2, 2, 2)); fig_log_mcc


## Add a panel for F-score
fig_log_elpd + fig_log_pred + fig_log_mcc + 
  plot_layout(nrow = 3, heights = c(1, 1, 1))


# 2B: Linear model -------------------------------------------------------------

## ELPD ------------------------------------------------------------------------

tbl_lin$stat_label <- "ELPD"

fig_lin_elpd <- ggplot(data = tbl_lin) +
  # Reference lines
  geom_hline(yintercept = tbl_lin$elpd[tbl_lin$n_predictors == 8], 
             color = "#5e468e", linetype = "dashed") +
  geom_vline(xintercept = 5, color = "#85B068", alpha = 0.7) +
  
  # Results points & lines
  geom_line(aes(x = n_predictors, y = elpd), color = "grey") +
  geom_point(data = subset(tbl_lin, n_predictors %notin% c(1, 5, 8)),
             aes(x = n_predictors, y = elpd), 
             size = 1.5, shape = 21, fill = "grey") +
  geom_point(data = subset(tbl_lin, n_predictors == 1),
             aes(x = n_predictors, y = elpd), 
             size = 2.3, shape = 21, fill = "#0F4C5E") +
  geom_point(data = subset(tbl_lin, n_predictors == 5),
             aes(x = n_predictors, y = elpd), 
             size = 2.3, shape = 21, fill = "#85B068") +
  geom_point(data = subset(tbl_lin, n_predictors == 8),
             aes(x = n_predictors, y = elpd), 
             size = 2.3, shape = 21, fill = "#5e468e") +
  
  # Facet for the convenient strip  label
  facet_wrap(.~ stat_label) +
  
  labs(y = element_blank(), x = element_blank()) +
  ggtitle("sgRNA Linear") +
  scale_y_continuous(limits = c(-1050, -810), 
                     breaks = seq(-1050, -800, 100)) +
  theme(legend.position = "none", 
        axis.text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        text = element_text(size = 11),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        plot.margin = margin(2, 2, 2, 2)); fig_lin_elpd  


## % within the 50%PI ----------------------------------------------------------

tbl_lin$stat_label <- "% within 50% PI"
fig_lin_pred <- ggplot(data = tbl_lin) +
  geom_hline(yintercept = tbl_lin$percent_50_test[tbl_lin$n_predictors == 8], 
             color = "#5e468e", linetype = "dashed") +
  geom_vline(xintercept = 5, color = "#85B068", alpha = 0.7) +
  
  # Results points & lines
  geom_line(aes(x = n_predictors, y = percent_50_test), color = "grey") +
  geom_point(data = subset(tbl_lin, n_predictors %notin% c(1, 5, 8)),
             aes(x = n_predictors, y = percent_50_test), 
             size = 1.5, shape = 21, fill = "grey") +
  geom_point(data = subset(tbl_lin, n_predictors == 1),
             aes(x = n_predictors, y = percent_50_test), 
             size = 2.3, shape = 21, fill = "#0F4C5E") +
  geom_point(data = subset(tbl_lin, n_predictors == 5),
             aes(x = n_predictors, y = percent_50_test), 
             size = 2.3, shape = 21, fill = "#85B068") +
  geom_point(data = subset(tbl_lin, n_predictors == 8),
             aes(x = n_predictors, y = percent_50_test), 
             size = 2.3, shape = 21, fill = "#5e468e") +
  
  # Facet for the convenient strip  label
  facet_wrap(.~ stat_label) +
  
  labs(y = element_blank(), x = element_blank()) +
  scale_x_continuous(breaks = 1:8,
                     labels = tbl_lin[order(tbl_lin$n_predictors), ]$Predictors) +
  scale_y_continuous(limits = c(47, 57.5), 
                     breaks = seq(48, 58, 4)) +
  theme(legend.position = "none", 
        axis.text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_blank(),
        text = element_text(size = 11),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor.x = element_blank(),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        plot.margin = margin(2, 2, 2, 2)); fig_lin_pred


## MAE -------------------------------------------------------------------------

tbl_lin$stat_label <- "MAE"
fig_lin_mae <- ggplot(data = tbl_lin) +
  geom_hline(yintercept = tbl_lin$mae_test[tbl_lin$n_predictors == 8], 
             color = "#5e468e", linetype = "dashed") +
  geom_vline(xintercept = 5, color = "#85B068", alpha = 0.7) +
  
  geom_line(aes(x = n_predictors, y = mae_test), color = "grey") +
  geom_point(data = subset(tbl_lin, n_predictors %notin% c(1, 5, 8)),
             aes(x = n_predictors, y = mae_test), 
             size = 1.5, shape = 21, fill = "grey") +
  geom_point(data = subset(tbl_lin, n_predictors == 1),
             aes(x = n_predictors, y = mae_test), 
             size = 2.3, shape = 21, fill = "#0F4C5E") +
  geom_point(data = subset(tbl_lin, n_predictors == 5),
             aes(x = n_predictors, y = mae_test), 
             size = 2.3, shape = 21, fill = "#85B068") +
  geom_point(data = subset(tbl_lin, n_predictors == 8),
             aes(x = n_predictors, y = mae_test), 
             size = 2.3, shape = 21, fill = "#5e468e") +
  
  # Facet for the convenient strip  label
  facet_wrap(.~ stat_label) +
  
  labs(y = element_blank(), x = element_blank()) +
  scale_x_continuous(breaks = 1:8,
                     labels = tbl_lin[order(tbl_lin$n_predictors), ]$predictors) +
  scale_y_continuous(limits = c(0.4, 0.6),
                     breaks = seq(0.4, 0.6, 0.1)) +
  theme(legend.position = "none", 
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                                   colour = c("#0F4C5E", rep("black", 3),
                                              "#85B068", rep("black", 2),
                                              "#5e468e"),
                                   face = c("bold", rep("plain", 3),
                                            "bold", rep("plain", 2),
                                            "bold")),
        axis.title = element_text(size = 10),
        axis.title.y = element_blank(),
        text = element_text(size = 11),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor.x = element_blank(),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        plot.margin = margin(2, 2, 2, 2)); fig_lin_mae

fig_lin_elpd + fig_lin_pred + fig_lin_mae + plot_layout(nrow = 3, heights = c(1, 1))


# 2C: Culture Logistic ---------------------------------------------------------

## ELPD ------------------------------------------------------------------------

tbl_cul$stat_label <- "ELPD"
fig_cul_elpd <- ggplot(data = tbl_cul) +
  geom_hline(yintercept = tbl_cul$elpd[tbl_cul$n_predictors == 10], 
             color = "#5e468e",
             linetype = "dashed") +
  geom_vline(xintercept = 8, color = "#85B068", alpha = 0.7) +
  geom_line(aes(x = n_predictors, y = elpd), color = "grey") +
  geom_point(data = subset(tbl_cul, n_predictors %notin% c(1, 8, 10)),
             aes(x = n_predictors, y = elpd), 
             size = 1.5, shape = 21, fill = "grey") +
  geom_point(data = subset(tbl_cul, n_predictors == 1),
             aes(x = n_predictors, y = elpd), 
             size = 2.3, shape = 21, fill = "#0F4C5E") +
  geom_point(data = subset(tbl_cul, n_predictors == 8),
             aes(x = n_predictors, y = elpd), 
             size = 2.3, shape = 21, fill = "#85B068") +
  geom_point(data = subset(tbl_cul, n_predictors == 10),
             aes(x = n_predictors, y = elpd), 
             size = 2.3, shape = 21, fill = "#5e468e") +
  
  # Facet for the convenient strip  label
  facet_wrap(.~ stat_label) +
  ggtitle("Culture Logistic") +
  labs(y = element_blank(), x = element_blank()) +
  scale_y_continuous(limits = c(-405, -339), breaks = seq(-400, -340, 20)) +
  theme(legend.position = "none", 
        axis.text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        text = element_text(size = 11),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        plot.margin = margin(2, 2, 2, 2)); fig_cul_elpd  


## Prediction Accuracy ---------------------------------------------------------

tbl_cul$stat_label <- "Prediction Accuracy (%)"
fig_cul_pred <- ggplot(data = tbl_cul) +
  # Reference lines
  geom_hline(yintercept = tbl_cul$percent_test_adj[tbl_cul$n_predictors == 10], 
             color = "#5e468e", linetype = "dashed") +
  geom_vline(xintercept = 8, color = "#85B068", alpha = 0.7) +
  
  # Results points & lines
  geom_line(aes(x = n_predictors, y = percent_test_adj), color = "grey") +
  geom_point(data = subset(tbl_cul, n_predictors %notin% c(1, 8, 10)),
             aes(x = n_predictors, y = percent_test_adj), 
             size = 1.5, shape = 21, fill = "grey") +
  geom_point(data = subset(tbl_cul, n_predictors == 1),
             aes(x = n_predictors, y = percent_test_adj), 
             size = 2.3, shape = 21, fill = "#0F4C5E") +
  geom_point(data = subset(tbl_cul, n_predictors == 8),
             aes(x = n_predictors, y = percent_test_adj), 
             size = 2.3, shape = 21, fill = "#85B068") +
  geom_point(data = subset(tbl_cul, n_predictors == 10),
             aes(x = n_predictors, y = percent_test_adj), 
             size = 2.3, shape = 21, fill = "#5e468e") +
  
  
  # Facet for the convenient strip  label
  facet_wrap(.~ stat_label) +
  
  labs(y = element_blank(), x = element_blank()) +
  scale_y_continuous(limits = c(81.5, 86.5), breaks = seq(82, 96, 2),
                     expand = c(0, 0)) +
  theme(legend.position = "none", 
        axis.text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_blank(),
        text = element_text(size = 11),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        plot.margin = margin(2, 2, 2, 2)); fig_cul_pred



## MCC -------------------------------------------------------------------------

tbl_cul$stat_label <- "MCC"
fig_cul_mcc <- ggplot(tbl_cul) + 
  # Reference lines
  geom_hline(yintercept = tbl_cul$MCC[tbl_cul$n_predictors == 10], 
             color = "#5e468e", linetype = "dashed") +
  geom_vline(xintercept = 8, color = "#85B068", alpha = 0.7) +
  
  # Results points & lines
  geom_line(aes(x = n_predictors, y = MCC), color = "grey") +
  geom_point(data = subset(tbl_cul, n_predictors %notin% c(1, 8, 10)),
             aes(x = n_predictors, y = MCC),
             size = 1.5, shape = 21, fill = "grey") +
  geom_point(data = subset(tbl_cul, n_predictors == 1),
             aes(x = n_predictors, y = MCC),
             size = 2.3, shape = 21, fill = "#0F4C5E") +
  geom_point(data = subset(tbl_cul, n_predictors == 8),
             aes(x = n_predictors, y = MCC),
             size = 2.3, shape = 21, fill = "#85B068") +
  geom_point(data = subset(tbl_cul, n_predictors == 10),
             aes(x = n_predictors, y = MCC),
             size = 2.3, shape = 21, fill = "#5e468e") +
  
  # Facet for the convenient strip  label
  facet_wrap(.~ stat_label) +
  
  scale_x_continuous(breaks = 1:10, 
                     labels = tbl_cul[order(tbl_cul$n_predictors), ]$predictors) +
  scale_y_continuous(limits = c(0.45, 0.62), breaks = seq(0.5, 0.6, 0.1)) +
  labs(x = element_blank(), y = element_blank()) +
  theme(legend.position = "none", 
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                                   colour = c("#0F4C5E", rep("black", 6),
                                              "#85B068", "black",
                                              "#5e468e"),
                                   face = c("bold", rep("plain", 6),
                                            "bold", "plain",
                                            "bold")),
        axis.title = element_text(size = 10),
        axis.title.y = element_blank(),
        text = element_text(size = 11),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor.x = element_blank(),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        plot.margin = margin(2, 2, 2, 2)); fig_cul_mcc


## Add a panel for F-score
fig_cul_elpd + fig_cul_pred + fig_cul_mcc + 
  plot_layout(nrow = 3, heights = c(1, 1, 1))



# Combine into 1 plot -------------------------------------------------------------


fig2 <- fig_log_elpd + fig_log_pred + fig_log_mcc + 
         fig_lin_elpd + fig_lin_pred + fig_lin_mae + 
         fig_cul_elpd + fig_cul_pred + fig_cul_mcc +
         plot_layout(nrow = 3, ncol = 3, byrow = FALSE, heights = c(1, 1, 1),
                     widths = c(8, 8, 9.5)); fig2



# Save -------------------------------------------------------------------------

ggsave('./outputs/figures/fig2-performance-stats.tiff',
       plot = fig2,
       device = 'tiff',
       height = 8,
       width = 7.5,
       units = 'in',
       bg = 'white')



