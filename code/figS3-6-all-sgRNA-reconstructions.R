# Generates supplemental figures 3 through 6
#   Plots all individual trajectories with sgRNA data, stratified by organ system
# No dependencies, this file will run without running any other scripts

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



# Prep data --------------------------------------------------------------------

dat.fig <- read.csv("./data/pred-sg-data.csv")

# Removes NA & positively censored data for this analysis
dat.fig <- subset(dat.fig, cens_sg %in% c(0, 1) & cens_total %in% c(0, 1)) 

# Convert columns to numerics
dat.fig$val_sg <- as.numeric(dat.fig$val_sg)
dat.fig$val_total <- as.numeric(dat.fig$val_total)

# Combine indiv names with sample rep to get the right connecting lines
dat.fig$indiv_sample <- paste0(dat.fig$indiv, "_", dat.fig$sample_rep, 
                               "_", dat.fig$tg_idx)

# Indicator to simplify plotting
dat.fig$sample_sg <- paste0("sgRNA: ", dat.fig$organ_system)
dat.fig$sample_total <- paste0("total RNA: ", dat.fig$organ_system)
dat.fig$sample_culture <- paste0("culture: ", dat.fig$organ_system)

# Set values for samples below LOD 
dat.fig$val_sg_pred[dat.fig$pos_sg_pred == 0] <- -1
dat.fig$val_sg[dat.fig$cens_sg == 1] <- -0.75
dat.fig$val_total[dat.fig$cens_total == 1] <- -0.25


# Calculate some basic stats ---------------------------------------------------

nrow(subset(dat.fig, pos_sg == pos_sg_pred))/nrow(dat.fig) # prediction accuracy once fit to all of the data



# URT --------------------------------------------------------------------

dat.fig.urt <- subset(dat.fig, sample_type == "Non-invasive" & organ_system == "URT")

indiv_list <- c()
for (indiv.ii in unique(dat.fig.urt$indiv_sample)) {
  indiv.sub <- subset(dat.fig.urt, indiv_sample == indiv.ii)
  if (nrow(indiv.sub) >= 2){
    indiv_list <- c(indiv_list, indiv.ii)
  }
}

dat.fig.urt <- subset(dat.fig.urt, indiv_sample %in% indiv_list)
dat.fig.urt$sample_rep[dat.fig.urt$sample_rep == "Nasopharyngeal Swab"] <- "Nasoph. Swab"
dat.fig.urt$panel_label <-paste0(dat.fig.urt$indiv, "\n", dat.fig.urt$sample_rep)

fig.sg.traj.urt <- ggplot(dat.fig.urt) + 
  
  geom_text(aes(label = panel_label, x = 10.5, y = 8), size = 1.8) +
  
  geom_hline(data = subset(dat.fig.urt, lodq_total == -9),
             aes(yintercept = 0), linetype = "dotted",
             size = 0.5, color = "grey22", alpha = 0.7) +
  geom_hline(data = subset(dat.fig.urt, lodq_total != -9),
             aes(yintercept = lodq_total), linetype = "dotted",
             size = 0.5, color = "#0A2472", alpha = 0.7) +
  geom_hline(data = subset(dat.fig.urt, lodq_total != -9),
             aes(yintercept = lodq_sg), linetype = "dotted",
             size = 0.5, color = zissou_pal[1], alpha = 0.7) +
  
  geom_line(aes(x = dpi, y = val_sg, group = indiv_sample,
                color = "sgRNA (observed)"),
            alpha = 0.8) +
  geom_point(aes(x = dpi, y = val_sg, 
                 group = indiv_sample, 
                 fill = "sgRNA (observed)"),
             alpha = 0.8, size = 1.5, shape = 23) +
  
  geom_line(aes(x = dpi, y = val_sg_pred, group = indiv_sample,
                color = "sgRNA (predicted)"),
            alpha = 0.7) +
  geom_point(aes(x = dpi, y = val_sg_pred, fill = "sgRNA (predicted)",
                 group = indiv_sample),
             alpha = 0.8, size = 1.5, shape = 24) +
  
  geom_line(aes(x = dpi, y = val_total, group = indiv_sample,
                color = "total RNA (observed)"), 
            alpha = 0.7) +
  geom_point(aes(x = dpi, y = val_total, 
                 group = indiv_sample,
                 fill = "total RNA (observed)"),
             alpha = 1, size = 1.5, shape = 21) +
  
  scale_fill_manual(values = c("total RNA (observed)" = "#0A2472",
                               "sgRNA (observed)" = zissou_pal[1], #"#3B9AB2",
                               "sgRNA (predicted)" = "#0A8549"),
                    breaks = c("total RNA (observed)", "sgRNA (observed)", 
                               "sgRNA (predicted)")) +
  scale_color_manual(values = c("total RNA (observed)" = "#0A2472",
                                "sgRNA (observed)" = zissou_pal[1], #"#3B9AB2",
                                "sgRNA (predicted)" = "#0A8549")) +
  facet_wrap(indiv_sample ~., ncol = 10) + 
  labs(y = "log10 RNA copies / sample", x = "day post infection",
       fill = "URT") +
  theme(
    legend.position = "bottom", 
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.background = element_rect(fill='white'),
    legend.key=element_blank(),
    legend.key.size = unit(0.4, 'cm'),
    legend.box.spacing = unit(1, "pt"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    text = element_text(size = 7),
    strip.background = element_rect(colour = "black"),
    strip.text = element_blank(),
    panel.spacing = unit(0.1, "lines"), 
    panel.background = element_rect(fill = "white",
                                    colour = "black",
                                    size = 1, linetype = "solid"),
    panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                    colour = "light grey"), 
    panel.grid.minor = element_line(size = 0.02, linetype = 'solid',
                                    colour = "light grey"),
    plot.margin = margin(2, 2, 0, 2)) + 
  guides(color = "none", 
         fill = guide_legend(override.aes = list(size = 1.8, shape = c(21, 23, 24)))) +
  scale_x_continuous(breaks = seq(0, 16, by = 4), limits = c(0, 16),
                     labels = c(seq(0, 12, by = 4), ""), 
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(from = 0, to = 10, by = 2), 
                     minor_breaks = seq(1, 9, 2),
                     limits = c(-1.5, 10),
                     expand = c(0, 0)); fig.sg.traj.urt

ggsave("./outputs/figures/figS3-sgRNA-URT-traj-reconstructions.tiff", 
       fig.sg.traj.urt, width = 7.5, height = 10, units = "in")



# LRT --------------------------------------------------------------------------

dat.fig.lrt <- subset(dat.fig, sample_type == "Non-invasive" & organ_system == "LRT")

indiv_list <- c()
for (indiv.ii in unique(dat.fig.lrt$indiv_sample)) {
  indiv.sub <- subset(dat.fig.lrt, indiv_sample == indiv.ii)
  if (nrow(indiv.sub) >= 2){
    indiv_list <- c(indiv_list, indiv.ii)
  }
}

dat.fig.lrt <- subset(dat.fig.lrt, indiv_sample %in% indiv_list)
dat.fig.lrt$panel_label <-paste0(dat.fig.lrt$indiv, "\n", dat.fig.lrt$sample_rep)


fig.sg.traj.lrt <- ggplot(dat.fig.lrt) + 
  
  geom_text(aes(label = panel_label, x = 10.5, y = 8), size = 1.8) +
  
  geom_hline(data = subset(dat.fig.lrt, lodq_total == -9),
             aes(yintercept = 0), linetype = "dotted",
             size = 0.5, color = "grey22", alpha = 0.7) +
  geom_hline(data = subset(dat.fig.lrt, lodq_total != -9),
             aes(yintercept = lodq_total), linetype = "dotted",
             size = 0.5, color = zissou_pal[12], alpha = 0.7) +
  geom_hline(data = subset(dat.fig.lrt, lodq_sg != -9),
             aes(yintercept = lodq_sg), linetype = "dotted",
             size = 0.5, color = zissou_pal[10], alpha = 0.7) +
  
  geom_line(aes(x = dpi, y = val_sg, group = indiv_sample,
                color = "sgRNA (observed)"),
            alpha = 0.8) +
  geom_point(aes(x = dpi, y = val_sg, 
                 group = indiv_sample, 
                 fill = "sgRNA (observed)"),
             alpha = 0.8, size = 1.5, shape = 23) +
  
  geom_line(aes(x = dpi, y = val_sg_pred, group = indiv_sample,
                color = "sgRNA (predicted)"),
            alpha = 0.7) +
  geom_point(aes(x = dpi, y = val_sg_pred, fill = "sgRNA (predicted)",
                 group = indiv_sample),
             alpha = 0.8, size = 1.5, shape = 24) +
  
  geom_line(aes(x = dpi, y = val_total, group = indiv_sample,
                color = "total RNA (observed)"), 
            alpha = 0.7) +
  geom_point(aes(x = dpi, y = val_total, 
                 group = indiv_sample,
                 fill = "total RNA (observed)"),
             alpha = 1, size = 1.5, shape = 21) +
  
  scale_fill_manual(values = c("total RNA (observed)" = zissou_pal[12],
                               "sgRNA (observed)" = zissou_pal[10], #"#3B9AB2",
                               "sgRNA (predicted)" = zissou_pal[7]),
                    breaks = c("total RNA (observed)", "sgRNA (observed)", 
                               "sgRNA (predicted)")) +
  scale_color_manual(values = c("total RNA (observed)" = zissou_pal[12],
                                "sgRNA (observed)" = zissou_pal[10],  #"#3B9AB2",
                                "sgRNA (predicted)" = zissou_pal[7])) +
  
  facet_wrap(indiv_sample ~., ncol = 10, nrow = 14) + 
  labs(y = "log10 RNA copies / sample", x = "day post infection",
       fill = "LRT") +
  theme(
    legend.position = "bottom", 
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.background = element_rect(fill='white'),
    legend.key=element_blank(),
    legend.key.size = unit(0.4, 'cm'),
    legend.box.spacing = unit(1, "pt"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    text = element_text(size = 7),
    strip.background = element_rect(colour = "black"),
    strip.text = element_blank(),
    panel.spacing = unit(0.1, "lines"), 
    panel.background = element_rect(fill = "white",
                                    colour = "black",
                                    size = 1, linetype = "solid"),
    panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                    colour = "light grey"), 
    panel.grid.minor = element_line(size = 0.02, linetype = 'solid',
                                    colour = "light grey"),
    plot.margin = margin(2, 2, 2, 2)) + 
  guides(color = "none", 
         fill = guide_legend(override.aes = list(size = 1.8, shape = c(21, 23, 24)))) +
  scale_x_continuous(breaks = seq(0, 16, by = 4), limits = c(0, 16),
                     labels = c(seq(0, 12, by = 4), ""), 
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(from = 0, to = 10, by = 2), 
                     minor_breaks = seq(1, 9, 2),
                     limits = c(-1.5, 10),
                     expand = c(0, 0)); fig.sg.traj.lrt

ggsave("./outputs/figures/figS4-sgRNA-LRT-traj-reconstructions.tiff", 
       fig.sg.traj.lrt, width = 7.5, height = 10 * 9/14, units = "in")


# GI & Other --------------------------------------------------------------------

dat.fig.other <- subset(dat.fig, sample_type == "Non-invasive" & organ_system %in% c("GI", "Other"))

indiv_list <- c()
for (indiv.ii in unique(dat.fig.other$indiv_sample)) {
  indiv.sub <- subset(dat.fig.other, indiv_sample == indiv.ii)
  if (nrow(indiv.sub) >= 2){
    indiv_list <- c(indiv_list, indiv.ii)
  }
}

dat.fig.other <- subset(dat.fig.other, indiv_sample %in% indiv_list)
dat.fig.other$panel_label <-paste0(dat.fig.other$indiv, "\n", dat.fig.other$sample_rep)


fig.sg.traj.other <- ggplot(dat.fig.other) + 
  
  geom_text(aes(label = panel_label, x = 9, y = 8), size = 1.8) +
  
  geom_hline(data = subset(dat.fig.other, lodq_total == -9),
             aes(yintercept = 0), linetype = "dotted",
             size = 0.5, color = "grey22", alpha = 0.7) +
  geom_hline(data = subset(dat.fig.other, lodq_total != -9),
             aes(yintercept = lodq_total), linetype = "dotted",
             size = 0.5, color = zissou_pal[12], alpha = 0.7) +
  geom_hline(data = subset(dat.fig.other, lodq_sg != -9),
             aes(yintercept = lodq_sg), linetype = "dotted",
             size = 0.5, color = zissou_pal[10], alpha = 0.7) +
  
  geom_line(aes(x = dpi, y = val_sg, group = indiv_sample,
                color = "sgRNA (observed)"),
            alpha = 0.8) +
  geom_point(aes(x = dpi, y = val_sg, 
                 group = indiv_sample, 
                 fill = "sgRNA (observed)"),
             alpha = 0.8, size = 1.5, shape = 23) +
  
  geom_line(aes(x = dpi, y = val_sg_pred, group = indiv_sample,
                color = "sgRNA (predicted)"),
            alpha = 0.7) +
  geom_point(aes(x = dpi, y = val_sg_pred, fill = "sgRNA (predicted)",
                 group = indiv_sample),
             alpha = 0.8, size = 1.5, shape = 24) +
  
  geom_line(aes(x = dpi, y = val_total, group = indiv_sample,
                color = "total RNA (observed)"), 
            alpha = 0.7) +
  geom_point(aes(x = dpi, y = val_total, 
                 group = indiv_sample,
                 fill = "total RNA (observed)"),
             alpha = 1, size = 1.5, shape = 21) +
  
  scale_fill_manual(values = c("total RNA (observed)" = "#710193",
                               "sgRNA (observed)" = "#C36F8E", #"#3B9AB2",
                               "sgRNA (predicted)" = "light pink"),
                    breaks = c("total RNA (observed)", "sgRNA (observed)", 
                               "sgRNA (predicted)")) +
  scale_color_manual(values = c("total RNA (observed)" = "#710193",
                                "sgRNA (observed)" = "#C36F8E",  #"#3B9AB2",
                                "sgRNA (predicted)" = "light pink")) +
  
  facet_wrap(indiv_sample ~., ncol = 10, nrow = 14) + 
  labs(y = "log10 RNA copies / sample", x = "day post infection",
       fill = "GI & Other") +
  theme(
    legend.position = "bottom", 
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.background = element_rect(fill='white'),
    legend.key=element_blank(),
    legend.key.size = unit(0.4, 'cm'),
    legend.box.spacing = unit(1, "pt"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    text = element_text(size = 7),
    strip.background = element_rect(colour = "black"),
    strip.text = element_blank(),
    panel.spacing = unit(0.1, "lines"), 
    panel.background = element_rect(fill = "white",
                                    colour = "black",
                                    size = 1, linetype = "solid"),
    panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                    colour = "light grey"), 
    panel.grid.minor = element_line(size = 0.02, linetype = 'solid',
                                    colour = "light grey"),
    plot.margin = margin(4, 2, 2, 2)) + 
  guides(color = "none", 
         fill = guide_legend(override.aes = list(size = 1.8, shape = c(21, 23, 24)))) +
  scale_x_continuous(breaks = seq(0, 16, by = 4), limits = c(0, 16),
                     labels = c(seq(0, 12, by = 4), ""), 
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(from = 0, to = 10, by = 2), 
                     minor_breaks = seq(1, 9, 2),
                     limits = c(-1.5, 10),
                     expand = c(0, 0)); fig.sg.traj.other

ggsave("./outputs/figures/figS5-sgRNA-Other-traj-reconstructions.tiff", 
       fig.sg.traj.other, width = 7.5, height = 10 * 3/14, units = "in")


# Invasives --------------------------------------------------------------------

dat.fig.inv <- subset(dat.fig, sample_type == "Invasive")

# Change tissue names to be shorter
dat.fig.inv$sample_rep <- str_replace_all(dat.fig.inv$sample_rep, "right upper lobe", "RUL")
dat.fig.inv$sample_rep <- str_replace_all(dat.fig.inv$sample_rep, "right lower lobe", "RLL")
dat.fig.inv$sample_rep <- str_replace_all(dat.fig.inv$sample_rep, "right middle lobe", "RML")
dat.fig.inv$sample_rep <- str_replace_all(dat.fig.inv$sample_rep, "right middle", "RM")
dat.fig.inv$sample_rep <- str_replace_all(dat.fig.inv$sample_rep, "right lower", "RL")
dat.fig.inv$sample_rep <- str_replace_all(dat.fig.inv$sample_rep, "right accessory lobe", "RAL")
dat.fig.inv$sample_rep <- str_replace_all(dat.fig.inv$sample_rep, "left upper lobe", "LUL")
dat.fig.inv$sample_rep <- str_replace_all(dat.fig.inv$sample_rep, "left lower lobe", "LLL")
dat.fig.inv$sample_rep <- str_replace_all(dat.fig.inv$sample_rep, "left middle lobe", "LML")
dat.fig.inv$sample_rep <- str_replace_all(dat.fig.inv$sample_rep, "upper lobe", "UL")
dat.fig.inv$sample_rep <- str_replace_all(dat.fig.inv$sample_rep, "lower lobe", "LL")
dat.fig.inv$sample_rep <- str_replace_all(dat.fig.inv$sample_rep, "average of 7 lobes", "avg.")
dat.fig.inv$sample_rep <- str_replace_all(dat.fig.inv$sample_rep, "left", "L")
dat.fig.inv$sample_rep <- str_replace_all(dat.fig.inv$sample_rep, "right", "R")
dat.fig.inv$sample_rep <- str_replace_all(dat.fig.inv$sample_rep, "Descending colon", "Colon")
dat.fig.inv$sample_rep <- str_replace_all(dat.fig.inv$sample_rep, "Intestine (large)", "L. intestine")
dat.fig.inv$sample_rep <- str_replace_all(dat.fig.inv$sample_rep, "Small intestine", "S. intestine")


# Indicator 
dat.fig.inv$organ_combined <- dat.fig.inv$organ_system
dat.fig.inv$organ_combined[dat.fig.inv$organ_combined %in% c("GI", "Other")] <- "GI & Other"

dat.fig.inv$system_sg_pred <- paste0(dat.fig.inv$organ_combined, ": sgRNA (predicted)")
dat.fig.inv$system_sg <- paste0(dat.fig.inv$organ_combined, ": sgRNA (observed)")
dat.fig.inv$system_total <- paste0(dat.fig.inv$organ_combined, ": total RNA (observed)")
dat.fig.inv$system_culture <- paste0(dat.fig.inv$organ_combined, "culture")
#
## Order samples by organ system
tissue_order <- c(sort(unique(dat.fig.inv$sample_rep[dat.fig.inv$organ_system == "URT"])),
                  sort(unique(dat.fig.inv$sample_rep[dat.fig.inv$organ_system == "LRT"])),
                  sort(unique(dat.fig.inv$sample_rep[dat.fig.inv$organ_system == "GI"])),
                  sort(unique(dat.fig.inv$sample_rep[dat.fig.inv$organ_system == "Other"])))

dat.fig.inv$sample_factor <- factor(dat.fig.inv$sample_rep,
                                    levels = tissue_order,
                                    ordered = TRUE)

dat.fig.inv$panel_label <- paste0(dat.fig.inv$dpi, ": ", dat.fig.inv$indiv)

dat.fig.inv$label_location <- NA
for (indiv.ii in unique(dat.fig.inv$indiv)) {
  dat.indiv <- subset(dat.fig.inv, indiv == indiv.ii)
  dat.indiv <- arrange(dat.indiv, match(dat.indiv$sample_rep, levels(dat.indiv$sample_factor)))
  first.tissue <- dat.indiv$sample_rep[1]
  dat.fig.inv$label_location[dat.fig.inv$indiv == indiv.ii] <- first.tissue
}

dat.fig.inv <- arrange(dat.fig.inv, match(dat.fig.inv$sample_rep, levels(dat.fig.inv$sample_factor)))

annot.df <- dat.fig.inv
annot.df$sample_factor <- dat.fig.inv$label_location


# Invasive Samples
fig.sg.inv <- ggplot(dat.fig.inv) + 
  
  geom_hline(data = subset(dat.fig.inv, lodq_total == -9),
             aes(yintercept = 0), linetype = "dotted",
             size = 0.5, color = "grey22", alpha = 0.7) +
  geom_hline(data = subset(dat.fig.inv, lodq_total != -9),
             aes(yintercept = lodq_total, color = sample_total), linetype = "dotted",
             size = 0.5, alpha = 0.7) +
  geom_hline(data = subset(dat.fig.inv, lodq_sg != -9),
             aes(yintercept = lodq_sg, color = sample_sg), linetype = "dotted",
             size = 0.5, alpha = 0.7) +
  
  
  geom_point(aes(x = sample_factor, 
                 y = val_sg_pred, 
                 group = indiv, 
                 fill = system_sg_pred),
             alpha = 1, size = 1.2, shape = 24) +
  geom_point(aes(x = sample_factor,
                 y = val_sg, 
                 group = indiv, 
                 fill = system_sg),
             alpha = 1, size = 1.2, shape = 23) +
  geom_point(aes(x = sample_factor, 
                 y = val_total, 
                 group = indiv_sample,
                 fill = system_total),
             alpha = 1, size = 1.2, shape = 21) +
  scale_fill_manual(values = c("URT: total RNA (observed)" = "#0A2472",
                               "URT: sgRNA (observed)" = zissou_pal[1],
                               "URT: sgRNA (predicted)" = "#0A8549",
                               "LRT: total RNA (observed)" = zissou_pal[12],
                               "LRT: sgRNA (observed)" = zissou_pal[10],
                               "LRT: sgRNA (predicted)" = zissou_pal[7],
                               "GI & Other: total RNA (observed)" = "#710193",
                               "GI & Other: sgRNA (observed)" = "#C36F8E",
                               "GI & Other: sgRNA (predicted)" = "light pink"),
                    breaks = c("URT: total RNA (observed)", "URT: sgRNA (observed)", "URT: sgRNA (predicted)",
                               "LRT: total RNA (observed)", "LRT: sgRNA (observed)", "LRT: sgRNA (predicted)",
                               "GI & Other: total RNA (observed)", "GI & Other: sgRNA (observed)",
                               "GI & Other: sgRNA (predicted)")) +
  scale_color_manual(values = c("URT: total RNA (observed)" = "#0A2472",
                                "URT: sgRNA (observed)" = zissou_pal[1],
                                "URT: sgRNA (predicted)" = "#0A8549",
                                "LRT: total RNA (observed)" = zissou_pal[12],
                                "LRT: sgRNA (observed)" = zissou_pal[10],
                                "LRT: sgRNA (predicted)" = zissou_pal[7],
                                "GI & Other: total RNA (observed)" = "#710193",
                                "GI & Other: sgRNA (observed)" = "#C36F8E",
                                "GI & Other: sgRNA (predicted)" = "light pink")) +
  facet_wrap(.~ dpi + indiv, 
             labeller = function (labels) {
               labels <- lapply(labels, as.character)
               list(do.call(paste, c(labels, list(sep = ": "))))
             },
             scales = "free_x", 
             ncol = 8) +
  labs(y = element_blank(), x = element_blank(),
       fill = "System: Assay") +
  theme(
    legend.position = "bottom", 
    legend.justification = "right",
    legend.box.spacing = unit(0, "pt"),
    legend.spacing.y = unit(0.2, 'cm'),
    legend.margin = margin(-4, 0, -4, -4),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.background = element_rect(fill='white'),
    legend.key=element_blank(),
    legend.key.size = unit(1, "lines"),
    axis.title.x = element_text(vjust = - 1),
    axis.text.x = element_text(size = 6, angle = 65, hjust = 1),
    axis.text.y = element_text(size = 7),
    axis.ticks.length.x = unit(1, "pt"),
    text = element_text(size = 9),
    strip.background = element_rect(colour = "black", fill = "white"),
    strip.text.x = element_text(size = 6.5, face = "bold",
                                margin = margin(.1, 0, .1, 0, "cm")),
    panel.spacing = unit(0.1, "lines"), 
    panel.background = element_rect(fill = "white",
                                    colour = "black",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                    colour = "light grey"), 
    panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                    colour = "light grey"),
    plot.margin = margin(10, 10, 10, 5)) + 
  guides(color = "none", 
         fill = guide_legend(nrow = 3, byrow = FALSE,
                             override.aes = list(size = 1.6,
                             shape = c(21, 23, 24, 21, 23, 24, 21, 23, 24)))) +
  scale_y_continuous(breaks = seq(from = 0, to = 10, by = 2), 
                     minor_breaks = seq(1, 9, 2),
                     limits = c(-1.5, 10.2),
                     expand = c(0, 0)); fig.sg.inv

ggsave("./outputs/figures/figS6-sgRNA-invasive-reconstructions.tiff", 
       fig.sg.inv, width = 7.5, height = 10, units = "in")
