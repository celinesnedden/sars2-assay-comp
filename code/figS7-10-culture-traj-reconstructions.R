# Generates supplemental figures 7 through 10
#   Plots all individual trajectories with culture data, stratified by organ system
# No dependencies, this file will run without running any other scripts

# Prep environment -------------------------------------------------------------

# Install & load color palette & ggplot 
req_pkgs <- c("wesanderson", "ggplot2", "ggforce", "stringr",
              "patchwork", "gridExtra")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependences = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))

# Extract desired color palette
zissou_pal <- wes_palette("Zissou1", 12, type = "continuous")

# Minor function for convenience
`%notin%` <- Negate(`%in%`)



# Prep the data ----------------------------------------------------------------

dat.fig <- read.csv("./data/pred-culture-data.csv")
dat.fig <- subset(dat.fig, article != "Johnston et al. 2020")

# Convert columns to numerics
dat.fig$val_sg <- as.numeric(dat.fig$val_sg)
dat.fig$val_total <- as.numeric(dat.fig$val_total)

# Combine indiv names with sample rep to get the right connecting lines
dat.fig$indiv_sample <- paste0(dat.fig$indiv, "_", dat.fig$sample_rep, 
                               "_", dat.fig$tg_idx)

# Get all individuals with at least two inf virus samples 
indiv_list <- c()
for (indiv.ii in unique(dat.fig$indiv_sample)) {
  dat.indiv <- subset(dat.fig, indiv_sample == indiv.ii)
  dat.indiv.inf <- subset(dat.indiv, !is.na(pos_inf))
  if (nrow(dat.indiv.inf) >= 2) {
    indiv_list <- c(indiv_list, indiv.ii)
  }
}
dat.fig <- subset(dat.fig, indiv_sample %in% indiv_list)


# Set RNA values for samples below LOD 
dat.fig$val_sg[dat.fig$cens_sg == 1] <- -1
dat.fig$val_total[dat.fig$cens_total == 1] <- 0
dat.fig$val_sg_pred[dat.fig$pos_sg_pred == 0 | dat.fig$val_sg_pred == -9] <- -1
dat.fig$val_sg[dat.fig$val_sg == -9] <- -1


# Indicator 
dat.fig$sample_sg <- paste0("sgRNA: ", dat.fig$organ_system)
dat.fig$sample_total <- paste0("total RNA: ", dat.fig$organ_system)
dat.fig$sample_culture <- paste0("culture: ", dat.fig$organ_system)


# Annotate names with symbols to correlate between figures
dat.fig$indiv_annot <- NA
for (indiv.ii in unique(dat.fig$indiv)){
  indiv.sub <- subset(dat.fig, indiv == indiv.ii & sample_type == "Non-invasive")
  samples <- unique(indiv.sub$organ_system)
  
  indiv_annot <- paste0(indiv.ii, " ")
  if ("URT" %in% samples){indiv_annot <- paste0(indiv_annot, "*")}
  if ("LRT" %in% samples){indiv_annot <- paste0(indiv_annot, "#")}
  if ("GI" %in% samples){indiv_annot <- paste0(indiv_annot, "^")}
  if ("Other" %in% samples){indiv_annot <- paste0(indiv_annot, "&")}
  
  dat.fig$indiv_annot[dat.fig$indiv == indiv.ii] <- indiv_annot
}


# For ease of visualization
dat.fig$organ_system <- factor(dat.fig$organ_system, levels = c("URT", "LRT", "GI", "Other"))
dat.fig$sample_rep[dat.fig$sample_rep == "Nasopharyngeal Swab"] <- "Nasoph. Swab"
dat.fig$sample_rep[dat.fig$sample_rep == "Oropharyngeal Swab"] <- "Oroph. Swab"

dat.fig$sample_rep[dat.fig$sample_rep == "Nasal Cavity Swab"] <- "Nasal Swab"

dat.fig$panel_label <-paste0(dat.fig$indiv, "\n", dat.fig$sample_rep)

# Non-invasive samples ---------------------------------------------------------

dat.fig.ni <- subset(dat.fig, sample_type == "Non-invasive")


### S7: URT only ---------------------------------------------------------------

dat.fig.ni.urt <- subset(dat.fig.ni, organ_system == "URT")

figS8 <- ggplot(dat.fig.ni.urt) + 
  geom_text(aes(label = panel_label, x = 10.5, y = 10.5), size = 1.8) +
  
  geom_hline(data = subset(dat.fig.ni.urt, lodq_total == -9),
             aes(yintercept = 0), linetype = "dotted",
             size = 0.5, color = "grey22", alpha = 0.7) +
  geom_hline(data = subset(dat.fig.ni.urt, lodq_total != -9),
             aes(yintercept = lodq_total), linetype = "dotted",
             size = 0.5, color = "#0A2472", alpha = 0.7) +
  geom_hline(data = subset(dat.fig.ni.urt, lodq_total != -9),
             aes(yintercept = lodq_sg), linetype = "dotted",
             size = 0.5, color = zissou_pal[1], alpha = 0.7) +
  
  geom_hline(aes(yintercept = 12.5), color = zissou_pal[8], size = 0.5) +

  geom_point(data = subset(dat.fig.ni.urt, pos_inf == 1 & pos_inf_pred == 1),
             aes(x = dpi, y = 13.5, 
                 fill = "Culture (True Positive)"),
             color = "black", 
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.fig.ni.urt, pos_inf == 1 & pos_inf_pred == 0),
             aes(x = dpi, y = 13.5, 
                 fill = "Culture (False Negative)"),
             color = "transparent",
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.fig.ni.urt, pos_inf == 1 & is.na(pos_inf_pred)),
             aes(x = dpi, y = 13.5, 
                 fill = "Culture (True Positive)"),
             color = zissou_pal[1], 
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.fig.ni.urt, pos_inf == 0 & is.na(pos_inf_pred)),
             aes(x = dpi, y = 13.5, 
                 fill = "Culture (True Negative)"),
             color = zissou_pal[1], 
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.fig.ni.urt, pos_inf == 1 & pos_inf_pred == 0),
             aes(x = dpi, y = 13.5, 
                 fill = "Culture (False Negative)"),
             color = "transparent",
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.fig.ni.urt, pos_inf == 0 & pos_inf_pred == 0),
             aes(x = dpi, y = 13.5, fill = "Culture (True Negative)"),  color = "black",
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.fig.ni.urt, pos_inf == 0 & pos_inf_pred == 1),
             aes(x = dpi, y = 13.5, fill = "Culture (False Positive)"),  color = "transparent",
             size = 1.5, shape = 22) +
  geom_line(aes(x = dpi, y = val_sg, group = indiv_sample,
                color = "sgRNA (observed)"),
            alpha = 0.8) +
  geom_point(aes(x = dpi, y = val_sg, 
                 group = indiv_sample, 
                 fill = "sgRNA (observed)"),
             alpha = 1, size = 1.5, shape = 23) +
  geom_line(data = subset(dat.fig.ni.urt, is.na(val_sg)),
            aes(x = dpi, y = val_sg_pred, group = indiv_sample,
                color = "sgRNA (predicted)"),
            alpha = 0.8) +
  geom_point(data = subset(dat.fig.ni.urt, is.na(val_sg)),
             aes(x = dpi, y = val_sg_pred, 
                 fill = "sgRNA (predicted)",
                 group = indiv_sample), 
             alpha = 0.8, size = 1.5, shape = 24) +
  geom_line(aes(x = dpi, y = val_total, group = indiv_sample,
                color = "total RNA (observed)"),
            alpha = 0.8) +
  geom_point(aes(x = dpi, y = val_total, 
                 group = indiv_sample,
                 fill = "total RNA (observed)"),
             alpha = 1, size = 1.5, shape = 21) +
  scale_fill_manual(values = c("total RNA (observed)" = "#0A2472",
                               "sgRNA (observed)" = zissou_pal[1],
                               "sgRNA (predicted)" = "#0A8549",
                               "Culture (True Positive)" = zissou_pal[8],
                               "Culture (True Negative)" = "grey62",
                               "Culture (False Positive)" = "grey62",
                               "Culture (False Negative)" = zissou_pal[8]),
                    breaks = c("Culture (True Positive)",
                               "Culture (True Negative)", "Culture (False Positive)",
                               "Culture (False Negative)",
                               "total RNA (observed)", "sgRNA (observed)",
                               "sgRNA (predicted)")) +
  scale_color_manual(values = c("total RNA (observed)" = "#0A2472",
                                "sgRNA (observed)" = zissou_pal[1],
                                "sgRNA (predicted)" = "#0A8549")) +
  facet_wrap(str_remove(indiv_annot, "[*]") + sample_rep ~., ncol = 9) +
  labs(y = "log10 RNA copies / sample", x = "day post infection",
       fill = "Assay") +
  theme(#legend.position = c(0.94, 0.5),
    legend.position = "bottom", 
    legend.justification = "right",
    legend.box.spacing = unit(0, "pt"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.background = element_rect(fill='white'),
    legend.key = element_blank(),
    legend.key.size = unit(1, "lines"),
    axis.text = element_text(size = 7),
    text = element_text(size = 9),
    strip.background = element_rect(colour = "black", fill = "white"),
    strip.text = element_blank(),
    panel.spacing = unit(0.1, "lines"), 
    panel.background = element_rect(fill = "white",
                                    colour = "black",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                    colour = "light grey"), 
    panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                    colour = "light grey"),
    plot.margin = margin(10, 10, 10, 5)) + 
  guides(color = "none", #, 
         fill = guide_legend(override.aes = list(size = 1.5,
                                                 shape = c( 22, 22, 22, 22, 21, 23, 24),
                                                 color = c("black", "black", "white", "white",
                                                           "black", "black", "black")),
                             nrow = 2)) +
  scale_x_continuous(breaks = seq(0, 16, by = 4), limits = c(0, 16),
                     labels = c(seq(0, 12, by = 4), ""), 
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(from = 0, to = 12, by = 2), 
                     minor_breaks = seq(1, 11, 2),
                     limits = c(-2, 14.5),
                     expand = c(0, 0)); figS8


ggsave("./outputs/figures/figS7-culture-URT-traj-reconstructions.tiff", figS8,
       width = 7.5, height = 10, units = "in")


### S8: LRT only -------------------------------------------------------------------

dat.fig.ni.lrt <- subset(dat.fig.ni, organ_system == "LRT")

figS9 <- ggplot(dat.fig.ni.lrt) + 
  geom_text(aes(label = panel_label, x = 10.5, y = 10.5), size = 1.8) +
  
  geom_hline(data = subset(dat.fig.ni.lrt, lodq_total == -9),
             aes(yintercept = 0), linetype = "dotted",
             size = 0.5, color = "grey22", alpha = 0.7) +
  geom_hline(data = subset(dat.fig.ni.lrt, lodq_total != -9),
             aes(yintercept = lodq_total), linetype = "dotted",
             size = 0.5, color = "#0A2472", alpha = 0.7) +
  geom_hline(data = subset(dat.fig.ni.lrt, lodq_total != -9),
             aes(yintercept = lodq_sg), linetype = "dotted",
             size = 0.5, color = zissou_pal[1], alpha = 0.7) +
  
  geom_hline(aes(yintercept = 12.5), color = zissou_pal[8], size = 0.5) +
  
  geom_point(data = subset(dat.fig.ni.lrt, pos_inf == 1 & pos_inf_pred == 1),
             aes(x = dpi, y = 13.5, 
                 fill = "Culture (True Positive)"),
             color = "black", 
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.fig.ni.lrt, pos_inf == 1 & pos_inf_pred == 0),
             aes(x = dpi, y = 13.5, 
                 fill = "Culture (False Negative)"),
             color = "transparent",
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.fig.ni.lrt, pos_inf == 1 & is.na(pos_inf_pred)),
             aes(x = dpi, y = 13.5, 
                 fill = "Culture (True Positive)"),
             color = zissou_pal[1], 
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.fig.ni.lrt, pos_inf == 0 & is.na(pos_inf_pred)),
             aes(x = dpi, y = 13.5, 
                 fill = "Culture (True Negative)"),
             color = zissou_pal[1], 
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.fig.ni.lrt, pos_inf == 1 & pos_inf_pred == 0),
             aes(x = dpi, y = 13.5, 
                 fill = "Culture (False Negative)"),
             color = "transparent",
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.fig.ni.lrt, pos_inf == 0 & pos_inf_pred == 0),
             aes(x = dpi, y = 13.5, fill = "Culture (True Negative)"),  color = "black",
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.fig.ni.lrt, pos_inf == 0 & pos_inf_pred == 1),
             aes(x = dpi, y = 13.5, fill = "Culture (False Positive)"),  color = "transparent",
             size = 1.5, shape = 22) +
  geom_line(aes(x = dpi, y = val_sg, group = indiv_sample,
                color = "sgRNA (observed)"),
            alpha = 0.8) +
  geom_point(aes(x = dpi, y = val_sg, 
                 group = indiv_sample, 
                 fill = "sgRNA (observed)"),
             alpha = 1, size = 1.5, shape = 23) +
  geom_line(data = subset(dat.fig.ni.lrt, is.na(val_sg)),
            aes(x = dpi, y = val_sg_pred, group = indiv_sample,
                color = "sgRNA (predicted)"),
            alpha = 0.8) +
  geom_point(data = subset(dat.fig.ni.lrt, is.na(val_sg)),
             aes(x = dpi, y = val_sg_pred, 
                 fill = "sgRNA (predicted)",
                 group = indiv_sample), 
             alpha = 0.8, size = 1.5, shape = 24) +
  geom_line(aes(x = dpi, y = val_total, group = indiv_sample,
                color = "total RNA (observed)"),
            alpha = 0.8) +
  geom_point(aes(x = dpi, y = val_total, 
                 group = indiv_sample,
                 fill = "total RNA (observed)"),
             alpha = 1, size = 1.5, shape = 21) +
  scale_fill_manual(values = c("total RNA (observed)" = zissou_pal[12],
                               "sgRNA (observed)" = zissou_pal[10],
                               "sgRNA (predicted)" = zissou_pal[7],
                               "Culture (True Positive)" = zissou_pal[8],
                               "Culture (True Negative)" = "grey62",
                               "Culture (False Positive)" = "grey62",
                               "Culture (False Negative)" = zissou_pal[8]),
                    breaks = c("Culture (True Positive)",
                               "Culture (True Negative)", "Culture (False Positive)",
                               "Culture (False Negative)",
                               "total RNA (observed)", "sgRNA (observed)",
                               "sgRNA (predicted)")) +
  scale_color_manual(values = c("total RNA (observed)" = zissou_pal[12],
                                "sgRNA (observed)" = zissou_pal[10],
                                "sgRNA (predicted)" = zissou_pal[7])) +
  facet_wrap(str_remove(indiv_annot, "[*]") + sample_rep ~., ncol = 9) +
  labs(y = "log10 RNA copies / sample", x = "day post infection",
       fill = "Assay") +
  theme(#legend.position = c(0.94, 0.5),
    legend.position = "bottom", 
    legend.justification = "right",
    legend.box.spacing = unit(0, "pt"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.background = element_rect(fill='white'),
    legend.key=element_blank(),
    legend.key.size = unit(1, "lines"),
    axis.text = element_text(size = 7),
    text = element_text(size = 9),
    strip.background = element_rect(colour = "black", fill = "white"),
    strip.text = element_blank(),
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
         fill = guide_legend(override.aes = list(size = 1.5,
                                                 shape = c( 22, 22, 22, 22, 21, 23, 24),
                                                 color = c("black", "black", "white", "white",
                                                           "black", "black", "black")),
                             nrow = 2)) +
  scale_x_continuous(breaks = seq(0, 16, by = 4), limits = c(0, 16),
                     labels = c(seq(0, 12, by = 4), ""), 
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(from = 0, to = 12, by = 2), 
                     minor_breaks = seq(1, 11, 2),
                     limits = c(-2, 14.5),
                     expand = c(0, 0)); figS9


ggsave("./outputs/figures/figS8-culture-LRT-traj-reconstructions.tiff", figS9,
       width = 7.5, height = 10 * 6/12, units = "in")


### S9: GI & other -----------------------------------------------------------------

dat.fig.ni.other <- subset(dat.fig.ni, organ_system %in% c("GI", "Other"))

figS10 <- ggplot(dat.fig.ni.other) + 
  geom_text(aes(label = panel_label, x = 10.5, y = 10.5), size = 1.8) +
  
  geom_hline(data = subset(dat.fig.ni.other, lodq_total == -9),
             aes(yintercept = 0), linetype = "dotted",
             size = 0.5, color = "grey22", alpha = 0.7) +
  geom_hline(data = subset(dat.fig.ni.other, lodq_total != -9),
             aes(yintercept = lodq_total), linetype = "dotted",
             size = 0.5, color = "#0A2472", alpha = 0.7) +
  geom_hline(data = subset(dat.fig.ni.other, lodq_total != -9),
             aes(yintercept = lodq_sg), linetype = "dotted",
             size = 0.5, color = zissou_pal[1], alpha = 0.7) +
  
  geom_hline(aes(yintercept = 12.5), color = zissou_pal[8], size = 0.5) +
  
  geom_point(data = subset(dat.fig.ni.other, pos_inf == 1 & pos_inf_pred == 1),
             aes(x = dpi, y = 13.5, 
                 fill = "Culture (True Positive)"),
             color = "black", 
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.fig.ni.other, pos_inf == 1 & pos_inf_pred == 0),
             aes(x = dpi, y = 13.5, 
                 fill = "Culture (False Negative)"),
             color = "transparent",
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.fig.ni.other, pos_inf == 1 & is.na(pos_inf_pred)),
             aes(x = dpi, y = 13.5, 
                 fill = "Culture (True Positive)"),
             color = "#710193", 
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.fig.ni.other, pos_inf == 0 & is.na(pos_inf_pred)),
             aes(x = dpi, y = 13.5, 
                 fill = "Culture (True Negative)"),
             color = "#710193", 
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.fig.ni.other, pos_inf == 1 & pos_inf_pred == 0),
             aes(x = dpi, y = 13.5, 
                 fill = "Culture (False Negative)"),
             color = "transparent",
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.fig.ni.other, pos_inf == 0 & pos_inf_pred == 0),
             aes(x = dpi, y = 13.5, fill = "Culture (True Negative)"),  color = "black",
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.fig.ni.other, pos_inf == 0 & pos_inf_pred == 1),
             aes(x = dpi, y = 13.5, fill = "Culture (False Positive)"),  color = "transparent",
             size = 1.5, shape = 22) +
  geom_line(aes(x = dpi, y = val_sg, group = indiv_sample,
                color = "sgRNA (observed)"),
            alpha = 0.8) +
  geom_point(aes(x = dpi, y = val_sg, 
                 group = indiv_sample, 
                 fill = "sgRNA (observed)"),
             alpha = 1, size = 1.5, shape = 23) +
  geom_line(data = subset(dat.fig.ni.other, is.na(val_sg)),
            aes(x = dpi, y = val_sg_pred, group = indiv_sample,
                color = "sgRNA (predicted)"),
            alpha = 0.8) +
  geom_point(data = subset(dat.fig.ni.other, is.na(val_sg)),
             aes(x = dpi, y = val_sg_pred, 
                 fill = "sgRNA (predicted)",
                 group = indiv_sample), 
             alpha = 0.8, size = 1.5, shape = 24) +
  geom_line(aes(x = dpi, y = val_total, group = indiv_sample,
                color = "total RNA (observed)"),
            alpha = 0.8) +
  geom_point(aes(x = dpi, y = val_total, 
                 group = indiv_sample,
                 fill = "total RNA (observed)"),
             alpha = 1, size = 1.5, shape = 21) +
  scale_fill_manual(values = c("total RNA (observed)" = "#710193",
                               "sgRNA (observed)" = "#C36F8E",
                               "sgRNA (predicted)" = "light pink",
                               "Culture (True Positive)" = zissou_pal[8],
                               "Culture (True Negative)" = "grey62",
                               "Culture (False Positive)" = "grey62",
                               "Culture (False Negative)" = zissou_pal[8]),
                    breaks = c("Culture (True Positive)",
                               "Culture (True Negative)", "Culture (False Positive)",
                               "Culture (False Negative)",
                               "total RNA (observed)", "sgRNA (observed)",
                               "sgRNA (predicted)")) +
  scale_color_manual(values = c("total RNA (observed)" = "#710193",
                                "sgRNA (observed)" = "#C36F8E",
                                "sgRNA (predicted)" = "light pink")) +
  facet_wrap(str_remove(indiv_annot, "[*]") + sample_rep ~., ncol = 9) +
  labs(y = "log10 RNA copies / sample", x = "day post infection",
       fill = "Assay") +
  theme(
    legend.position = "bottom", 
    legend.justification = "right",
    legend.box.spacing = unit(0, "pt"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.background = element_rect(fill='white'),
    legend.key=element_blank(),
    legend.key.size = unit(1, "lines"),
    axis.text = element_text(size = 7),
    text = element_text(size = 9),
    strip.background = element_rect(colour = "black", fill = "white"),
    strip.text = element_blank(),
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
         fill = guide_legend(override.aes = list(size = 1.5,
                                                 shape = c( 22, 22, 22, 22, 21, 23, 24),
                                                 color = c("black", "black", "white", "white",
                                                           "black", "black", "black")),
                             nrow = 2)) +
  scale_x_continuous(breaks = seq(0, 16, by = 4), limits = c(0, 16),
                     labels = c(seq(0, 12, by = 4), ""), 
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(from = 0, to = 12, by = 2), 
                     minor_breaks = seq(1, 11, 2),
                     limits = c(-2, 14.5),
                     expand = c(0, 0)); figS10


ggsave("./outputs/figures/figS9-culture-Other-traj-reconstructions.tiff", figS10,
       width = 7.5, height = 10 * 9/12, units = "in")


# S10: Invasive samples --------------------------------------------------------

dat.fig <- read.csv("./data/pred-culture-data.csv")
dat.fig <- subset(dat.fig, article != "Johnston et al. 2020")

# Removes data without culture samples
dat.inv <- subset(dat.fig, !is.na(pos_inf) & sample_type == "Invasive") 

# Convert columns to numerics
dat.inv$val_sg <- as.numeric(dat.inv$val_sg)
dat.inv$val_total <- as.numeric(dat.inv$val_total)

# Set RNA values for samples below LOD 
dat.inv$val_sg[dat.inv$cens_sg == 1] <- -1
dat.inv$val_total[dat.inv$cens_total == 1] <- 0
dat.inv$val_sg_pred[dat.inv$pos_sg_pred == 0 | dat.inv$val_sg_pred == -9] <- -1
dat.inv$val_sg[dat.inv$val_sg == -9] <- -1

# Change tissue names to be shorter
dat.inv$sample_rep <- str_replace_all(dat.inv$sample_rep, "right upper lobe", "RUL")
dat.inv$sample_rep <- str_replace_all(dat.inv$sample_rep, "right lower lobe", "RLL")
dat.inv$sample_rep <- str_replace_all(dat.inv$sample_rep, "right middle lobe", "RML")
dat.inv$sample_rep <- str_replace_all(dat.inv$sample_rep, "right middle", "RM")
dat.inv$sample_rep <- str_replace_all(dat.inv$sample_rep, "right lower", "RL")
dat.inv$sample_rep <- str_replace_all(dat.inv$sample_rep, "right accessory lobe", "RAL")
dat.inv$sample_rep <- str_replace_all(dat.inv$sample_rep, "left upper lobe", "LUL")
dat.inv$sample_rep <- str_replace_all(dat.inv$sample_rep, "left lower lobe", "LLL")
dat.inv$sample_rep <- str_replace_all(dat.inv$sample_rep, "left middle lobe", "LML")
dat.inv$sample_rep <- str_replace_all(dat.inv$sample_rep, "upper lobe", "UL")
dat.inv$sample_rep <- str_replace_all(dat.inv$sample_rep, "lower lobe", "LL")
dat.inv$sample_rep <- str_replace_all(dat.inv$sample_rep, "average of 7 lobes", "avg.")
dat.inv$sample_rep <- str_replace_all(dat.inv$sample_rep, "left", "L")
dat.inv$sample_rep <- str_replace_all(dat.inv$sample_rep, "right", "R")
dat.inv$sample_rep <- str_replace_all(dat.inv$sample_rep, "Descending colon", "Colon")
dat.inv$sample_rep <- str_replace_all(dat.inv$sample_rep, "Intestine (large)", "L. intestine")
dat.inv$sample_rep <- str_replace_all(dat.inv$sample_rep, "Small intestine", "S. intestine")

dat.inv$id_tissue <- paste0(dat.inv$indiv, ": ", dat.inv$sample_rep, 
                            "_", dat.inv$tg_idx)

dat.inv$organ_combined <- dat.inv$organ_system
dat.inv$organ_combined[dat.inv$organ_combined %in% c("GI", "Other")] <- "GI & Other"

dat.inv$system_sg_pred <- paste0(dat.inv$organ_combined, ": sgRNA (predicted)")
dat.inv$system_sg <- paste0(dat.inv$organ_combined, ": sgRNA (observed)")
dat.inv$system_total <- paste0(dat.inv$organ_combined, ": total RNA (observed)")
dat.inv$system_culture <- paste0(dat.inv$organ_combined, "culture")
#
## Order samples by organ system
tissue_order <- c(sort(unique(dat.inv$sample_rep[dat.inv$organ_system == "URT"])),
                  sort(unique(dat.inv$sample_rep[dat.inv$organ_system == "LRT"])),
                  sort(unique(dat.inv$sample_rep[dat.inv$organ_system == "GI"])),
                  sort(unique(dat.inv$sample_rep[dat.inv$organ_system == "Other"])))

dat.inv$sample_factor <- factor(dat.inv$sample_rep,
                                    levels = tissue_order,
                                    ordered = TRUE)

dat.inv$panel_label <- paste0(dat.inv$dpi, ": ", dat.inv$indiv)

dat.inv$label_location <- NA
for (indiv.ii in unique(dat.inv$indiv)) {
  dat.indiv <- subset(dat.inv, indiv == indiv.ii)
  dat.indiv <- arrange(dat.indiv, match(dat.indiv$sample_rep, levels(dat.indiv$sample_factor)))
  first.tissue <- dat.indiv$sample_rep[2]
  dat.inv$label_location[dat.inv$indiv == indiv.ii] <- first.tissue
}

dat.inv <- arrange(dat.inv, match(dat.inv$sample_rep, levels(dat.inv$sample_factor)))

# Invasive Samples
figS11 <- ggplot(dat.inv) + 
  geom_hline(aes(yintercept = 10), color = zissou_pal[8], size = 0.5) +
  
  geom_hline(data = subset(dat.inv, lodq_total == -9),
             aes(yintercept = 0), linetype = "dotted",
             size = 0.5, color = "grey22", alpha = 0.7) +
  geom_hline(data = subset(dat.inv, lodq_total != -9),
             aes(yintercept = lodq_total, color = system_total), linetype = "dotted",
             size = 0.5, alpha = 0.7) +
  geom_hline(data = subset(dat.inv, lodq_sg != -9),
             aes(yintercept = lodq_sg, color = system_sg), linetype = "dotted",
             size = 0.5, alpha = 0.7) +
  
  geom_point(data = subset(dat.inv, pos_inf == 1 & pos_inf_pred == 1),
             aes(x = sample_factor, y = 11, 
                 fill = "Culture (True Positive)"),
             color = "black", 
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.inv, pos_inf == 1 & pos_inf_pred == 0),
             aes(x = sample_factor, y = 11, 
                 fill = "Culture (False Negative)"),
             color = "transparent",
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.inv, pos_inf == 1 & is.na(pos_inf_pred)),
             aes(x = sample_factor, y = 11, 
                 fill = "Culture (True Positive)"),
             color = zissou_pal[1], 
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.inv, pos_inf == 0 & is.na(pos_inf_pred)),
             aes(x = sample_factor, y = 11, 
                 fill = "Culture (True Negative)"),
             color = zissou_pal[1], 
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.inv, pos_inf == 1 & pos_inf_pred == 0),
             aes(x = sample_factor, y = 11, 
                 fill = "Culture (False Negative)"),
             color = "transparent",
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.inv, pos_inf == 0 & pos_inf_pred == 0),
             aes(x = sample_factor, y = 11, fill = "Culture (True Negative)"),  
             color = "black",
             size = 1.5, shape = 22) +
  geom_point(data = subset(dat.inv, pos_inf == 0 & pos_inf_pred == 1),
             aes(x = sample_factor, y = 11, fill = "Culture (False Positive)"),  
             color = "transparent",
             size = 1.5, shape = 22) +
  
  geom_point(data = subset(dat.inv, is.na(val_sg)),
             aes(x = sample_factor, 
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
                 group = indiv,
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
                               "GI & Other: sgRNA (predicted)" = "light pink",
                               "Culture (True Positive)" = zissou_pal[8],
                               "Culture (True Negative)" = "grey62",
                               "Culture (False Positive)" = "grey62",
                               "Culture (False Negative)" = zissou_pal[8]),
                    breaks = c("URT: total RNA (observed)", "URT: sgRNA (observed)",
                               "URT: sgRNA (predicted)",
                               "LRT: total RNA (observed)", "LRT: sgRNA (observed)",
                               "LRT: sgRNA (predicted)",
                               "GI & Other: total RNA (observed)", 
                               "GI & Other: sgRNA (observed)",
                               "GI & Other: sgRNA (predicted)",
                               "Culture (True Positive)",
                               "Culture (True Negative)",
                               "Culture (False Positive)",
                               "Culture (False Negative)")) +
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
       fill = "System:Assay") +
  theme(
    legend.position = "bottom", 
    legend.justification = "right",
    legend.box.spacing = unit(0, "pt"),
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
         fill = guide_legend(nrow = 4, byrow = FALSE,
                             override.aes = list(shape = c(21, 23, 21, 23, 24, 21, 23, 24, 22, 22, 22, 22),
                                                 size = 1.5,
                                                 color = c(rep("black", 10), "white", "white")))) +
  scale_y_continuous(breaks = seq(from = 0, to = 10, by = 2), 
                     minor_breaks = seq(1, 9, 2),
                     labels = c(0, 2, 4, 6, 8, 10),
                     limits = c(-1.5, 12),
                     expand = c(0, 0)); figS11

ggsave("./outputs/figures/figS10-culture-invasive-reconstructions.tiff", figS11,
       width = 7.5, height = 9, units = "in")
