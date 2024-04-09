# Generates figure 4, some example reconstructions of sgRNA trajectories
#   See the supplemental figures for all reconstructed trajectories


# Load & prep data ------------------------------------------------------------

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


# URT samples  -----------------------------------------------------------------

dat.fig.urt <- subset(dat.fig, sample_type == "Non-invasive" & organ_system == "URT")


indiv_list <- c()
for (indiv.ii in unique(dat.fig.urt$indiv_sample)) {
  indiv.sub <- subset(dat.fig.urt, indiv_sample == indiv.ii & as.numeric(dpi) <= 12)
  if (nrow(indiv.sub) >= 3 & 1 %in% indiv.sub$pos_total){
    indiv_list <- c(indiv_list, indiv.ii)
  }
}
set.seed(11)
dat.fig.urt.small <- subset(dat.fig.urt, indiv_sample %in% sample(indiv_list, 12))


fig.sg.traj.urt <- ggplot(dat.fig.urt.small) + 
  
  geom_hline(data = subset(dat.fig.urt.small, lodq_total == -9),
             aes(yintercept = 0), linetype = "dotted",
             size = 0.5, color = "grey22", alpha = 0.7) +
  geom_hline(data = subset(dat.fig.urt.small, lodq_total != -9),
             aes(yintercept = lodq_total), linetype = "dotted",
             size = 0.5, color = "#0A2472", alpha = 0.7) +
  geom_hline(data = subset(dat.fig.urt.small, lodq_total != -9),
             aes(yintercept = lodq_sg), linetype = "dotted",
             size = 0.5, color = zissou_pal[1], alpha = 0.7) +
  
  geom_line(aes(x = dpi, y = val_sg, group = indiv_sample,
                color = "sgRNA (observed)"),
            alpha = 0.8) +
  geom_point(aes(x = dpi, y = val_sg, 
                 group = indiv_sample, 
                 fill = "sgRNA (observed)"),
             alpha = 0.8, size = 1.8, shape = 23) +
  
  geom_line(aes(x = dpi, y = val_sg_pred, group = indiv_sample,
                color = "sgRNA (predicted)"),
            alpha = 0.7) +
  geom_point(aes(x = dpi, y = val_sg_pred, fill = "sgRNA (predicted)",
                 group = indiv_sample),
             alpha = 0.8, size = 1.8, shape = 24) +
  
  geom_ribbon(aes(x = dpi, ymin = val_sg, ymax = val_sg_pred,
                  group = indiv_sample),
              fill = zissou_pal[1],
                  #fill = str_remove(sample_sg, "RNA: URT")
              alpha = 0.17, show.legend = FALSE) +
  
  geom_line(aes(x = dpi, y = val_total, group = indiv_sample,
                color = "total RNA"), 
            alpha = 0.7) +
  geom_point(aes(x = dpi, y = val_total, 
                 group = indiv_sample,
                 fill = "total RNA"),
             alpha = 1, size = 1.8, shape = 21) +
  
  scale_fill_manual(values = c("total RNA" = "#0A2472",
                               "sgRNA (observed)" = zissou_pal[1], #"#3B9AB2",
                               "sgRNA (predicted)" = "#0A8549"),
                    breaks = c("total RNA", "sgRNA (observed)", 
                               "sgRNA (predicted)")) +
  scale_color_manual(values = c("total RNA" = "#0A2472",
                                "sgRNA (observed)" = zissou_pal[1], #"#3B9AB2",
                                "sgRNA (predicted)" = "#0A8549")) +
  facet_wrap(indiv ~., ncol = 6, nrow = 2) + 
  labs(y = "log10 RNA copies / sample", x = "day post infection",
       fill = "URT") +
  theme(#legend.position = c(0.94, 0.5),
    legend.position = "right", 
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.background = element_rect(fill='white'),
    legend.key=element_blank(),
    legend.key.size = unit(0.4, 'cm'),
    legend.box.spacing = unit(1, "pt"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    axis.ticks.length.x = unit(1, "pt"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
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
  scale_x_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 12),
                     labels = c("0", seq(2, 10, by = 2)), 
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(from = 0, to = 10, by = 2), 
                     minor_breaks = seq(1, 9, 2),
                     limits = c(-1.5, 10),
                     expand = c(0, 0)); fig.sg.traj.urt


# LRT samples  -----------------------------------------------------------------

dat.fig.lrt <- subset(dat.fig, sample_type == "Non-invasive" & organ_system == "LRT")

indiv_list <- c()
for (indiv.ii in unique(dat.fig.lrt$indiv_sample)) {
  indiv.sub <- subset(dat.fig.lrt, indiv_sample == indiv.ii & as.numeric(dpi) <= 11)
  if (nrow(indiv.sub) >= 3 & 1 %in% indiv.sub$pos_total){
    indiv_list <- c(indiv_list, indiv.ii)
  }
}
set.seed(111)
dat.fig.lrt.small <- subset(dat.fig.lrt, indiv_sample %in% sample(indiv_list, 12) &
                              dpi <= 12)


fig.sg.traj.lrt <- ggplot(dat.fig.lrt.small) + 
  
  geom_hline(data = subset(dat.fig.lrt.small, lodq_total == -9),
             aes(yintercept = 0), linetype = "dotted",
             size = 0.5, color = "grey22", alpha = 0.7) +
  geom_hline(data = subset(dat.fig.lrt.small, lodq_total != -9),
             aes(yintercept = lodq_total), linetype = "dotted",
             size = 0.5, color = zissou_pal[12], alpha = 0.7) +
  geom_hline(data = subset(dat.fig.lrt.small, lodq_sg != -9),
             aes(yintercept = lodq_sg), linetype = "dotted",
             size = 0.5, color = zissou_pal[10], alpha = 0.7) +
  
  geom_line(aes(x = dpi, y = val_sg, group = indiv_sample,
                color = "sgRNA (observed)"),
            alpha = 0.8) +
  geom_point(aes(x = dpi, y = val_sg, 
                 group = indiv_sample, 
                 fill = "sgRNA (observed)"),
             alpha = 0.8, size = 1.8, shape = 23) +
  
  geom_line(aes(x = dpi, y = val_sg_pred, group = indiv_sample,
                color = "sgRNA (predicted)"),
            alpha = 0.7) +
  geom_point(aes(x = dpi, y = val_sg_pred, fill = "sgRNA (predicted)",
                 group = indiv_sample),
             alpha = 0.8, size = 1.8, shape = 24) +
  
  geom_ribbon(aes(x = dpi, ymin = val_sg, ymax = val_sg_pred,
                  group = indiv_sample),
              fill = zissou_pal[10],
              #fill = str_remove(sample_sg, "RNA: URT")
              alpha = 0.17, show.legend = FALSE) +
  
  geom_line(aes(x = dpi, y = val_total, group = indiv_sample,
                color = "total RNA"), 
            alpha = 0.7) +
  geom_point(aes(x = dpi, y = val_total, 
                 group = indiv_sample,
                 fill = "total RNA"),
             alpha = 1, size = 1.8, shape = 21) +

  scale_fill_manual(values = c("total RNA" = zissou_pal[12],
                               "sgRNA (observed)" = zissou_pal[10], #"#3B9AB2",
                               "sgRNA (predicted)" = zissou_pal[7]),
                    breaks = c("total RNA", "sgRNA (observed)", 
                               "sgRNA (predicted)")) +
  scale_color_manual(values = c("total RNA" = zissou_pal[12],
                                "sgRNA (observed)" = zissou_pal[10],  #"#3B9AB2",
                                "sgRNA (predicted)" = zissou_pal[7])) +
  
  facet_wrap(indiv ~., ncol = 6, nrow = 2) + 
  labs(y = "log10 RNA copies / sample", x = "day post infection",
       fill = "LRT") +
  theme(#legend.position = c(0.94, 0.5),
    legend.position = "right", 
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.background = element_rect(fill='white'),
    legend.key=element_blank(),
    legend.key.size = unit(0.4, 'cm'),
    legend.box.spacing = unit(1, "pt"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    axis.title.y = element_blank(),
    axis.ticks.length.x = unit(1, "pt"),
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
    plot.margin = margin(0.5, 2, 2, 2)) + 
  guides(color = "none", 
         fill = guide_legend(override.aes = list(size = 1.8, shape = c(21, 23, 24)))) +
  scale_x_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 12),
                     labels = c("0", seq(2, 10, by = 2)), 
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(from = 0, to = 10, by = 2), 
                     minor_breaks = seq(1, 9, 2),
                     limits = c(-1.5, 10),
                     expand = c(0, 0)); fig.sg.traj.lrt

# Combine -----------------------------------------------------------------------

fig.sg.traj <- fig.sg.traj.urt + fig.sg.traj.lrt + plot_layout(nrow = 2); fig.sg.traj

# Generate text-only figure for joint labels
y_title <- ggplot() + 
  annotate(geom = "text", x = 1, y = 1, size = 3,
           label = "log10 RNA copies / sample", angle = 90) +
  theme_void()

# Combine into one figure
fig <- y_title + fig.sg.traj + plot_layout(ncol = 2, widths = c(0.035, 1.1)); fig

ggsave("./outputs/figures/fig4-sgRNA-traj-reconstructions.tiff", 
       fig, width = 174, height = 80, units = "mm")


