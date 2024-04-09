# Generates Figure 1, which includes:
#   - Example individual trajectories
#   - Sample distribution venn diagram

# NOTE: It will not be possible to create the venn diagram (Fig1B) using this
#   script until the full database is published in a subsequent manuscript


# Prep environment -------------------------------------------------------------

# Install & load color palette & ggplot 
req_pkgs <- c("wesanderson", "patchwork", "eulerr", "ggpubr", "stringr")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependences = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))

# Desired color palette
zissou_pal <- wes_palette("Zissou1", 12, type = "continuous")

# Function for convenience
`%notin%` <- Negate(`%in%`)




# 1A: Example trajectories -----------------------------------------------------

### Prep the data --------------------------------------------------------------

dat.fig <- read.csv("./data/clean-data.csv")

# Subset to chosen example individuals & non-invasive samples
ex_indiv <- c("AP_7521", "ES_AGM7", "ES_AGM8")
dat.fig <- subset(dat.fig, indiv %in% ex_indiv & sample_type == "Non-invasive" &
                    organ_system %in% c("URT", "LRT"))

# Convert columns to numerics
dat.fig$val_sg <- as.numeric(dat.fig$val_sg)
dat.fig$val_total <- as.numeric(dat.fig$val_total)

# Combines IDs with sample info to get correct connecting lines
dat.fig$indiv_sample <- paste0(dat.fig$indiv, "_", dat.fig$sample_rep, 
                               "_", dat.fig$tg_idx)

# Set values for samples below LOD 
dat.fig$val_sg[dat.fig$cens_sg == 1] <- -0.5
dat.fig$val_total[dat.fig$cens_total == 1] <- 0

# Indicator to simplify plotting
dat.fig$sample_sg <- paste0("sgRNA: ", dat.fig$organ_system)
dat.fig$sample_total <- paste0("total RNA: ", dat.fig$organ_system)
dat.fig$sample_culture <- paste0("culture: ", dat.fig$organ_system)

# Update LOD columns when unknown for indication on plots
dat.fig$lodq_total_idx[dat.fig$lodq_total == -9] <- "Unknown"
dat.fig$lodq_total_idx[dat.fig$lodq_total != -9] <- "Known"
dat.fig$lodq_total[dat.fig$lodq_total == -9] <- 0
dat.fig$lodq_sg_idx[dat.fig$lodq_sg == -9] <- "Unknown"
dat.fig$lodq_sg_idx[dat.fig$lodq_sg != -9] <- "Known"
dat.fig$lodq_sg[dat.fig$lodq_sg == -9] <- 0


### URT ------------------------------------------------------------------------

# Only including nasal samples for this plot
dat.fig.urt <- subset(dat.fig, organ_system == "URT" & sample_grp == "Nasal Sample")

# Make the figure
fig1_URT <- ggplot(dat.fig.urt) + 
  
  # LOD lines
  geom_hline(data = subset(dat.fig.urt, lodq_total_idx == "Known"),
             aes(yintercept = lodq_total), linetype = "dotted",
             size = 0.7, color = "#0A2472", alpha = 0.7) +
  geom_hline(data = subset(dat.fig.urt, lodq_sg_idx == "Known"),
             aes(yintercept = lodq_sg), linetype = "dotted",
             size = 0.7, color = zissou_pal[1], alpha = 0.7) +
  geom_hline(data = subset(dat.fig.urt, lodq_total_idx == "Unknown"),
             aes(yintercept = lodq_total), linetype = "dotted",
             size = 0.7, color = "grey36", alpha = 0.7) +
  geom_hline(data = subset(dat.fig.urt, lodq_sg_idx == "Unknown"),
             aes(yintercept = lodq_sg), linetype = "dotted",
             size = 0.7, color = "grey36", alpha = 0.7) +
  geom_hline(aes(yintercept = 10), color = zissou_pal[8], size = 0.5) +
  
  # totRNA & sgRNA trajectories
  geom_ribbon(aes(x = dpi, ymin = val_sg, ymax = val_total,
                  fill = str_remove(sample_sg, "RNA: URT"), 
                  group = indiv_sample),
              alpha = 0.17, show.legend = FALSE) +
  geom_line(aes(x = dpi, y = val_sg, group = indiv_sample,
                color = str_remove(sample_sg, "RNA: URT")), 
            alpha = 0.5) +
  geom_point(aes(x = dpi, y = val_sg, 
                 group = indiv_sample, 
                 fill = str_remove(sample_sg, "RNA: URT")),
             alpha = 1, size = 2, shape = 23) +
  geom_line(aes(x = dpi, y = val_total, group = indiv_sample,
                color = str_remove(sample_total, " RNA: URT")),
            alpha = 0.5) +
  geom_point(aes(x = dpi, y = val_total, 
                 group = indiv_sample,
                 fill = str_remove(sample_total, " RNA: URT")),
             alpha = 1, size = 2, shape = 21) +
  
  # Culture data points
  geom_point(data = subset(dat.fig.urt, pos_inf == 1),
             aes(x = dpi, y = 10.75, 
                 fill = str_remove(sample_culture, ": URT")), 
             size = 2, shape = 22) +
  geom_point(data = subset(dat.fig.urt, pos_inf == 0),
             aes(x = dpi, y = 10.75), fill = "grey88",
             size = 2, shape = 22) +
  
  # Manually set color palette
  scale_fill_manual(values = c("total" = "#0A2472",
                               "sg" =    zissou_pal[1],
                               "culture" = zissou_pal[8]),
                    breaks = c("total", "sg", "culture")) +  
  scale_color_manual(values = c("total" = "#0A2472",
                                "sg" = zissou_pal[1],
                                "culture" = zissou_pal[8])) + 
  
  # Axes
  scale_x_continuous(breaks = seq(0, 15, by = 2), limits = c(-0.01, 15),
                     labels = seq(0, 15, by = 2), 
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(from = 0, to = 10, by = 2), 
                     minor_breaks = seq(1, 9, 2),
                     limits = c(-1.5, 11.5),
                     expand = c(0, 0)) +
  
  # Other aesthetics
  facet_grid(. ~ indiv) +
  labs(y = "", x = "", 
       fill = "URT") +
  guides(color = "none", 
         fill = guide_legend(override.aes = list(shape = c(21, 23, 22)))) +
  theme(legend.position = "right", 
        legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, 'cm'),
        legend.background = element_rect(fill='white'),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 7),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        strip.background  = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(size = 7.5),
        strip.text.y = element_blank(),
        plot.margin = margin(2, 2, 2, 2)); fig1_URT


### LRT ------------------------------------------------------------------------

dat.fig.lrt <- subset(dat.fig, organ_system == "LRT")

fig1_LRT <- ggplot(dat.fig.lrt) + 
  
  # LOD lines
  geom_hline(data = subset(dat.fig.lrt, lodq_total_idx == "Known"),
             aes(yintercept = lodq_total), linetype = "dotted",
             size = 0.7, color = zissou_pal[12], alpha = 0.7) +
  geom_hline(data = subset(dat.fig.lrt, lodq_total_idx == "Known"),
             aes(yintercept = lodq_sg), linetype = "dotted",
             size = 0.7, color = zissou_pal[10], alpha = 0.7) +
  geom_hline(data = subset(dat.fig.lrt, lodq_total_idx == "Unknown"),
             aes(yintercept = lodq_total), linetype = "dotted",
             size = 0.7, color = "grey36", alpha = 0.7) +
  geom_hline(data = subset(dat.fig.lrt, lodq_sg_idx == "Unknown"),
             aes(yintercept = lodq_sg), linetype = "dotted",
             size = 0.7, color = "grey36", alpha = 0.7) +
  geom_hline(aes(yintercept = 10), color = zissou_pal[8], size = 0.5) +
  
  # totRNA & sgRNA trajectories
  geom_ribbon(aes(x = dpi, ymin = val_sg, ymax = val_total,
                  fill = str_remove(sample_sg, "RNA: LRT"), 
                  group = indiv_sample),
              alpha = 0.17, show.legend = FALSE) +
  geom_line(aes(x = dpi, y = val_sg, group = indiv_sample,
                color = str_remove(sample_sg, "RNA: LRT")), 
            alpha = 0.5) +
  geom_point(aes(x = dpi, y = val_sg, 
                 group = indiv_sample, 
                 fill = str_remove(sample_sg, "RNA: LRT")),
             alpha = 1, size = 2, shape = 23) +
  geom_line(aes(x = dpi, y = val_total, group = indiv_sample,
                color = str_remove(sample_total, " RNA: LRT")),
            alpha = 0.5) +
  geom_point(aes(x = dpi, y = val_total, 
                 group = indiv_sample,
                 fill = str_remove(sample_total, " RNA: LRT")),
             alpha = 1, size = 2, shape = 21) +
  
  # Culture data points: none are negative so no code for those samples
  geom_point(data = subset(dat.fig.lrt, pos_inf == 1),
             aes(x = dpi, y = 10.75, 
                 fill = str_remove(sample_culture, ": LRT")), 
             size = 2, shape = 22) +

  # Manually set color palette
  scale_fill_manual(values = c("total" = zissou_pal[12],
                               "sg" = zissou_pal[10],
                               "culture" = zissou_pal[8]),
                    breaks = c("total", "sg", "culture")) +  
  scale_color_manual(values = c("total" = zissou_pal[12],
                                "sg" = zissou_pal[10],
                                "culture" = zissou_pal[8])) + 
  
  # Axes
  scale_x_continuous(breaks = seq(0, 15, by = 2), limits = c(-0.01, 15),
                     labels = seq(0, 15, by = 2), 
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(from = 0, to = 10, by = 2), 
                     minor_breaks = seq(1, 9, 2),
                     limits = c(-1.5, 11.5),
                     expand = c(0, 0)) +
  
  # Other aesthetics
  facet_grid(. ~ indiv) +
  labs(y = "log10 RNA copies / sample", x = "day post infection",
       fill = "LRT") +
  guides(color = "none", 
         fill = guide_legend(override.aes = list(shape = c(21, 23, 22)))) +
  theme(legend.position = "right", 
        legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, 'cm'),
        legend.background = element_rect(fill='white'),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        axis.title.y = element_blank(),
        text = element_text(size = 11),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        strip.background  = element_rect(colour = "white", fill = NA),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        plot.margin = margin(2, 2, 2, 2)); fig1_LRT

### Combine ------------------------------------------------------------------------

fig1A <- fig1_URT + fig1_LRT + 
  plot_layout(nrow = 2, ncol = 1, heights = c(2, 2)); fig1A

# Generate text-only figure for joint labels
y_title <- ggplot() + 
  annotate(geom = "text", x = 1, y = 1, size = 2.8,
           label = "log10 RNA copies / sample", angle = 90) +
  theme_void()

# Combine into one figure
fig1A <- y_title + fig1A + plot_layout(ncol = 2, widths = c(0.035, 1.1)); fig1A




# 1B: Sample distribution -----------------------------------------------------

## Get sample distributions from full dataset ---------------------------------

# Load data (NOTE: dat.full will not load until the full database has been
#   published, so the venn diagram will not be constructable. Once published
#   we will update this repo so the code below runs)
dat <- read.csv("./data/clean-data.csv")
dat.full <- read.csv("./data/full-data.csv") 


## Samples with all three measurements available
dat.3 <- subset(dat, !is.na(val_total) & !is.na(val_sg) & !is.na(val_inf))
n3_total <- nrow(dat.3)


## Samples with only two measurements available
dat.2 <- subset(dat, (!is.na(val_total) & !is.na(val_sg) & is.na(val_inf)) |
                  (!is.na(val_total) & is.na(val_sg) & !is.na(val_inf)) |
                  (is.na(val_total) & !is.na(val_sg) & !is.na(val_inf)))
n2_total <- nrow(dat.2)
n2_culture_tRNA <- nrow(subset(dat.2, !is.na(val_total) & is.na(val_sg) & !is.na(val_inf)))
n2_culture_sg <- nrow(subset(dat.2, is.na(val_total) & !is.na(val_sg) & !is.na(val_inf)))
n2_tRNA_sg <- nrow(subset(dat.2, !is.na(val_total) & !is.na(val_sg) & is.na(val_inf)))


## Samples with only one measurement available, either because:
##      1. ID names were not provided so linking samples is not possible, or
##      2. Multiple assays were not run
## These are calculated for convenience as the number of samples that do not
##      fall into the previous two categories
n1_total <- nrow(subset(dat.full, rna_type %notin% c("Unknown", "gRNA", ""))) - 3 * n3_total - 2 * n2_total 
n1_tRNA <- nrow(subset(dat.full, rna_type == "total RNA")) - n2_tRNA_sg - n2_culture_tRNA - n3_total
n1_sgRNA <- nrow(subset(dat.full, rna_type == "sgRNA")) - n2_tRNA_sg - n2_culture_sg - n3_total
n1_culture <- nrow(subset(dat.full, is.na(rna_type))) - n2_culture_tRNA - n2_culture_sg - n3_total


## Plot ------------------------------------------------------------------------

combo_venn <- c("total RNA" = n1_tRNA,
                "sgRNA" = n1_sgRNA,
                "culture" = n1_culture,
                "total RNA&sgRNA" = n2_tRNA_sg,
                "total RNA&culture" = n2_culture_tRNA,
                "sgRNA&culture" = n2_culture_sg,
                "sgRNA&culture&total RNA" = n3_total) 


# Use pacakge to get and plot venn diagram
venn_fit <- euler(combo_venn)
fig1B <- plot(venn_fit,
              fills = c("#0A2472", zissou_pal[1], zissou_pal[8]),
              edges = FALSE,
              labels = list(fontsize = 8),
              quantities = list(fontsize = 6),
              alpha = 0.8); fig1B



# Combine and save -------------------------------------------------------------

fig1 <- ggarrange(fig1A, fig1B, ncol = 3, widths = c(1.5, 0.4, 0.05),
                  labels = c("A", "B"),
                  font.label = list(size = 10, color = "black", face = "plain")) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        panel.background = element_rect(fill = "white", color = "white")); fig1


ggsave("./outputs/figures/fig1-indiv-traj-sample-dist.tiff", fig1,
       width = 174, height = 72, units = "mm")


