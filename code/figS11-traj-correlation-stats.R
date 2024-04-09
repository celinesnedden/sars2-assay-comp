# Generates supplemental figure 11, which includes various statistics 
#   correlating individual PCR and culture trajectories


# Prep environment -------------------------------------------------------------

# Install & load color palette & ggplot 
req_pkgs <- c("wesanderson", "ggplot2", "ggforce", "stringr",
              "patchwork", "gridExtra", "ltm")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependences = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))

# Extract desired coloredpalette
zissou_pal <- wes_palette("Zissou1", 12, type = "continuous")

# Minor function for convenience
`%notin%` <- Negate(`%in%`)


# Prep data ---------------------------------------------------------------------

dat.fig <- read.csv("./data/clean-data.csv")

# Convert columns to numerics
dat.fig$val_sg <- as.numeric(dat.fig$val_sg)
dat.fig$val_total <- as.numeric(dat.fig$val_total)

# Combine indiv names with sample rep to get the right connecting lines
dat.fig$indiv_sample <- paste0(dat.fig$indiv, "_", dat.fig$sample_rep, 
                               "_", dat.fig$tg_idx)


# Set values for samples below LOD 
dat.fig$val_sg[dat.fig$cens_sg == 1] <- -0.5
dat.fig$val_total[dat.fig$cens_total == 1] <- 0


# Indicator 
dat.fig$sample_sg <- paste0("sgRNA: ", dat.fig$organ_system)
dat.fig$sample_total <- paste0("total RNA: ", dat.fig$organ_system)
dat.fig$sample_culture <- paste0("culture: ", dat.fig$organ_system)


# PCR Statistics ---------------------------------------------------------------

dat.fig.pcr <- subset(dat.fig, pos_total %in% c(0, 1) & pos_sg %in% c(0, 1))

# note: various statistics reported in the manuscript are also calculated below
percent_tpos_sgneg <- nrow(subset(dat.fig.pcr, pos_total == 1 & pos_sg == 0))/nrow(subset(dat.fig.pcr, pos_total == 1)) * 100
cat("The percent of totRNA postiive samples that are sgRNA negative is: ", percent_tpos_sgneg)


## 11A: Absolute difference histogram -------------------------------------------

dat.fig.pcr$abs_diff <- NA
dat.fig.pcr$abs_diff[dat.fig.pcr$cens_sg == 1 & dat.fig.pcr$cens_total != 1] <- dat.fig.pcr$val_total[dat.fig.pcr$cens_sg == 1 & dat.fig.pcr$cens_total != 1]
dat.fig.pcr$abs_diff[dat.fig.pcr$cens_sg == 1 & dat.fig.pcr$cens_total == 1] <- 0
dat.fig.pcr$abs_diff[is.na(dat.fig.pcr$abs_diff)] <- dat.fig.pcr$val_total[is.na(dat.fig.pcr$abs_diff)] - dat.fig.pcr$val_sg[is.na(dat.fig.pcr$abs_diff)]

dat.fig.pcr$tg_name <- NA
dat.fig.pcr$tg_name[dat.fig.pcr$tg_idx == 1] <- "T↑ SG↑"
dat.fig.pcr$tg_name[dat.fig.pcr$tg_idx == 2] <- "T↓ SG↑"
dat.fig.pcr$tg_name[dat.fig.pcr$tg_idx == 3] <- "T↑ SG↓"
dat.fig.pcr$tg_name[dat.fig.pcr$tg_idx == 4] <- "T↓ SG↓"
dat.fig.pcr$tg_name <- factor(dat.fig.pcr$tg_name)

figS11A <- ggplot(data = subset(dat.fig.pcr, cens_sg == 0 & cens_total == 0)) +
  geom_histogram(aes(x = abs_diff, fill = tg_name),
                 color = "black",
                 position = "stack", alpha = 1,
                 breaks = seq(from = -2, to = 7, by = 0.5)) +
  geom_vline(xintercept = median(subset(dat.fig.pcr, cens_sg == 0 & cens_total == 0)$abs_diff),
             color = "#5e468e", linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("T↑ SG↑" = "#456187" ,
                               "T↓ SG↑" =  "#F39237",
                               "T↑ SG↓" = "#52A4BC" ,
                               "T↓ SG↓" = "#A1C084"),
                    drop = FALSE) +
  labs(fill = "Target Gene", x = "log10 totRNA - log10 sgRNA") +
  facet_wrap(.~ tg_name, nrow = 4, drop = FALSE) +
  scale_x_continuous(breaks = seq(-2, 7, 1),
                     limits = c(-2, 7), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 80, 20), 
                     limits = c(0, 85), expand = c(0, 0)) +
  theme(legend.position = c(0.765, 0.832), 
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, 'cm'),
        legend.background = element_rect(fill = 'white', color = "black"),
        legend.key = element_blank(),
        text = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(2, 2, 2, 2)); figS11A


# Median difference
med_diff <- median(dat.fig.pcr$abs_diff[dat.fig.pcr$cens_sg == 0 & dat.fig.pcr$cens_total == 0])
cat("The median difference between totRNA and sgRNA, when both are detectable, is: ", med_diff)


## 11B: tRNA, when sgRNA < LOD  -------------------------------------------

dat.fig.pcrS11B <- subset(dat.fig.pcr, pos_sg == 0 & cens_total == 0)

figS11B <- ggplot(dat.fig.pcrS11B) +
  geom_histogram(aes(x = val_total, 
                     fill = as.character(tg_name)),
                 breaks = seq(-2, 7, 0.5),
                 color = "black",
                 position = "stack", alpha = 1) +
  geom_vline(xintercept = median(dat.fig.pcrS11B$val_total),
             color = "#5e468e", linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("T↑ SG↑" = "#456187" ,
                               "T↓ SG↑" =  "#F39237",
                               "T↑ SG↓" = "#52A4BC" ,
                               "T↓ SG↓" = "#A1C084")) +
  labs(x = "log10 totRNA (sgRNA negative)") +
  facet_wrap(.~ tg_name, nrow = 4) +
  scale_x_continuous(breaks = seq(-2, 7, 1),
                     limits = c(-2, 7), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 32, 10), 
                     limits = c(0, 32), expand = c(0, 0)) +
  theme(legend.position = "none", 
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, 'cm'),
        legend.background = element_rect(fill = 'white', color = "black"),
        legend.key = element_blank(),
        text = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(2, 2, 2, 2)); figS11B

# Range statistic
cat("The range of totRNA values when sgRNA < LOD is: ", range(dat.fig.pcrS11B$val_total))


## 11C: Pearson Correlation for Positive vals -----------------------------------

dat.fig.pcrS11C <- subset(dat.fig.pcr, sample_type == "Non-invasive")

corr.df <- data.frame(indiv = NA, pearsonr = NA, r_sq = NA, slope = NA, tg_idx = NA)
for (indiv.ii in unique(dat.fig.pcrS11C$indiv_sample)) {
  indiv.subset <- subset(dat.fig.pcrS11C, indiv_sample == indiv.ii)
  if (TRUE %in% (indiv.subset$cens_total == 0 & indiv.subset$cens_sg == 0)) {
    sample.subset <- subset(indiv.subset, cens_total == 0 & cens_sg == 0)
    if (nrow(sample.subset) > 2) {
      pearsonr <- cor(sample.subset$val_total, sample.subset$val_sg)
      model <- lm(val_sg ~ val_total, data = sample.subset)
      rsquared <- summary(model)$r.squared
      slope <- model$coefficients[2]
      if (pearsonr < 0) {print(indiv.ii)}
    }
    corr.df <- rbind(corr.df, c(indiv.ii, pearsonr, rsquared, slope, unique(sample.subset$tg_idx)))
  }
}

# Individual correlations
corr.df <- corr.df[-1, ]


figS11C <- ggplot(data = corr.df) +
  geom_histogram(aes(x = as.numeric(pearsonr)),
                 fill = zissou_pal[8], color = "black",
                 position = "stack", alpha = 1,
                 breaks = seq(from = 0, to = 1, by = 0.1)) +
  geom_vline(xintercept = median(as.numeric(corr.df$pearsonr)),
             color = "#5e468e", linetype = "dashed", size = 1) +
  labs(fill = "Target Gene", x = "Pearson R") +
  scale_x_continuous(breaks = seq(from = -1, to = 1, by = 0.1),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(from = 0, to = 100, by = 20),
                     limits = c(0, 100), expand = c(0, 0)) +
  theme(text = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor.y = element_line(size = 0.1, linetype = 'solid',
                                          colour = "light grey"),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(2, 2, 2, 2)); figS11C

cat("The median correlation between sgRNA & totRNA quantities, when both 
    are detectable and for individuals with at least three sample times, is: ", 
    median(as.numeric(corr.df$pearsonr), na.rm = TRUE))


## 11D: First negative comparisons  ---------------------------------------------

# loops over all individual trajectories and finds the first RNA negatives
first.neg.df <- data.frame(indiv = NA, total_dpi = NA, sg_dpi = NA)
inf_indiv <- c()
for (indiv.ii in unique(dat.fig.pcr$indiv_sample[dat.fig.pcr$sample_type == "Non-invasive"])) {
  indiv.subset <- subset(dat.fig.pcr, indiv_sample == indiv.ii)
  indiv.subset <- indiv.subset[order(indiv.subset$dpi), ]
  
  # Check total RNA first
  total_cens <- indiv.subset$cens_total[indiv.subset$dpi != 1]
  total_dpi_cens <- indiv.subset$dpi[indiv.subset$cens_total == 1 & indiv.subset$dpi != 1]
  total_switch <- str_detect(paste0(total_cens, collapse = ""), "10")
  
  if (length(total_dpi_cens) == 0){
    # No samples < LOD
    print(paste0(indiv.ii, " has no total RNA samples < LOD", sep = ""))
    total_dpi <- NA
  }
  else if (total_switch == FALSE){
    total_dpi <- min(total_dpi_cens)
  }
  else if (total_switch == TRUE) {
    first_dpi <- min(indiv.subset$dpi[indiv.subset$cens_total == 1 & indiv.subset$dpi > 1])
    
    if (length(indiv.subset$dpi[indiv.subset$cens_total == 1 & 
                                indiv.subset$dpi > first_dpi]) == 0) {
      # No other negative sample besides switch
      print(paste0(indiv.ii, " has no total RNA samples < LOD after the first switch", sep = ""))
      total_dpi <- NA
    }
    else {
      next_dpi <- min(indiv.subset$dpi[indiv.subset$cens_total == 1 & 
                                         indiv.subset$dpi > first_dpi])
      total_dpi <- next_dpi
    }
  }
  
  # Check sgRNA
  sg_cens <- indiv.subset$cens_sg[indiv.subset$dpi != 1]
  sg_dpi_cens <- indiv.subset$dpi[indiv.subset$cens_sg == 1 & indiv.subset$dpi != 1]
  sg_switch <- str_detect(paste0(sg_cens, collapse = ""), "10")
  
  if (length(sg_dpi_cens) == 0) {
    # No samples < LOD
    print(paste0(indiv.ii, " has no sgRNA samples < LOD", sep = ""))
    sg_dpi <- NA
  }
  else if (sg_switch == FALSE){
    sg_dpi <- min(sg_dpi_cens)
  }
  else if (sg_switch == TRUE) {
    first_dpi <- min(indiv.subset$dpi[indiv.subset$cens_sg == 1 & indiv.subset$dpi != 1])
    
    if (length(indiv.subset$dpi[indiv.subset$cens_sg == 1 & 
                                indiv.subset$dpi > first_dpi]) == 0) {
      print(paste0(indiv.ii, " has no sgRNA samples < LOD after the first switch", sep = ""))
      sg_dpi <- NA
    }
    else {
      next_dpi <- min(indiv.subset$dpi[indiv.subset$cens_sg == 1 & 
                                         indiv.subset$dpi > first_dpi])
      sg_dpi <- next_dpi
    }
  }
  first.neg.df <- rbind(first.neg.df, c(indiv.ii, total_dpi, sg_dpi))
}

# Calculates the max observed first RNA negative, to set figure axes labels
max_dpi <- max(as.numeric(c(subset(first.neg.df, sg_dpi != Inf)$sg_dpi, 
                            subset(first.neg.df, total_dpi != Inf)$total_dpi))) + 2

# Samples that are never observed < LOD are plotted at the max dpi 
first.neg.df$total_dpi[is.na(first.neg.df$total_dpi) & !is.na(first.neg.df$sg_dpi)] <- max_dpi

# Converting the matrix to the correct form for ggplot
first.neg.matrix <- data.frame(sg_dpi = NA,
                               total_dpi = NA,
                               n = NA)
for (row_num in 1:max_dpi) {
  for (col_num in 1:max_dpi) {
    entry.subset <- subset(first.neg.df, sg_dpi == row_num & total_dpi == col_num)
    first.neg.matrix <- rbind(first.neg.matrix, 
                              c(row_num, col_num, nrow(entry.subset)))
  }
}

first.neg.matrix$n[first.neg.matrix$n == 0] <- NA

figS11D <- ggplot(first.neg.matrix) +
  annotate(geom = "rect", xmin = 29, xmax = 31, ymin = -Inf, ymax = Inf,
           alpha = 0.5, fill = "light grey", color = "black", size = 0.4) +
  geom_abline(linetype = "dashed", color = "grey", alpha = 0.5) +
  geom_point(aes(x = as.numeric(total_dpi), y = as.numeric(sg_dpi),
                 size = n),
             fill = zissou_pal[10], shape = 21) +
  scale_size_area(breaks = c(1, 10, 20)) +
  scale_x_continuous(breaks = c(seq(1, max_dpi - 2, 2), 30), 
                     labels = c(seq(1, max_dpi - 2, 2),  "None"),
                     limits = c(0, 31), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(1, max_dpi - 2, 2),
                     limits = c(1, max_dpi -2)) +
  labs(x = "First totRNA Negative (dpi)",
       y = "First sgRNA Negative (dpi)",
       size = "Number of\nIndividuals") +
  theme(legend.position =  c(0.15, 0.80),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.05, 'cm'),
        legend.background = element_rect(fill = 'white', color = "black"),
        legend.key = element_blank(),
        text = element_text(size = 7),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(2, 2, 2, 2)); figS11D



# Culture Statistics -----------------------------------------------------------

## Simple reported stats ------------------------------------------------------

dat.fig.cul <- subset(dat.fig, !is.na(pos_inf) & !is.na(pos_total)) 

# Percent of samples where total RNA and inf virus results disagree 
cul.tot.disagree <- sum(dat.fig.cul$pos_total != dat.fig.cul$pos_inf)/nrow(dat.fig.cul) * 100
cat("TotRNA and culture positivity results disagree for ", cul.tot.disagree,
    "% of all samples")

# Percent of totRNA positive samples where results disagree
dat.cul.tot.pos <- subset(dat.fig.cul, pos_total == 1)
cul.totpos.disagree <- sum(dat.cul.tot.pos$pos_total != dat.cul.tot.pos$pos_inf) / nrow(dat.cul.tot.pos) * 100
cat("TotRNA and culture positivity results disagree for ", cul.totpos.disagree,
    "% of totRNA-positive samples")

# Ranges for sgRNA data
dat.fig.cul.sg <- subset(dat.fig, !is.na(pos_inf) & !is.na(pos_sg)) 
cat("The max sgRNA values for culture positive samples is: ", 
    max(dat.fig.cul$val_sg[dat.fig.cul$pos_inf == 1], na.rm = TRUE))
cat("The max of sgRNA values for culture negative samples is: ", 
    max(dat.fig.cul$val_sg[dat.fig.cul$pos_inf == 0], na.rm = TRUE))


## 11E: tRNA values when culture pos  ------------------------------------------

dat.fig.cul$val_total[dat.fig.cul$pos_total == 0] <- 0

cat("The maximum totRNA value for culture positive samples is: ",
    max(subset(dat.fig.cul, pos_inf == 1 & pos_total == 1)$val_total))
cat("The maximum totRNA value for culture negative samples is: ",
    max(subset(dat.fig.cul, pos_inf == 0 & pos_total == 1)$val_total))
cat("The minimum totRNA value for culture positive samples is: ",
    min(subset(dat.fig.cul, pos_inf == 1 & pos_total == 1)$val_total))

n_neg <- nrow(subset(dat.fig.cul, pos_inf == 1 & pos_total == 0))
figS11E <- ggplot(data = subset(dat.fig.cul, pos_inf == 1 & pos_total == 1)) +
  geom_histogram(aes(x = val_total, fill = culture_assay), 
                 #fill = zissou_pal[7],
                 color = "black",
                 position = "stack", alpha = 1,
                 breaks = seq(0, 12.5, 0.5)) +
  geom_vline(aes(xintercept = median(val_total)),
             color = "#5e468e", linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("plaque" = "#52A4BC",
                               "TCID50" = "#ABC699")) +
  labs(fill = "Assay", x = "log10 totRNA copies (culture pos.)") + 
  facet_wrap(.~ culture_assay, nrow = 2) +
  scale_x_continuous(breaks = seq(-2, 12, 2),
                     limits = c(-2, 12.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 30, 10), 
                     limits = c(0, 38), expand = c(0, 0)) +
  theme(legend.position = c(0.2, 0.88), 
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, 'cm'),
        legend.background = element_rect(fill = 'white', color = "black"),
        legend.spacing.y = unit(0.05, "cm"),
        legend.key = element_blank(),
        text = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(2, 2, 2, 2)); figS11E



## 11F: tRNA values when culture neg  ------------------------------------------

figS11F <- ggplot(data = subset(dat.fig.cul, pos_inf == 0 & pos_total == 1)) +
  geom_histogram(aes(x = val_total, fill = culture_assay), 
                 color = "black",
                 position = "stack", alpha = 1,
                 breaks = seq(-2, 12.5, 0.5)) +
  geom_vline(xintercept = median(subset(dat.fig.cul, pos_inf == 0 & pos_total == 1)$val_total, na.rm = TRUE),
             color = "#5e468e", linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("plaque" = "#52A4BC",
                               "TCID50" = "#ABC699")) +
  labs(x = "log10 totRNA copies (culture neg.)") + 
  facet_wrap(.~ culture_assay, nrow = 2) +
  scale_x_continuous(breaks = seq(-2, 12, 2),
                     limits = c(-2, 12.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 40, 10), 
                     limits = c(0, 38), expand = c(0, 0)) +
  theme(legend.position = "none", 
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, 'cm'),
        legend.spacing.y = unit(0.05, "cm"),
        legend.background = element_rect(fill = 'white', color = "black"),
        legend.key = element_blank(),
        text = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(2, 2, 2, 2)); figS11F



## 11G: First negative comparisons  ---------------------------------------------

dat.fig.cul$cens_inf[dat.fig.cul$pos_inf == 0] <- 1
dat.fig.cul$cens_inf[dat.fig.cul$pos_inf == 1] <- 0

# loops over all individual trajectories and finds the first negatives
first.neg.df <- data.frame(indiv = NA, total_dpi = NA, inf_dpi = NA)
inf_indiv <- c()
for (indiv.ii in unique(dat.fig.cul$indiv_sample[dat.fig.cul$sample_type == "Non-invasive"])) {
  indiv.subset <- subset(dat.fig.cul, indiv_sample == indiv.ii)
  indiv.subset <- indiv.subset[order(indiv.subset$dpi), ]
  
  # Check total RNA first
  total_cens <- indiv.subset$cens_total[indiv.subset$dpi != 1]
  total_dpi_cens <- indiv.subset$dpi[indiv.subset$cens_total == 1 & indiv.subset$dpi != 1]
  total_switch <- str_detect(paste0(total_cens, collapse = ""), "10")
  
  if (length(total_dpi_cens) == 0){
    # No samples < LOD
    print(paste0(indiv.ii, " has no total RNA samples < LOD", sep = ""))
    total_dpi <- NA
  }
  else if (total_switch == FALSE){
    total_dpi <- min(total_dpi_cens)
  }
  else if (total_switch == TRUE) {
    first_dpi <- min(indiv.subset$dpi[indiv.subset$cens_total == 1 & indiv.subset$dpi > 1])
    
    if (length(indiv.subset$dpi[indiv.subset$cens_total == 1 & 
                                indiv.subset$dpi > first_dpi]) == 0) {
      # No other negative sample besides switch
      print(paste0(indiv.ii, " has no total RNA samples < LOD after the first switch", sep = ""))
      total_dpi <- NA
    }
    else {
      next_dpi <- min(indiv.subset$dpi[indiv.subset$cens_total == 1 & 
                                         indiv.subset$dpi > first_dpi])
      total_dpi <- next_dpi
    }
  }
  
  # Check sgRNA
  inf_cens <- indiv.subset$cens_inf[indiv.subset$dpi != 1]
  inf_dpi_cens <- indiv.subset$dpi[indiv.subset$cens_inf == 1 & indiv.subset$dpi != 1]
  inf_switch <- str_detect(paste0(inf_cens, collapse = ""), "10")
  
  if (length(inf_dpi_cens) == 0) {
    # No samples < LOD
    print(paste0(indiv.ii, " has no negative culure samples", sep = ""))
    inf_dpi <- NA
  }
  else if (inf_switch == FALSE){
    inf_dpi <- min(inf_dpi_cens)
  }
  else if (inf_switch == TRUE) {
    first_dpi <- min(indiv.subset$dpi[indiv.subset$cens_inf == 1 & indiv.subset$dpi != 1])
    
    if (length(indiv.subset$dpi[indiv.subset$cens_inf == 1 & 
                                indiv.subset$dpi > first_dpi]) == 0) {
      print(paste0(indiv.ii, " has no negative culture samples after the first switch", sep = ""))
      inf_dpi <- NA
    }
    else {
      next_dpi <- min(indiv.subset$dpi[indiv.subset$cens_inf == 1 & 
                                         indiv.subset$dpi > first_dpi])
      inf_dpi <- next_dpi
    }
  }
  first.neg.df <- rbind(first.neg.df, c(indiv.ii, total_dpi, inf_dpi))
}

# Calculates the max observed first RNA negative, to set figure axes labels
max_dpi <- max(as.numeric(c(subset(first.neg.df, inf_dpi != Inf)$inf_dpi, 
                            subset(first.neg.df, total_dpi != Inf)$total_dpi))) + 2

# Samples that are never observed < LOD are plotted at the max dpi 
first.neg.df$total_dpi[is.na(first.neg.df$total_dpi) & !is.na(first.neg.df$inf_dpi)] <- max_dpi

# Converting the matrix to the correct form for ggplot
first.neg.matrix <- data.frame(inf_dpi = NA,
                               total_dpi = NA,
                               n = NA)
for (row_num in 1:max_dpi) {
  for (col_num in 1:max_dpi) {
    entry.subset <- subset(first.neg.df, inf_dpi == row_num & total_dpi == col_num)
    first.neg.matrix <- rbind(first.neg.matrix, 
                              c(row_num, col_num, nrow(entry.subset)))
  }
}

first.neg.matrix$n[first.neg.matrix$n == 0] <- NA

figS11G <- ggplot(first.neg.df) +
  annotate(geom = "rect", xmin = 29, xmax = 31, ymin = -Inf, ymax = Inf,
           alpha = 0.5, fill = "light grey", color = "black", size = 0.4) +
  geom_abline(linetype = "dashed", color = "grey", alpha = 0.5) +
  geom_count(aes(x = as.numeric(total_dpi), y = as.numeric(inf_dpi)),
             fill = zissou_pal[10], shape = 21) +
  scale_size_area(limits = c(1, 40), breaks = c(1, 20, 40)) +
  scale_x_continuous(breaks = c(seq(1, 27, 2), 30), 
                     labels = c(seq(1, 27, 2), "None"),
                     limits = c(0, max_dpi + 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(1, max_dpi - 2, 2),
                     limits = c(1, max_dpi -2)) +
  labs(x = "First totRNA Negative (dpi)",
       y = "First Culture Negative (dpi)",
       size = "Number of\nIndividuals") +
  theme(legend.position = c(0.15, 0.8),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.size = unit(0.05, 'cm'),
        legend.background = element_rect(fill = 'white', color = "black"),
        legend.key = element_blank(),
        text = element_text(size = 7),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor = element_blank()); figS11G

# 1 individual where culture RNA first negative is after totRNA
subset(first.neg.df, as.numeric(total_dpi) < as.numeric(inf_dpi))


## 11H: First positive comparisons  ---------------------------------------------

# loops over all individual trajectories and finds the first positives
first.pos.df <- data.frame(indiv = NA, total_dpi = NA, inf_dpi = NA)
inf_indiv <- c()
for (indiv.ii in unique(dat.fig.cul$indiv_sample)) {
  indiv.sub <- subset(dat.fig.cul, indiv_sample == indiv.ii)
  
  if (unique(indiv.sub$sample_type) == "Non-invasive" & (1 %in% indiv.sub$pos_total | 1 %in% indiv.sub$pos_inf)) {
    first_total <- min(indiv.sub$dpi[indiv.sub$pos_total == 1])
    first_inf <- min(indiv.sub$dpi[indiv.sub$pos_inf == 1])
    if (first_inf == Inf){
      inf_indiv <- c(inf_indiv, indiv.ii)
      first_inf <- "Never inf pos."
    }
    if (first_total == Inf){
      inf_indiv <- c(inf_indiv, indiv.ii)
      first_total <- "Never total pos."
    }
    first.pos.df <- rbind(first.pos.df, c(indiv.ii, first_total, first_inf))
  }
}

# Individuals where culture is positive before totRNA
subset(first.pos.df, as.numeric(total_dpi) > as.numeric(inf_dpi))
cat("The number of individuals where culture positivity preceded totRNA
        positivity is: ", nrow(subset(first.pos.df, as.numeric(total_dpi) > as.numeric(inf_dpi))))


# Subset where culture is positive and totRNA is never positive
subset(first.pos.df, total_dpi == "Never total pos." & inf_dpi != "Never inf. pos")
cat("The number of individuals where culture positivity preceded totRNA
        positivity is: ", nrow(subset(first.pos.df, total_dpi == "Never total pos." & inf_dpi != "Never inf. pos")))

# Set never positives to right hand side of plot
first.pos.df$total_dpi[first.pos.df$total_dpi == "Never total pos."] <- 12
first.pos.df$inf_dpi[first.pos.df$inf_dpi == "Never inf. pos"] <- 12


figS11H <- ggplot(first.pos.df) +
  annotate(geom = "rect", xmin = 11, xmax = 13, ymin = -Inf, ymax = Inf,
           alpha = 0.5, fill = "light grey", color = "black", size = 0.4) +
  geom_abline(linetype = "dashed", color = "grey", alpha = 0.5) +
  geom_count(aes(x = as.numeric(total_dpi), y = as.numeric(inf_dpi)),
             fill = zissou_pal[1], shape = 21) +
  scale_size_area(breaks = c(1, 30, 60)) +
  scale_x_continuous(breaks = c(seq(1, 10, 1), 12), 
                     labels = c(seq(1, 10, 1), "None"),
                     limits = c(0, 13), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(1, 11, 1),
                     limits = c(1, 11)) +
  labs(x = "First totRNA Positive (dpi)",
       y = "First Culture Positive (dpi)",
       size = "Number of\nIndividuals") +
  theme(legend.position = c(0.7, 0.2),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.size = unit(0.05, 'cm'),
        legend.background = element_rect(fill = 'white', color = "black"),
        legend.key = element_blank(),
        text = element_text(size = 7),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor = element_blank()); figS11H




# Combine into 1 plot ----------------------------------------------------------

figS11 <- figS11A + figS11B + figS11C + figS11D + figS11E + figS11F + 
                   figS11G + figS11H + plot_layout(ncol = 4, widths = c(0.8, 0.8, 1, 1)) +
  plot_annotation(tag_levels = "A"); figS11


# Save -------------------------------------------------------------------------

ggsave("./outputs/figures/figS11-traj-correlation-stats.tiff", figS11,
       width = 10, height = 6, units = "in")

