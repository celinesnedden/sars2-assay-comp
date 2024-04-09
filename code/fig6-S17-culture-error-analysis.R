# Generates figure 6, which includes the culture error analysis stratified by type
#   And also figure 17, with individual trajectories for certain errors

# Prep environment -------------------------------------------------------------

# Install & load packages
req_pkgs <- c("wesanderson", "ggplot2", "cmdstanr", "tidyverse", "ggridges",
              "ggpubr", "patchwork", "ggExtra")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependences = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))

`%notin%` <- Negate(`%in%`)  # for convenience

# set color palette
zissou_pal <- wes_palette("Zissou1", 12, type = "continuous")
traj_palette <- data.frame(assay = c("total", "total", 
                                     "sg_true", "sg_true", 
                                     "sg_pred", "sg_pred",
                                     "culture", "culture"),
                           pos_neg = c("all", "all", 
                                       "all", "all", 
                                       "all", "all",
                                       "pos", "neg"),
                           organ = c("URT", "LRT",
                                     "URT", "LRT",
                                     "URT", "LRT",
                                     "all", "all"),
                           hex = c("#0A2472", # total: URT
                                   zissou_pal[12], # total: LRT
                                   zissou_pal[1], # sg_true: URT
                                   zissou_pal[10], # sg_true: LRT
                                   "light blue", # sg_pred : URT
                                   "#FFBC47", # sg_pred : LRT,
                                   "#E4BC11", # culture: positive
                                   "grey62" # culture negative
                           )) 

# Load data with predictions ---------------------------------------------------

dat <- read.csv("./data/pred-culture-data.csv")
dat <- subset(dat, !is.na(pos_total))


# Add sample indicator for facetting -------------------------------------------

## Combine indiv names with sample rep to get the right connecting lines
dat$indiv_sample <- paste0(dat$indiv, "_", dat$sample_rep, 
                           "_", dat$tg_idx)


# Add TP/TN/FP/FN classifier ---------------------------------------------------

dat$pos_inf_pred_class <- NA
dat$pos_inf_pred_class[dat$pos_inf == 1 & dat$pos_inf_pred == 1] <- "TP"
dat$pos_inf_pred_class[dat$pos_inf == 0 & dat$pos_inf_pred == 0] <- "TN"
dat$pos_inf_pred_class[dat$pos_inf == 1 & dat$pos_inf_pred == 0] <- "FN"
dat$pos_inf_pred_class[dat$pos_inf == 0 & dat$pos_inf_pred == 1] <- "FP"

# Add correct/incorrect classifier ---------------------------------------------------

dat$pos_inf_pred_correct <- NA
dat$pos_inf_pred_correct[dat$pos_inf == 1 & dat$pos_inf_pred == 1] <- "Correct"
dat$pos_inf_pred_correct[dat$pos_inf == 0 & dat$pos_inf_pred == 0] <- "Correct"
dat$pos_inf_pred_correct[dat$pos_inf == 1 & dat$pos_inf_pred == 0] <- "Incorrect"
dat$pos_inf_pred_correct[dat$pos_inf == 0 & dat$pos_inf_pred == 1] <- "Incorrect"


# Calculate some statistics ----------------------------------------------------

# Percent of correct samples
dat.culture <- subset(dat, !is.na(pos_inf) & pos_total %in% c(0, 1))
nrow(subset(dat.culture, pos_inf_pred_correct == "Correct")) / nrow(dat.culture) * 100

# Percent that are false negatives
nrow(subset(dat.culture, pos_inf_pred_class == "FN")) / nrow(subset(dat.culture, pos_inf_pred_correct == "Incorrect")) * 100
nrow(subset(dat.culture, pos_inf_pred_class == "FN"))
nrow(subset(dat.culture, pos_inf_pred_correct == "Incorrect"))

# Percent of false negatives that are totRNA negative
nrow(subset(dat.culture, pos_inf_pred_class == "FN" & pos_total == 0 ))

# For distribution comparisons, please see 09-error-analysis-distribution-comparisons.R


# Classify error types ---------------------------------------------------------

dat.ni <- subset(dat, sample_type == "Non-invasive" & !is.na(pos_inf))
dat.inv <- subset(dat, sample_type == "Invasive" & !is.na(pos_inf))
dat.ni$error_type <- NA
dat.inv$error_type <- "NA"
dat.ni$error_type_edge <- NA
dat.inv$error_type_edge <- NA

for (indiv.ii in unique(dat.ni$indiv_sample)) {
  dat.indiv <- subset(dat.ni, indiv_sample == indiv.ii)
  dat.indiv <- dat.indiv %>% arrange(dpi)
  
  for (row_num in 1:nrow(dat.indiv)) {
    if (dat.indiv$pos_inf_pred_class[row_num] %in% c("FN", "FP")) {
      # Prediction is incorrect, need to screen which type of error
      
      # If an edge sample, don't need to consult surroundings
      if (dat.indiv$dpi[row_num] == min(dat.indiv$dpi) | dat.indiv$dpi[row_num] == max(dat.indiv$dpi)) {
        # Edge sample: first dpi
        dat.ni$error_type[dat.ni$indiv_sample == indiv.ii & dat.ni$dpi == dat.indiv$dpi[row_num]] <- "Edge"
        
        if (nrow(dat.indiv) > 1 & row_num != nrow(dat.indiv)) {
          if (dat.indiv$pos_inf[row_num] != dat.indiv$pos_inf[row_num + 1]) {
            dat.ni$error_type_edge[dat.ni$indiv_sample == indiv.ii & dat.ni$dpi == dat.indiv$dpi[row_num]] <- "Possible Transition"
            
          }
          
        }
        
        if (nrow(dat.indiv) > 1 & row_num == nrow(dat.indiv)) {
          if (dat.indiv$pos_inf[row_num] != dat.indiv$pos_inf[row_num - 1]) {
            dat.ni$error_type_edge[dat.ni$indiv_sample == indiv.ii & dat.ni$dpi == dat.indiv$dpi[row_num]] <- "Possible Transition"
          }
        }
        
      }
      
      else {
        # Get surrounding 1 or 2 values, if available
        # If it's not an edge, it has to have at least three samples, so we can get the 1 surrounding
        one_surrounding <- paste0(dat.indiv$pos_inf[(row_num - 1):(row_num + 1)], collapse = "")
        
        # Check if there are two previous samples
        if (row_num >= 3) {
          two_previous <- paste0(dat.indiv$pos_inf[(row_num - 3):(row_num)], collapse = "")
          two_prev_one_after <- paste0(dat.indiv$pos_inf[(row_num - 3):(row_num + 1)], collapse = "")
        }
        else {
          two_previous <- "None"
          two_prev_one_after <- "None"}
        
        # Check if there are two samples afterwards
        if (row_num <= (nrow(dat.indiv) - 2)) {
          two_after <- paste0(dat.indiv$pos_inf[(row_num):(row_num + 3)], collapse = "")
          one_prev_two_after <- paste0(dat.indiv$pos_inf[(row_num - 1):(row_num + 3)], collapse = "")
        }
        else {
          two_after <- "None"}
        

        # Now screen for blips, transitions, mid-strings,...
        if (one_surrounding %in% c("010", "101")) {
          # Blip!
          dat.ni$error_type[dat.ni$indiv_sample == indiv.ii & dat.ni$dpi == dat.indiv$dpi[row_num]] <- "Blip"
        }
        
        else if (two_previous %in% c("001", "110") | two_after %in% c("011", "100") |
                 one_surrounding %in% c("011", "110", "100", "001") |
                 two_prev_one_after %in% c("1100", "0011") |
                 one_prev_two_after %in% c("1100", "1000")) {
          # Transition
          dat.ni$error_type[dat.ni$indiv_sample == indiv.ii & dat.ni$dpi == dat.indiv$dpi[row_num]] <- "Transition"
        }
        
        else if (one_surrounding %in% c("000", "111")) {
          # Transition
          dat.ni$error_type[dat.ni$indiv_sample == indiv.ii & dat.ni$dpi == dat.indiv$dpi[row_num]] <- "Mid-String"
        }

        else {
          dat.ni$error_type[dat.ni$indiv_sample == indiv.ii & dat.ni$dpi == dat.indiv$dpi[row_num]] <- "Other"
        }
      }
    }
  }
}

# All individuals with errors --------------------------------------------------

indiv_list <- c()
indiv_onesample <- c()
for (indiv.ii in unique(dat.ni$indiv_sample)) {
  indiv.sub <- subset(dat.ni, indiv_sample == indiv.ii)
  if (str_detect(paste0(indiv.sub$error_type, collapse = ""), "Edge|Blip|Transition|Mid") &
      nrow(indiv.sub >= 2)){
    indiv_list <- c(indiv_list, indiv.ii)
  }
  else if (nrow(indiv.sub) == 1) {
    indiv_onesample <- c(indiv_onesample, indiv.ii)
  }
}


all_indivs <- unique(dat.ni$indiv_sample[dat.ni$indiv_sample %notin% indiv_onesample])

dat.ni <- subset(dat.ni, indiv_sample %notin% indiv_onesample) # remove individuals with only one sample

dat.indiv <- subset(dat.ni, indiv_sample %in% all_indivs) # all individuals
dat.none <- subset(dat.ni, indiv_sample %notin% indiv_list) # individuals with no errors

# The following should be TRUE, if not there's an issue
length(unique(c(dat.indiv$indiv_sample, dat.none$indiv_sample))) == length(all_indivs)


# Calculate more statistics ----------------------------------------------------

# Number of trajectories
length(unique(dat.none$indiv_sample)); length(all_indivs); length(unique(dat.none$indiv_sample))/length(all_indivs)

# Subset to data with only errors
dat.wrong <- subset(dat.indiv, !is.na(error_type))

# The number and percent of edge errors
nrow(subset(dat.wrong, error_type == "Edge")); nrow(dat.wrong)
nrow(subset(dat.wrong, error_type == "Edge")) / nrow(dat.wrong) * 100

nrow(subset(dat.wrong, error_type == "Edge" & error_type_edge == "Possible Transition"))

# The number and percent of transition errors
nrow(subset(dat.wrong, error_type == "Transition")); nrow(dat.wrong)
nrow(subset(dat.wrong, error_type == "Transition")) / nrow(dat.wrong) * 100

# The number and percent of data blip errors
nrow(subset(dat.wrong, error_type == "Blip")); nrow(dat.wrong)
nrow(subset(dat.wrong, error_type == "Blip")) / nrow(dat.wrong) * 100
subset(dat.wrong, error_type == "Blip")$pos_inf_pred_class

# The number and percent of prediction blip errors
nrow(subset(dat.wrong, error_type == "Mid-String")); nrow(dat.wrong)
nrow(subset(dat.wrong, error_type == "Mid-String")) / nrow(dat.wrong) * 100
subset(dat.wrong, error_type == "Mid-String")$pos_inf_pred_class
subset(dat.wrong, error_type == "Mid-String" & pos_inf_pred_class == "FP")$culture_assay

# Duration exploration ---------------------------------------------------------

indiv_duration_correct <- c()
indiv_duration_incorrect <- c()
indiv_duration_no_pos <- c()
for (indiv.ii in unique(dat.ni$indiv_sample)) {
  dat.sub <- subset(dat.ni, indiv_sample == indiv.ii)
  
  if (1 %in% dat.sub$pos_inf) {
    dat.sub <- subset(dat.sub, dpi == max(as.numeric(dat.sub$dpi)[dat.sub$pos_inf == 1]))
    if (is.na(dat.sub$error_type)) {
      indiv_duration_correct <- c(indiv_duration_correct, indiv.ii)
    }
    else {
      indiv_duration_incorrect <- c(indiv_duration_incorrect, indiv.ii)
    }
  }
  else {
    indiv_duration_no_pos <- c(indiv_duration_no_pos, indiv.ii)
  }
}

length(indiv_duration_correct) / (length(indiv_duration_correct) + length(indiv_duration_incorrect)) * 100

dat.ni$duration_correct[dat.ni$indiv_sample %in% indiv_duration_correct] <- "Correct"
dat.ni$duration_correct[dat.ni$indiv_sample %in% indiv_duration_incorrect] <- "Incorrect"


fig.duration <- ggplot(dat.ni) + 
  
  geom_line(aes(x = dpi, y = indiv_sample),
            color = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex, 
            alpha = 0.3) +
  geom_point(aes(x = dpi, y = indiv_sample, fill = as.character(pos_inf),
                 color = pos_inf_pred_correct,
                 group = indiv_sample),
             shape = 22, size = 4) +
  
  # Set colors
  scale_fill_manual(values = c("0" = subset(traj_palette, assay == "culture" & pos_neg == "neg")$hex,
                               "1" = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex,
                               "Edge" = "#BA141D",
                               "Transition" = "transparent",
                               "Blip" = "transparent",
                               "Other" = "transparent")) +  
  scale_color_manual(values = c("Correct" = "black",
                                "Incorrect" = "transparent")) +  
  # Other plot features
  facet_wrap(. ~ duration_correct, scales = "free_y") +
  labs(y = "individual", x = "", 
       fill = "Pos/Neg",
       tag = "A") +
  ggtitle("Edge") +
  theme(legend.position = "none", 
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill='white'),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        axis.text = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.ticks.y= element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        text = element_text(size = 14),
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
  scale_x_continuous(breaks = seq(0, 25, by = 3), limits = c(-0.01, 22),
                     expand = c(0, 0)); fig.duration






# Plot by error type -----------------------------------------------------------

edge_list <- c()
for (indiv.ii in unique(dat.ni$indiv_sample)) {
  indiv.sub <- subset(dat.ni, indiv_sample == indiv.ii)
  if (str_detect(paste0(indiv.sub$error_type, collapse = ""), "Edge") &
      nrow(indiv.sub) > 1){
    edge_list <- c(edge_list, indiv.ii)
  }
}

blip_list <- c()
for (indiv.ii in unique(dat.ni$indiv_sample)) {
  indiv.sub <- subset(dat.ni, indiv_sample == indiv.ii)
  if (str_detect(paste0(indiv.sub$error_type, collapse = ""), "Blip")){
    blip_list <- c(blip_list, indiv.ii)
  }
}

trans_list <- c()
for (indiv.ii in unique(dat.ni$indiv_sample)) {
  indiv.sub <- subset(dat.ni, indiv_sample == indiv.ii)
  if (str_detect(paste0(indiv.sub$error_type, collapse = ""), "Transition")){
    trans_list <- c(trans_list, indiv.ii)
  }
}

mid_list <- c()
for (indiv.ii in unique(dat.ni$indiv_sample)) {
  indiv.sub <- subset(dat.ni, indiv_sample == indiv.ii)
  if (str_detect(paste0(indiv.sub$error_type, collapse = ""), "Mid-String")){
    mid_list <- c(mid_list, indiv.ii)
  }
}


other_list <- c()
for (indiv.ii in unique(dat.ni$indiv_sample)) {
  indiv.sub <- subset(dat.ni, indiv_sample == indiv.ii)
  if (str_detect(paste0(indiv.sub$error_type, collapse = ""), "Other")){
    other_list <- c(other_list, indiv.ii)
  }
}



## Edge ------------------------------------------------------------------------

dat.edge <- subset(dat.ni, indiv_sample %in% edge_list)

figSX_edge <- ggplot(dat.edge) + 
  
  # tRNA & sgRNA trajectories
  #geom_point(data = subset(dat.edge, is.na(error_type)),
  #           aes(x = dpi, y = indiv_sample,
  #               group = indiv_sample),
  #           color = "transparent",
  #           fill = "black",
  #           shape = 22, size = 6, alpha = 1) +
  geom_line(aes(x = dpi, y = indiv_sample),
            color = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex, 
            alpha = 0.3) +
  geom_point(data = subset(dat.edge, error_type == "Edge"),
             aes(x = dpi, y = indiv_sample, fill = as.character(error_type),
                 group = indiv_sample),
             color = "transparent",
             shape = 22, size = 6, alpha = 1) +
  geom_point(aes(x = dpi, y = indiv_sample, fill = as.character(pos_inf),
                 color = pos_inf_pred_correct,
                 group = indiv_sample),
             shape = 22, size = 4) +
  
  # Set colors
  scale_fill_manual(values = c("0" = subset(traj_palette, assay == "culture" & pos_neg == "neg")$hex,
                               "1" = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex,
                               "Edge" = "#BA141D",
                               "Transition" = "transparent",
                               "Blip" = "transparent",
                               "Other" = "transparent")) +  
  scale_color_manual(values = c("Correct" = "black",
                                "Incorrect" = "transparent")) +  
  # Other plot features
  labs(y = "individual", x = "", 
       fill = "Pos/Neg",
       tag = "A") +
  ggtitle("Edge") +
  theme(legend.position = "none", 
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill='white'),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        axis.text = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.ticks.y= element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        text = element_text(size = 14),
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
  scale_x_continuous(breaks = seq(0, 25, by = 3), limits = c(-0.01, 22),
                     expand = c(0, 0)); figSX_edge



## Data Blip -------------------------------------------------------------------

### Main figure panel ----------------------------------------------------------

dat.blip <- subset(dat.ni, indiv_sample %in% blip_list)

figSX_blip <- ggplot(dat.blip) + 
  
  
  # tRNA & sgRNA trajectories
  geom_line(aes(x = dpi, y = indiv_sample),
            color = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex, 
            alpha = 0.3) +
  geom_point(data = subset(dat.blip, error_type == "Blip"),
             aes(x = dpi, y = indiv_sample, fill = as.character(error_type),
                 group = indiv_sample),
             color = "transparent",
             shape = 22, size = 6, alpha = 1) +
  geom_point(aes(x = dpi, y = indiv_sample, fill = as.character(pos_inf),
                 color = pos_inf_pred_correct,
                 group = indiv_sample),
             shape = 22, size = 4) +
  
  # Set colors
  scale_fill_manual(values = c("0" = subset(traj_palette, assay == "culture" & pos_neg == "neg")$hex,
                               "1" = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex,
                               "Edge" = "transparent",
                               "Transition" = "transparent",
                               "Blip" = "#BA141D", #"blue"
                               "Other" = "transparent")) +  
  scale_color_manual(values = c("Correct" = "black",
                                "Incorrect" = "transparent")) +  
  # Other plot features
  ggtitle("Data Blip") +
  labs(y = "individual", x = "", 
       fill = "Pos/Neg") +
  
  theme(legend.position = "none", 
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill='white'),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks= element_blank(),
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
        
        strip.text.y = element_blank(),
        plot.margin = margin(2, 0, 0, 2)) + 
  guides(color = "none") +
  scale_x_continuous(breaks = seq(0, 25, by = 3), limits = c(-0.01, 22),
                     expand = c(0, 0)); figSX_blip


nrow(subset(dat.blip, pos_inf_pred_class == "FN" & error_type == "Blip")) 
nrow(subset(dat.blip, pos_inf_pred_class == "FP" & error_type == "Blip"))

dat.blip$culture_assay[dat.blip$pos_inf_pred_class == "FN" & dat.blip$error_type == "Blip"]
dat.blip$culture_assay[dat.blip$pos_inf_pred_class == "FP" & dat.blip$error_type == "Blip"]


### Supplemental trajectories --------------------------------------------------

dat.blip.traj <- dat.blip
dat.blip.traj$val_total[dat.blip.traj$val_total == -9] <- -0.5
dat.blip.traj$sample_rep[dat.blip.traj$sample_rep == "Nasopharyngeal Swab"] <- "Nasoph. Swab"
dat.blip.traj$sample_rep[dat.blip.traj$sample_rep == "Nasal Cavity Swab"] <- "Nasal Swab"
dat.blip.traj$sample_rep[dat.blip.traj$sample_rep == "Oropharyngeal Swab"] <- "Oroph. Swab"


dat.blip.traj$panel_label <-paste0(dat.blip.traj$indiv, "\n", dat.blip.traj$sample_rep)

dat.blip.traj <- subset(dat.blip.traj, article != "Johnston et al. 2020")

figSX_blip_traj <- ggplot(dat.blip.traj) + 
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dotted", color = "grey") +
  

  # tRNA & sgRNA trajectories
  geom_line(aes(x = dpi, y = val_total),alpha = 0.3) +
  geom_point(data = subset(dat.blip.traj, error_type == "Blip"),
             aes(x = dpi, y = val_total, fill = as.character(error_type),
                 group = indiv_sample),
             color = "transparent",
             shape = 22, size = 5, alpha = 1) +
  geom_point(aes(x = dpi, y = val_total, fill = as.character(pos_inf),
                 color = pos_inf_pred_correct,
                 #size = pos_inf,
                 group = indiv_sample),
             shape = 22, size = 3) +
  geom_text(aes(label = panel_label, x = 15, y = 10), size = 2.2) +
  
  # Set colors
  scale_fill_manual(values = c("0" = subset(traj_palette, assay == "culture" & pos_neg == "neg")$hex,
                               "1" = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex,
                               "Edge" = "transparent",
                               "Transition" = "transparent",
                               "Blip" = "#BA141D", #"blue"
                               "Other" = "transparent")) +  
  scale_color_manual(values = c("Correct" = "black",
                                "Incorrect" = "transparent")) +  
  # Other plot features
  facet_wrap(indiv_sample ~., ncol = 9) +
  ggtitle("Data Blip") +
  labs(y = "individual", x = "day post infection", 
       fill = "Pos/Neg", tag = "A") +
  theme(legend.position = "none", 
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill='white'),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),

        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        strip.text = element_blank(),
        plot.margin = margin(2, 0, 0, 2)) + 
  guides(color = "none") +
  scale_y_continuous(breaks = seq(0, 11, by = 2), limits = c(-1, 12),
                     expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 25, by = 6), limits = c(-0.01, 21),
                     expand = c(0, 0)); figSX_blip_traj


## Transition ------------------------------------------------------------------

dat.trans <- subset(dat.ni, indiv_sample %in% trans_list)
figSX_trans <- ggplot(dat.trans) + 
  
  # tRNA & sgRNA trajectories
  geom_line(aes(x = dpi, y = indiv_sample),
            color = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex, 
            alpha = 0.3) +
  geom_point(data = subset(dat.trans, error_type == "Transition"),
             aes(x = dpi, y = indiv_sample, fill = as.character(error_type),
                 group = indiv_sample),
             color = "transparent",
             shape = 22, size = 6, alpha =1) +
  geom_point(aes(x = dpi, y = indiv_sample, fill = as.character(pos_inf),
                 color = pos_inf_pred_correct,
                 #size = pos_inf,
                 group = indiv_sample),
             shape = 22, size = 4) +
  
  # Set colors
  scale_fill_manual(values = c("0" = subset(traj_palette, assay == "culture" & pos_neg == "neg")$hex,
                               "1" = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex,
                               "Edge" = "transparent",
                               "Transition" = "#BA141D", #"green", #"#99D977",
                               "Blip" = "transparent",
                               "Other" = "transparent")) +  
  scale_color_manual(values = c("Correct" = "black",
                                "Incorrect" = "transparent")) +  
  # Other plot features
  labs(y = "individual", x = "", 
       fill = "Pos/Neg") +
  ggtitle("Transition") +
  theme(legend.position = "none", 
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill='white'),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        axis.text = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y= element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
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
        
        strip.text.y = element_blank(),
        plot.margin = margin(2, 0, 0, 2)) + 
  guides(color = "none") +
  scale_x_continuous(breaks = seq(0, 25, by = 3), limits = c(-0.01, 22),
                     expand = c(0, 0)); figSX_trans



## Prediction Blip -------------------------------------------------------------

dat.mid <- subset(dat.ni, indiv_sample %in% mid_list)
figSX_mid <- ggplot(dat.mid) + 
  
  # tRNA & sgRNA trajectories
  geom_line(aes(x = dpi, y = indiv_sample),
            color = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex, 
            alpha = 0.3) +
  geom_point(data = subset(dat.mid, error_type == "Mid-String"),
             aes(x = dpi, y = indiv_sample, fill = as.character(error_type),
                 group = indiv_sample),
             color = "transparent",
             shape = 22, size = 6, alpha = 1) +
  geom_point(aes(x = dpi, y = indiv_sample, fill = as.character(pos_inf),
                 color = pos_inf_pred_correct,
                 group = indiv_sample),
             shape = 22, size = 4) +
  
  # Set colors
  scale_fill_manual(values = c("0" = subset(traj_palette, assay == "culture" & pos_neg == "neg")$hex,
                               "1" = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex,
                               "Mid-String" =  "#BA141D")) + #"#B46FD1")) +  
  scale_color_manual(values = c("Correct" = "black",
                                "Incorrect" = "transparent")) +  
  # Other plot features
  ggtitle("Prediction Blip") +
  labs(y = "individual", x = "", 
       fill = "Pos/Neg") +
  
  theme(legend.position = "none", 
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill='white'),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        axis.text = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y= element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
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
        
        strip.text.y = element_blank(),
        plot.margin = margin(2, 0, 0, 2)) + 
  guides(color = "none") +
  scale_x_continuous(breaks = seq(0, 25, by = 3), limits = c(-0.01, 22),
                     expand = c(0, 0)); figSX_mid



nrow(subset(dat.mid, pos_inf_pred_class == "FN" & error_type == "Mid-String")) 
nrow(subset(dat.mid, pos_inf_pred_class == "FP" & error_type == "Mid-String"))


dat.mid$culture_assay[dat.mid$pos_inf_pred_class == "FN" & dat.mid$error_type == "Mid-String"]
dat.mid$culture_assay[dat.mid$pos_inf_pred_class == "FP" & dat.mid$error_type == "Mid-String"]



### Supplemental trajectories --------------------------------------------------

dat.mid.traj <- dat.mid
dat.mid.traj$val_total[dat.mid.traj$val_total == -9] <- -0.5
dat.mid.traj$sample_rep[dat.mid.traj$sample_rep == "Nasopharyngeal Swab"] <- "Nasoph. Swab"
dat.mid.traj$sample_rep[dat.mid.traj$sample_rep == "Nasal Cavity Swab"] <- "Nasal Swab"
dat.mid.traj$sample_rep[dat.mid.traj$sample_rep == "Oropharyngeal Swab"] <- "Oroph. Swab"


dat.mid.traj$panel_label <-paste0(dat.mid.traj$indiv, "\n", dat.mid.traj$sample_rep)
dat.mid.traj <- subset(dat.mid.traj, article != "Johnston et al. 2020")


figSX_mid_traj <- ggplot(dat.mid.traj) + 
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dotted", color = "grey") +
  
  
  # tRNA & sgRNA trajectories
  geom_line(aes(x = dpi, y = val_total),alpha = 0.3) +
  geom_point(data = subset(dat.mid.traj, error_type == "Mid-String"),
             aes(x = dpi, y = val_total, fill = as.character(error_type),
                 group = indiv_sample),
             color = "#BA141D",
             shape = 22, size = 4.5, alpha = 1) +
  geom_point(aes(x = dpi, y = val_total, fill = as.character(pos_inf),
                 color = pos_inf_pred_correct,
                 group = indiv_sample),
             shape = 22, size = 3) +
  geom_text(aes(label = panel_label, x = 15, y = 10), size = 2.2) +
  # Set colors
  scale_fill_manual(values = c("0" = subset(traj_palette, assay == "culture" & pos_neg == "neg")$hex,
                               "1" = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex,
                               "Edge" = "transparent",
                               "Transition" = "transparent",
                               "Mid-String" = "#BA141D", #"blue"
                               "Other" = "transparent")) +  
  scale_color_manual(values = c("Correct" = "black",
                                "Incorrect" = "transparent")) +  
  # Other plot features
  facet_wrap(indiv_sample ~., ncol = 9) +
  ggtitle("Prediction Blip") +
  labs(y = "individual", x = "day post infection", 
       fill = "Pos/Neg", tag = "B") +
  theme(legend.position = "none", 
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill='white'),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        strip.text = element_blank(),
        plot.margin = margin(2, 0, 0, 2)) + 
  guides(color = "none") +
  scale_y_continuous(breaks = seq(0, 11, by = 2), limits = c(-1, 12),
                     expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 25, by = 6), limits = c(-0.01, 21),
                     expand = c(0, 0)); figSX_mid_traj



## None ------------------------------------------------------------------------

figSX_none <- ggplot(dat.none) + 
  
  # tRNA & sgRNA trajectories
  geom_line(aes(x = dpi, y = indiv_sample),
            color = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex, 
            alpha = 0.3) +
  geom_point(aes(x = dpi, y = indiv_sample, fill = as.character(pos_inf),
                 color = pos_inf_pred_correct,
                 #size = pos_inf,
                 group = indiv_sample),
             shape = 22, size = 2.5) +
  
  # Set colors
  scale_fill_manual(values = c("0" = subset(traj_palette, assay == "culture" & pos_neg == "neg")$hex,
                               "1" = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex)) +  
  scale_color_manual(values = c("Correct" = "black",
                                "Incorrect" = "transparent")) +  
  # Other plot features
  labs(y = "individual", x = "", 
       fill = "Pos/Neg") +
  ggtitle("No Errors") +
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
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y= element_blank(),
        text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
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
        
        strip.text.y = element_blank(),
        plot.margin = margin(2, 0, 0, 2)) + 
  guides(color = "none") +
  scale_x_continuous(breaks = seq(0, 25, by = 3), limits = c(-0.01, 22),
                     expand = c(0, 0)); figSX_none


# TotRNA values & dpis for incorrect samples -----------------------------------


dat.wrong <- subset(dat.ni, !is.na(val_total))
dat.wrong$error_type[is.na(dat.wrong$error_type)] <- "No Error"

dat.inv$error_type <- "All Invasive"
dat.inv <- subset(dat.inv, !is.na(val_total))
dat.wrong <- rbind(dat.wrong, dat.inv)

dat.wrong$val_total[dat.wrong$val_total == -9] <- -1

dat.wrong$error_type[dat.wrong$error_type == "Blip"] <- "Data Blip"
dat.wrong$error_type[dat.wrong$error_type == "Mid-String"] <- "Prediction Blip"

dat.wrong$val_total[dat.wrong$val_total == -1] <- runif(length(dat.wrong$val_total[dat.wrong$val_total == -1]),
                                                        -1.5, -0.5)

dat.wrong$pos_inf_pred_class[dat.wrong$pos_inf_pred_class == "TP"] <- "True Positive"
dat.wrong$pos_inf_pred_class[dat.wrong$pos_inf_pred_class == "TN"] <- "True Negative"
dat.wrong$pos_inf_pred_class[dat.wrong$pos_inf_pred_class == "FN"] <- "False Negative"
dat.wrong$pos_inf_pred_class[dat.wrong$pos_inf_pred_class == "FP"] <- "False Positive"




figSXA <- ggplot(dat.wrong) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey90",
            alpha = 0.5, color = "black") +
  geom_point(aes(x = dpi,
                 y = as.numeric(val_total), 
                 fill = pos_inf_pred_class,
                 color = pos_inf_pred_class),
             alpha = 0.7,
             shape = 22, 
             size = 2) +
  scale_fill_manual(values = c("False Positive" = subset(traj_palette, assay == "culture" & pos_neg == "neg")$hex,
                               "False Negative" = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex,
                               "True Positive" = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex,
                               "True Negative" = subset(traj_palette, assay == "culture" & pos_neg == "neg")$hex),
                    breaks = c("True Positive", "False Negative", 
                               "True Negative", "False Positive")) + 
  scale_color_manual(values = c("False Positive" = "#BA141D",
                                "False Negative" = "#BA141D",
                                "True Positive" = "black",
                                "True Negative" = "black")) + 
  facet_wrap(.~ factor(error_type, 
                       levels = c("Edge", "Transition", "Data Blip", "Prediction Blip", 
                                  "All Invasive", "No Error")),
             nrow = 1) +
  labs(x = "day post infection",
       y = "log10 totRNA copies / sample", fill = "Error Type",
       tag = "B") +
  scale_x_continuous(breaks = seq(0, 23, 6), limits = c(0, 22)) +
  scale_y_continuous(breaks = c(seq(0, 12, by = 2)), limits = c(-2.3, 12),
                     labels = c(seq(0, 12, by = 2)), 
                     expand = c(0, 0)) +
  guides(color = "none", 
         fill = guide_legend(override.aes = list(color = c("black", "#BA141D",
                                                           "black", "#BA141D"),
                                                 size = 3))) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill='white'),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        text = element_text(size = 14),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        strip.background  = element_blank(),
        strip.text.x = element_text(face = "bold", size = 11),
        plot.margin = margin(2, 4, 0, 2)); figSXA

# Combine into Figure 6 --------------------------------------------------------

figSX_smallpanels <- figSX_blip + figSX_mid + 
  plot_layout(nrow = 2, heights = c(length(blip_list), length(mid_list))); figSX_smallpanels

figSX <- figSX_edge + figSX_trans + figSX_smallpanels + 
  plot_layout(ncol = 3); figSX

x_label <- ggplot() +
  annotate("text", x = 1, y = 1, size = 4, 
           label = "day post infection") + 
  theme_void(); x_label

figSX <- ggarrange(figSX, x_label, nrow = 2, heights = c(1, 0.03)) +
  bgcolor("white")+
  border(color = "white"); figSX

fig5 <- ggarrange(figSX, figSXA, nrow = 2, heights = c(1, 0.35)); fig5

ggsave("./outputs/figures/fig6-culture-error-analysis.tiff", 
       fig5, width = 11, height = 11, units = "in")


# Combine into Supplemental Figure 17 ------------------------------------------

figSX <- figSX_blip_traj + figSX_mid_traj + plot_layout(nrow = 2)

ggsave("./outputs/figures/figS17-culture-blip-trajectories.tiff", 
       figSX, width = 12, height = 6.2, units = "in")
