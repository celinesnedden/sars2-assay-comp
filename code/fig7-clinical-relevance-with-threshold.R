# Generates figure 7, including analysis of clinical relevance of predictions

# Prep environment -------------------------------------------------------------

# Install & load packages
req_pkgs <- c("wesanderson", "ggplot2", "cmdstanr", "tidyverse", "ggridges",
              "ggpubr", "patchwork", "ggExtra", "ggbeeswarm")
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

# Best model
dat.best <- read.csv("./data/pred-culture-data.csv")
dat.best <- subset(dat.best, !is.na(pos_total))
dat.best$model <- "Best"

# Simple model
dat.simple <- read.csv("./data/pred-culture-data-simple.csv")
dat.simple <- subset(dat.simple, !is.na(pos_total))
dat.simple$model <- "Simple"

# Combine for classifying
dat <- rbind(dat.best, dat.simple)


# Prep data --------------------------------------------------------------------

## Add sample indicator for facetting -------------------------------------------

## Combine indiv names with sample rep to get the right connecting lines
dat$indiv_sample <- paste0(dat$indiv, "_", dat$sample_rep, 
                           #"_", dat$tg_idx)
                           "_", dat$tg_total)



## Subset to non-invasives only ------------------------------------------------

dat.ni <- subset(dat, sample_type == "Non-invasive" & !is.na(pos_inf) &
                   sample_grp %notin% c("Anal / Rectal Sample", "Ocular Sample"))


## Subset to individuals that ever test positive -------------------------------

for (indiv.ii in unique(dat.ni$indiv_sample)) {
  dat.sub <- subset(dat.ni, indiv_sample == indiv.ii)
  if (1 %notin% dat.sub$pos_inf) {
    dat.ni <- subset(dat.ni, indiv_sample != indiv.ii)
  }
}



# Find ALL relevant times ------------------------------------------------------

find_times <- function(dat.ni, threshold, return = "one per") {
  
  if (threshold != 0){
    dat.ni$pos_inf_pred[dat.ni$pos_inf_prob_pred < threshold] <- 0
    dat.ni$pos_inf_pred[dat.ni$pos_inf_prob_pred >= threshold] <- 1
  }
  else {
    dat.ni$pos_inf_pred[dat.ni$pos_inf_prob_pred <= threshold] <- 0
    dat.ni$pos_inf_pred[dat.ni$pos_inf_prob_pred > threshold] <- 1
  }
  
  # Set up df to return
  df.indiv.traj <- data.frame(indiv_sample = character(),
                              model = character(), # Best vs Simple model
                              
                              # Basic known times
                              dpi_first_pos_test = numeric(), # day of first positive test by any assay
                              last_pos_true_dpfp = numeric(), # last culture positive
                              next_neg_true_dpfp = numeric(), # next culture negative after last positive
                              end_inf_period_dpfp = numeric(), # the end of the infectious period (midpoint of previous 2)
                              
                              # First ever predicted negatives
                              first_ever_neg_pred_dpfp = numeric(),
                              first_ever_neg_pred_prob = numeric(),
                              first_ever_neg_true_dpfp = numeric(),
                              first_ever_neg_true_prob = numeric(),
                              
                              # Known and predicted consecutive negative times, plus probabilities
                              first_cons_neg_pred_dpfp = numeric(),
                              first_cons_neg_pred_prob = numeric(),
                              first_cons_neg_true_dpfp = numeric(),
                              first_cons_neg_true_prob = numeric(),
                              second_cons_neg_pred_dpfp = numeric(),
                              second_cons_neg_pred_prob = numeric(),
                              second_cons_neg_true_dpfp = numeric(),
                              second_cons_neg_true_prob = numeric(),
                              
                              
                              # Days still infectious
                              first_ever_neg_inf_days = numeric(), # infectious days using the 1st ever culture negative
                              second_cons_neg_inf_days = numeric(), # infectious days using the 2nd consecutive negative
                              
                              # Excess days isolated
                              first_ever_neg_excess_isol_days = numeric(),
                              second_cons_neg_excess_isol_days = numeric(),
                              
                              # Rebound class
                              rebound_class = numeric()
  )
  
  df.indiv.all <- dat.ni[1, ]
  df.indiv.all$days_post_first_pos <- NA
  
  for (indiv.ii in unique(dat.ni$indiv_sample)) {
    
    dpi_first_pos_test <- NA # day of first positive test by any assay
    last_pos_true_dpfp <- NA # last culture positive
    next_neg_true_dpfp <- NA # next culture negative after last positive
    end_inf_period_dpfp <- NA # the end of the infectious period (midpoint of previous 2)
    
    dat.indiv.sub <- subset(dat.ni, indiv_sample == indiv.ii)
    
    # Calculate true event times (not model dependent)
    dat.true.days <- subset(dat.indiv.sub, model == "Best")
    
    # First positive by any test
    first_pos <- min(as.numeric(dat.indiv.sub$dpi[dat.indiv.sub$pos_total == 1]))
    first_inf_pos <- min(as.numeric(dat.indiv.sub$dpi[dat.indiv.sub$pos_inf == 1]))
    if (first_pos > first_inf_pos) {
      dpi_first_pos_test <- first_inf_pos
    } else {
      dpi_first_pos_test <- first_pos
    }
    
    # Recalibrate all other dpis to this time 
    dat.true.days$days_post_first_pos <- dat.true.days$dpi - dpi_first_pos_test
    dat.indiv.sub$days_post_first_pos <- dat.indiv.sub$dpi - dpi_first_pos_test
    
    # Remove all days before first positive by any test
    dat.true.days <- subset(dat.true.days, days_post_first_pos >= 0)
    dat.indiv.sub <- subset(dat.indiv.sub, days_post_first_pos >= 0)
    
    # Add this to data with all sample times
    df.indiv.all <- rbind(df.indiv.all, dat.indiv.sub)
    
    # First ever culture negative time, besides day 0
    first_ever_neg_true_dpfp <- min(as.numeric(dat.true.days$days_post_first_pos[dat.true.days$pos_inf == 0 &
                                                                                   dat.true.days$days_post_first_pos > 0]))
    if (is.na(first_ever_neg_true_dpfp) | first_ever_neg_true_dpfp == Inf) {
      first_ever_neg_true_dpfp <- "Not Observed"
    }
    
    # Last culture positive time 
    last_pos_true_dpfp <- max(as.numeric(dat.true.days$days_post_first_pos[dat.true.days$pos_inf == 1]))
    
    if (is.na(last_pos_true_dpfp)) {
      print("NA LAST POSITIVE TIME")
    }
    if (last_pos_true_dpfp == Inf) {
      print("INF LAST POSITIVE TIME")
    }
    
    # Next culture negative time after last culture positive
    next_neg_true_dpfp <- min(as.numeric(dat.true.days$days_post_first_pos[
      dat.true.days$pos_inf == 0 & dat.true.days$days_post_first_pos > last_pos_true_dpfp]))
    if (is.na(next_neg_true_dpfp) | next_neg_true_dpfp == Inf) {
      next_neg_true_dpfp <- "Not Observed"
    }
    
    # Calculate end of the infectious period as midpoint between the previous two days
    if (next_neg_true_dpfp != "Not Observed") {
      end_inf_period_dpfp <- ceiling(last_pos_true_dpfp + (next_neg_true_dpfp - last_pos_true_dpfp)/2)
    } else {
      end_inf_period_dpfp <- ceiling(last_pos_true_dpfp + (10 - last_pos_true_dpfp)/2)
    }
    
    for (model.ii in c("Best", "Simple")) {
      dat.model.sub <- subset(dat.indiv.sub, model == model.ii)
      dat.model.sub <- dat.model.sub[order(as.numeric(dat.model.sub$days_post_first_pos)), ]
      print(dat.model.sub$days_post_first_pos)
      
      # Reset everything so it doesn't carry over from previous
      first_ever_neg_pred_dpfp <- NA
      first_ever_neg_pred_prob <- NA
      first_ever_neg_true_prob <- NA
      second_cons_neg_true_dpfp <- NA
      second_cons_neg_true_prob <- NA
      first_cons_neg_true_dpfp <- NA
      first_cons_neg_true_prob <- NA
      second_cons_neg_pred_dpfp <- NA
      second_cons_neg_pred_prob <- NA
      first_cons_neg_pred_dpfp <- NA
      first_cons_neg_pred_prob <- NA
      second_cons_neg_inf_days <- NA
      second_cons_neg_excess_isol_days <- NA
      first_ever_neg_inf_days <- NA
      first_ever_neg_excess_isol_days <- NA
      
      # Must have at least two samples beyond the first positive day to analyze
      if (nrow(dat.model.sub) >= 3) {
        
        # Initialize some columns to get store info from previous sample days
        dat.model.sub$prev_inf_pred <- NA
        dat.model.sub$prev_inf_pred_prob <- NA
        dat.model.sub$prev_pos_inf <- NA
        dat.model.sub$prev_dpfp <- NA
        
        # Start with the second row to get quantities on previous days
        for (row_num in 2:nrow(dat.model.sub)) {
          dat.model.sub$prev_inf_pred[row_num] <- dat.model.sub$pos_inf_pred[row_num - 1] 
          dat.model.sub$prev_inf_pred_prob[row_num] <- dat.model.sub$pos_inf_prob_pred[row_num - 1] 
          dat.model.sub$prev_pos_inf[row_num] <- dat.model.sub$pos_inf[row_num - 1]
          dat.model.sub$prev_dpfp[row_num] <- dat.model.sub$days_post_first_pos[row_num - 1]
        }
        
        # First predicted negative after the first positive day
        first_ever_neg_pred_dpfp <- min(dat.model.sub$days_post_first_pos[dat.model.sub$pos_inf_pred == 0 &
                                                                            dat.model.sub$days_post_first_pos > 0])
        if (is.na(first_ever_neg_pred_dpfp) | first_ever_neg_pred_dpfp == Inf) {
          first_ever_neg_pred_dpfp <- "Not Observed"
        }
        if (first_ever_neg_pred_dpfp == "Not Observed") {
          first_ever_neg_pred_prob <- "Not Observed"
        } else {
          first_ever_neg_pred_prob <- dat.model.sub$pos_inf_prob_pred[dat.model.sub$days_post_first_pos == first_ever_neg_pred_dpfp]
        }
        
        # Probability of true first negative, if ever happens
        if (first_ever_neg_true_dpfp == "Not Observed") {
          first_ever_neg_true_prob <- "Not Observed"
        } else {
          first_ever_neg_true_prob <- dat.model.sub$pos_inf_prob_pred[dat.model.sub$days_post_first_pos == first_ever_neg_true_dpfp]
        }
        
        # Find consecutive negative times (predicted and true)
        dat.model.sub.true <- subset(dat.model.sub, pos_inf == 0 & prev_pos_inf == 0)
        dat.model.sub.pred <- subset(dat.model.sub, pos_inf_pred == 0 & prev_inf_pred == 0)
        
        # If it ever truly happens, get the relevant probabilties / times
        if (nrow(dat.model.sub.true) > 0){
          
          # Take the first consecutive negative time  
          dat.model.sub.true <- subset(dat.model.sub.true, days_post_first_pos == min(as.numeric(dat.model.sub.true$days_post_first_pos)))
          
          second_cons_neg_true_dpfp <- unique(dat.model.sub.true$days_post_first_pos)
          second_cons_neg_true_prob <- unique(dat.model.sub.true$pos_inf_prob_pred)
          first_cons_neg_true_dpfp <- unique(dat.model.sub.true$prev_dpfp)
          first_cons_neg_true_prob <- unique(dat.model.sub.true$prev_inf_pred_prob)
        } else {
          second_cons_neg_true_dpfp <- "Not Observed"
          second_cons_neg_true_prob <- "Not Observed"
          first_cons_neg_true_dpfp <- "Not Observed"
          first_cons_neg_true_prob <- "Not Observed"
        }
        
        # If it's ever predicted, get the relevant probabilties / times
        if (nrow(dat.model.sub.pred) >= 1) {
          
          # Take the first consecutive negative time  
          dat.model.sub.pred <- subset(dat.model.sub.pred, days_post_first_pos == min(as.numeric(dat.model.sub.pred$days_post_first_pos)))
          second_cons_neg_pred_dpfp <- unique(dat.model.sub.pred$days_post_first_pos)
          second_cons_neg_pred_prob <- unique(dat.model.sub.pred$pos_inf_prob_pred)
          first_cons_neg_pred_dpfp <- unique(dat.model.sub.pred$prev_dpfp)
          first_cons_neg_pred_prob <- unique(dat.model.sub.pred$prev_inf_pred_prob)
          
          # Calculate the numbers of infectious days & excess isolation days
          second_cons_neg_inf_days <- end_inf_period_dpfp - second_cons_neg_pred_dpfp
          second_cons_neg_excess_isol_days <- second_cons_neg_pred_dpfp - end_inf_period_dpfp
          if (second_cons_neg_inf_days < 0){
            second_cons_neg_inf_days <- 0
          }
          else if (second_cons_neg_excess_isol_days < 0){
            second_cons_neg_excess_isol_days <- 0
          }
          
        } else {
          second_cons_neg_pred_dpfp <- "Not Observed"
          second_cons_neg_pred_prob <- "Not Observed"
          first_cons_neg_pred_dpfp <- "Not Observed"
          first_cons_neg_pred_prob <- "Not Observed"
          second_cons_neg_inf_days <- "Not Observed"
          second_cons_neg_excess_isol_days <- "Not Observed"
        }
        
        # Calculate the numbers of infectious days & excess infectious days for 1st ever negative
        if (first_ever_neg_pred_dpfp == Inf | first_ever_neg_pred_dpfp == "Not Observed") {
          first_ever_neg_inf_days <- "Negative Never Predicted"
          first_ever_neg_excess_isol_days <- "Negative Never Predicted"
        } else {
          first_ever_neg_inf_days <- end_inf_period_dpfp - first_ever_neg_pred_dpfp 
          first_ever_neg_excess_isol_days <- first_ever_neg_pred_dpfp - end_inf_period_dpfp
        }
        # Can't have negatives -- set them to zero
        if (first_ever_neg_inf_days < 0){
          first_ever_neg_inf_days <- 0
        }
        if (first_ever_neg_excess_isol_days < 0){
          first_ever_neg_excess_isol_days <- 0
        }
        
        # Figure out if theres evidence of an early / late rebound 
        indiv_pos_string <- paste0(dat.model.sub$pos_inf, collapse = "")
        if (str_detect(indiv_pos_string, "1001|101|10001")) {
          if (last_pos_true_dpfp < 5) {
            rebound_class <- "Early Rebound"
          }
          else if (last_pos_true_dpfp >= 5) {
            rebound_class <- "Late Rebound"
          }
          
        }
        else {
          rebound_class <- "None"
        }
        
        
        next_row <- data.frame(indiv_sample = indiv.ii,
                               model = model.ii, # Best vs Simple model
                               
                               # Basic known times
                               dpi_first_pos_test = dpi_first_pos_test, # day of first positive test by any assay
                               last_pos_true_dpfp = last_pos_true_dpfp, # last culture positive
                               next_neg_true_dpfp = next_neg_true_dpfp, # next culture negative after last positive
                               end_inf_period_dpfp = end_inf_period_dpfp, # the end of the infectious period (midpoint of previous 2)
                               
                               # First ever predicted negatives
                               first_ever_neg_pred_dpfp = first_ever_neg_pred_dpfp ,
                               first_ever_neg_pred_prob = first_ever_neg_pred_prob,
                               first_ever_neg_true_dpfp = first_ever_neg_true_dpfp,
                               first_ever_neg_true_prob = first_ever_neg_true_prob,
                               
                               # Known and predicted consecutive negative times, plus probabilities
                               first_cons_neg_pred_dpfp = first_cons_neg_pred_dpfp,
                               first_cons_neg_pred_prob = first_cons_neg_pred_prob,
                               first_cons_neg_true_dpfp = first_cons_neg_true_dpfp,
                               first_cons_neg_true_prob = first_cons_neg_true_prob,
                               second_cons_neg_pred_dpfp = second_cons_neg_pred_dpfp,
                               second_cons_neg_pred_prob = second_cons_neg_pred_prob,
                               second_cons_neg_true_dpfp = second_cons_neg_true_dpfp,
                               second_cons_neg_true_prob = second_cons_neg_true_prob,
                               
                               
                               # Days still infectious
                               first_ever_neg_inf_days = first_ever_neg_inf_days, # infectious days using the 1st ever culture negative
                               second_cons_neg_inf_days = second_cons_neg_inf_days, # infectious days using the 2nd consecutive negative
                               
                               # Excess days isolated
                               first_ever_neg_excess_isol_days = first_ever_neg_excess_isol_days,
                               second_cons_neg_excess_isol_days = second_cons_neg_excess_isol_days,
                               
                               # Rebound class
                               rebound_class = rebound_class
        )
        
        df.indiv.traj <- rbind(df.indiv.traj, next_row)
      }
    }
  }
  
  if (return == "one per") {
    return(df.indiv.traj)
  }
  else if (return == "all days") {
    df.indiv.all <- df.indiv.all[-1, ]
    return(df.indiv.all)
  }
}


df.indiv.traj <- find_times(dat.ni, 0.5)
df.indiv.all <- find_times(dat.ni, 0.5, "all days")
df.indiv.traj.0 <- find_times(dat.ni, 0.01)
df.indiv.traj.1 <- find_times(dat.ni, 0.1)
df.indiv.traj.2 <- find_times(dat.ni, 0.2)
df.indiv.traj.3 <- find_times(dat.ni, 0.3)
df.indiv.traj.4 <- find_times(dat.ni, 0.4)


# Remove individuals not classified by either model ----------------------------

# Find them for each and their overlap
best_unclass_indivs <- unique(subset(df.indiv.traj, model == "Best" & second_cons_neg_pred_dpfp == "Not Observed")$indiv_sample)
simple_unclass_indivs <- unique(subset(df.indiv.traj, model == "Simple" & second_cons_neg_pred_dpfp == "Not Observed")$indiv_sample)

remove_indivs <- unique(df.indiv.traj$indiv_sample[df.indiv.traj$indiv_sample %in% best_unclass_indivs &
                                                     df.indiv.traj$indiv_sample %in% simple_unclass_indivs])

# Remove them
df.indiv.traj <- subset(df.indiv.traj, indiv_sample %notin% remove_indivs)


# Change their isolation end times to be 10 days & adjust infectious / excess isolation times
df.indiv.traj$unclassified <- "No"
df.indiv.traj$unclassified[df.indiv.traj$second_cons_neg_pred_dpfp == "Not Observed"] <- "Yes"


# Conservatively have the unclassified individuals stop isolating at day 10
#   and take the largest observed end of the infectious period as their
#   infectious period end time
df.indiv.traj$second_cons_neg_pred_dpfp[df.indiv.traj$unclassified == "Yes"] <- 10
df.indiv.traj$second_cons_neg_excess_isol_days[df.indiv.traj$unclassified == "Yes"] <- 10 - max(df.indiv.traj$end_inf_period_dpfp, na.rm = TRUE)
df.indiv.traj$second_cons_neg_inf_days[df.indiv.traj$unclassified == "Yes"] <-  0





# Subset to model specific dfs --------------------------------------------------

df.indiv.traj.best <- subset(df.indiv.traj, model == "Best")
df.indiv.traj.simple <- subset(df.indiv.traj, model == "Simple")

# Create ten / five day 
df.indiv.traj.ten <- df.indiv.traj.best
df.indiv.traj.ten$isol_days <- 10
df.indiv.traj.ten$excess_isol_days <- 10 - df.indiv.traj.ten$end_inf_period_dpfp 
df.indiv.traj.ten$excess_isol_days[df.indiv.traj.ten$excess_isol_days < 0] <- 0
df.indiv.traj.ten$inf_days <- df.indiv.traj.ten$end_inf_period_dpfp - 10
df.indiv.traj.ten$inf_days[df.indiv.traj.ten$inf_days < 0] <- 0

df.indiv.traj.ten$model <- "Ten Day"


df.indiv.traj.five <- df.indiv.traj.best
df.indiv.traj.five$isol_days <- 5
df.indiv.traj.five$excess_isol_days <- 5 - df.indiv.traj.five$end_inf_period_dpfp 
df.indiv.traj.five$excess_isol_days[df.indiv.traj.five$excess_isol_days < 0] <- 0
df.indiv.traj.five$model <- "Five Day"
df.indiv.traj.five$inf_days <- df.indiv.traj.five$end_inf_period_dpfp  - 5
df.indiv.traj.five$inf_days[df.indiv.traj.five$inf_days < 0] <- 0

df.indiv.traj.threshold <- rbind(df.indiv.traj.ten, df.indiv.traj.five)


# Perfect model
df.indiv.traj.perfect <- df.indiv.traj.best 
df.indiv.traj.perfect$second_cons_neg_true_dpfp[df.indiv.traj.perfect$second_cons_neg_true_dpfp == "Not Observed"] <- NA
df.indiv.traj.perfect$excess_isol_days <- as.numeric(df.indiv.traj.perfect$second_cons_neg_true_dpfp) - as.numeric(df.indiv.traj.perfect$end_inf_period_dpfp)
df.indiv.traj.perfect$excess_isol_days[df.indiv.traj.perfect$excess_isol_days < 0] <- 0
df.indiv.traj.perfect$model <- "Perfect"
df.indiv.traj.perfect$inf_days <- as.numeric(df.indiv.traj.perfect$end_inf_period_dpfp) - 
  as.numeric(df.indiv.traj.perfect$second_cons_neg_true_dpfp)
df.indiv.traj.perfect$inf_days[df.indiv.traj.perfect$inf_days < 0] <- 0



# FIG S18: Trajectories ---------------------------------------------------------

## S18A: Best ---------------------------------------------------

df.indiv.traj.best <- df.indiv.traj.best[order(as.numeric(df.indiv.traj.best$end_inf_period_dpfp),
                                               as.numeric(df.indiv.traj.best$second_cons_neg_pred_dpfp)), ]

df.indiv.traj.best$indiv_sample_factor <- factor(df.indiv.traj.best$indiv_sample, 
                                                 levels = df.indiv.traj.best$indiv_sample)

df.indiv.all <- subset(df.indiv.all, indiv_sample %in% unique(df.indiv.traj.best$indiv_sample))
df.indiv.all$indiv_sample_factor <- factor(df.indiv.all$indiv_sample, 
                                           levels = df.indiv.traj.best$indiv_sample)

df.indiv.all$pos_inf[df.indiv.all$pos_inf == 0] <- "Observed Culture Neg."
df.indiv.all$pos_inf[df.indiv.all$pos_inf == 1] <- "Observed Culture Pos."


# Get the # of rebound indivdiuals
table(df.indiv.traj.best$rebound_class)


fig.traj.best <- ggplot(df.indiv.traj.best) + 
  
  geom_vline(xintercept = 10, linewidth = 1, color = "deeppink4") +
  geom_vline(xintercept = 5,  linewidth = 1, color = "orchid4") +
  
  geom_segment(aes(x = 0, xend = as.numeric(second_cons_neg_pred_dpfp), 
                   y = indiv_sample_factor,
                   yend = indiv_sample_factor),
               color = "grey77",
               linewidth = 0.05,
               alpha = 1) +
  
  geom_segment(aes(x = 0, xend = as.numeric(end_inf_period_dpfp), 
                   y = indiv_sample_factor,
                   yend = indiv_sample_factor,
                   color = "Infectious Period"),
               linewidth = 0.5,
               alpha = 0.5) +
  
  geom_segment(data = subset(df.indiv.traj.best, rebound_class != "None"),
               aes(x = 0, xend = as.numeric(end_inf_period_dpfp), 
                   y = indiv_sample_factor,
                   yend = indiv_sample_factor),
               linetype = "dotted",
               color = "red",
               linewidth = 0.5,
               alpha = 0.8) +
  
  geom_point(aes(x = as.numeric(first_ever_neg_pred_dpfp), 
                 y = indiv_sample_factor,
                 color = "First Predicted Neg."),
             fill = "transparent",
             shape = 23,
             stroke = 1,
             size = 2.5,
             alpha = 1) +
  
  geom_point(aes(x =  as.numeric(second_cons_neg_pred_dpfp), 
                 y = indiv_sample_factor,
                 fill = "Second Predicted Neg."),
             shape = 23,
             size = 2.5,
             alpha = 1) +
  
  geom_point(aes(x =  as.numeric(second_cons_neg_true_dpfp), 
                 y = indiv_sample_factor, fill = "True Second Consecutive Neg."),
             shape = 21,
             stroke = 0,
             size = 1,
             alpha = 1) +
  
  geom_point(data = df.indiv.all,
             aes(x = days_post_first_pos,  
                 y = indiv_sample_factor,
                 fill = as.character(pos_inf)),
             alpha = 0.3, shape = 22) +
  
  
  geom_point(aes(x = as.numeric(last_pos_true_dpfp), y = indiv_sample_factor,
                 fill = "Last Culture Pos."),
             shape = 22,
             alpha = 1) +
  geom_point(aes(x = as.numeric(next_neg_true_dpfp), y = indiv_sample_factor,
                 fill = "Next Culture Neg."),
             shape = 22,
             alpha = 1) +
  
  # Set colors
  scale_fill_manual(values = c("Second Predicted Neg." = "#85B068",
                               "True Second Consecutive Neg." = "red",
                               "Last Culture Pos." = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex,
                               "Next Culture Neg." = "grey62",
                               "Observed Culture Pos." = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex,
                               "Observed Culture Neg." = "grey62")) +  
  scale_color_manual(values = c("Infectious Period" = "#E4BC11",
                                "First Predicted Neg." = "#85B068")) +  
  
  scale_x_continuous(breaks = seq(0, 14, by = 2), limits = c(0, 14)) +
  coord_cartesian(xlim = c(0, 14), clip = "off") +
  guides(fill = guide_legend(override.aes = list(shape = c(22, 22, 22, 22, 23),
                                                 size = c(1.5, 1.5, 1.5, 1.5, 2),
                                                 alpha = c(0.1, 0.1, 1, 1, 1)),
                             order = 1),
         color = guide_legend(override.aes = list(shape = c(23, NA)),
                              order = 2)) +
  
  # Other plot features
  labs(y = "individual", x = "days after first positive by any test", 
       fill = "Pos/Neg",
       tag = "B", title = "Best Model") +
  theme(legend.position = "none", 
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill='white', color = "transparent"),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.spacing.y = unit(-0.03, "cm"),
        legend.margin=margin(c(0,3,3,3)),
        axis.text = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 11),
        axis.ticks.y= element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        text = element_text(size = 11),
        strip.text = element_text(face = "italic"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.2, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background  = element_rect(colour = "white", fill = "white"),
        plot.margin = margin(2, 10, 2, 2)); fig.traj.best


## S18B: Simple -----------------------------------------------------------------

df.indiv.traj.simple$indiv_sample_factor <- factor(df.indiv.traj.simple$indiv_sample, 
                                                   levels = df.indiv.traj.best$indiv_sample_factor)

df.indiv.all <- subset(df.indiv.all, indiv_sample %in% unique(df.indiv.traj.simple$indiv_sample))
df.indiv.all$indiv_sample_factor <- factor(df.indiv.all$indiv_sample, 
                                           levels = df.indiv.traj.simple$indiv_sample)

fig.traj.simple <- ggplot(df.indiv.traj.simple) + 
  
  geom_vline(xintercept = 10, linewidth = 1, color = "deeppink4") +
  geom_vline(xintercept = 5,  linewidth = 1, color = "orchid4") +
  
  geom_segment(aes(x = 0, xend = as.numeric(second_cons_neg_pred_dpfp), 
                   y = indiv_sample_factor,
                   yend = indiv_sample_factor),
               color = "grey77",
               linewidth = 0.05,
               alpha = 1) +
  
  geom_segment(aes(x = 0, xend = as.numeric(end_inf_period_dpfp), 
                   y = indiv_sample_factor,
                   yend = indiv_sample_factor,
                   color = "Infectious Period"),
               linewidth = 0.5,
               alpha = 0.5) +
  
  
  geom_segment(data = subset(df.indiv.traj.simple, rebound_class != "None"),
               aes(x = 0, xend = as.numeric(end_inf_period_dpfp), 
                   y = indiv_sample_factor,
                   yend = indiv_sample_factor,
                   color = "Rebound Region"),
               linetype = "dotted",
               linewidth = 0.5,
               alpha = 0.8) +
  
  geom_point(aes(x = as.numeric(first_ever_neg_pred_dpfp), 
                 y = indiv_sample_factor,
                 color = "First Predicted Neg."),
             fill = "transparent",
             shape = 23,
             stroke = 1,
             size = 2.5,
             alpha = 1) +
  
  geom_point(aes(x =  as.numeric(second_cons_neg_pred_dpfp), 
                 y = indiv_sample_factor,
                 fill = "Second Predicted Neg."),
             shape = 23,
             size = 2.5,
             alpha = 1) +
  
  geom_point(aes(x =  as.numeric(second_cons_neg_true_dpfp), 
                 y = indiv_sample_factor, fill = "True Second Consecutive Neg."),
             shape = 21,
             stroke = 0,
             size = 1,
             alpha = 1) +
  
  geom_point(data = df.indiv.all,
             aes(x = days_post_first_pos,  
                 y = indiv_sample_factor,
                 fill = as.character(pos_inf)),
             alpha = 0.3, shape = 22) +
  
  
  geom_point(aes(x = as.numeric(last_pos_true_dpfp), y = indiv_sample_factor,
                 fill = "Last Culture Pos."),
             shape = 22,
             alpha = 1) +
  geom_point(aes(x = as.numeric(next_neg_true_dpfp), y = indiv_sample_factor,
                 fill = "Next Culture Neg."),
             shape = 22,
             alpha = 1) +
  
  # Set colors
  scale_fill_manual(values = c("Second Predicted Neg." = "#0F4C5E",
                               "True Second Consecutive Neg." = "red",
                               "Last Culture Pos." = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex,
                               "Next Culture Neg." = "grey62",
                               "Observed Culture Pos." = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex,
                               "Observed Culture Neg." = "grey62")) +  
  scale_color_manual(values = c("Infectious Period" = "#E4BC11",
                                "First Predicted Neg." = "#0F4C5E",
                                "Rebound Region" = zissou_pal[11])) +  
  
  scale_x_continuous(breaks = seq(0, 14, by = 2), limits = c(0, 14)) +
  coord_cartesian(xlim = c(0, 14), clip = "off") +
  guides(fill = guide_legend(override.aes = list(shape = c(22, 22, 22, 22, 23, 21),
                                                 size = c(1.5, 1.5, 1.5, 1.5, 2, 1),
                                                 alpha = c(0.1, 0.1, 1, 1, 1, 1)),
                             order = 1, nrow = 3),
         color = guide_legend(override.aes = list(shape = c(23, NA, NA)),
                              order = 2, nrow = 3)) +
  
  # Other plot features
  labs(y = "individual", x = "days after first positive by any test", 
       fill = "Pos/Neg",
       tag = "A", title = "Simple Model") +
  theme(legend.position = "none", 
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill='white', color = "transparent"),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.spacing.y = unit(-0.03, "cm"),
        legend.margin=margin(c(0,3,3,3)),
        axis.text = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 11),
        axis.ticks.y= element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        text = element_text(size = 11),
        strip.text = element_text(face = "italic"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.2, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background  = element_rect(colour = "white", fill = "white"),
        plot.margin = margin(2, 10, 2, 2)); fig.traj.simple

fig.traj.legend <- get_legend(fig.traj.simple)


## S18: Legend -----------------------------------------------------------------


fig.traj.legend <- ggplot(df.indiv.traj.simple) + 
  
  geom_segment(aes(x = 0, xend = as.numeric(end_inf_period_dpfp), 
                   y = indiv_sample_factor,
                   yend = indiv_sample_factor,
                   color = "Infectious Period"),
               linewidth = 0.5,
               alpha = 0.5) +
  
  
  geom_segment(data = subset(df.indiv.traj.simple, rebound_class != "None"),
               aes(x = 0, xend = as.numeric(end_inf_period_dpfp), 
                   y = indiv_sample_factor,
                   yend = indiv_sample_factor,
                   color = "Rebound Region"),
               linetype = "dotted",
               linewidth = 0.5,
               alpha = 0.8) +
  
  geom_point(aes(x = as.numeric(first_ever_neg_pred_dpfp), 
                 y = indiv_sample_factor,
                 color = "First Predicted Neg."),
             fill = "transparent",
             shape = 23,
             stroke = 1,
             size = 2.5,
             alpha = 1) +
  
  geom_point(aes(x =  as.numeric(second_cons_neg_pred_dpfp), 
                 y = indiv_sample_factor,
                 fill = "Second Predicted Neg."),
             shape = 23,
             size = 2.5,
             alpha = 1) +
  
  geom_point(aes(x =  as.numeric(second_cons_neg_true_dpfp), 
                 y = indiv_sample_factor, fill = "True Second Consecutive Neg."),
             shape = 21,
             stroke = 0,
             size = 1,
             alpha = 1) +
  
  geom_point(data = df.indiv.all,
             aes(x = days_post_first_pos,  
                 y = indiv_sample_factor,
                 fill = as.character(pos_inf)),
             alpha = 0.3, shape = 22) +
  
  
  geom_point(aes(x = as.numeric(last_pos_true_dpfp), y = indiv_sample_factor,
                 fill = "Last Culture Pos."),
             shape = 22,
             alpha = 1) +
  geom_point(aes(x = as.numeric(next_neg_true_dpfp), y = indiv_sample_factor,
                 fill = "Next Culture Neg."),
             shape = 22,
             alpha = 1) +
  
  # Set colors
  scale_fill_manual(values = c("Second Predicted Neg." = "#0F4C5E",
                               "True Second Consecutive Neg." = "red",
                               "Last Culture Pos." = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex,
                               "Next Culture Neg." = "grey62",
                               "Observed Culture Pos." = subset(traj_palette, assay == "culture" & pos_neg == "pos")$hex,
                               "Observed Culture Neg." = "grey62"),
                    breaks = c("Observed Culture Pos.", "Observed Culture Neg.",
                               "Last Culture Pos.", "Next Culture Neg.",
                               "Second Predicted Neg.", "True Second Consecutive Neg."),
                    labels = c("Observed Culture Pos.", "Observed Culture Neg.",
                               "Last Culture Pos.", "Next Culture Neg.",
                               "Isolation Ends (2nd Pred. Cons. Neg.)", "True Second Consecutive Neg.")) +  
  
  scale_color_manual(values = c("Infectious Period" = "#E4BC11",
                                "First Predicted Neg." = "#0F4C5E",
                                "Rebound Region" = zissou_pal[11]),
                     breaks = c("First Predicted Neg.",
                                "Infectious Period",
                                "Rebound Region")) +  
  
   guides(fill = guide_legend(override.aes = list(shape = c(22, 22, 22, 22, 23, 21),
                                                 size = c(1.5, 1.5, 1.5, 1.5, 2, 1),
                                                 alpha = c(0.1, 0.1, 1, 1, 1, 1)),
                             order = 1, nrow = 2),
         color = guide_legend(override.aes = list(shape = c(23, NA, NA)),
                              order = 2, nrow = 2)) +
  
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill='white', color = "transparent"),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.spacing.y = unit(-0.03, "cm"),
        legend.margin=margin(c(0,3,3,3)),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.2, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background  = element_rect(colour = "white", fill = "white"),
        plot.margin = margin(2, 10, 2, 2))

fig.traj.legend <- get_legend(fig.traj.legend)

## Combine -----------------------------------------------------------------

fig.traj <- fig.traj.simple + fig.traj.best + plot_layout(ncol = 2); fig.traj

fig.traj <- ggarrange(fig.traj, fig.traj.legend, nrow = 2, heights = c(1, 0.1))

ggsave('./outputs/figures/figS18-clinical-time-series.tiff',
       plot = fig.traj,
       device = 'tiff',
       height = 10,
       width = 7.5,
       units = 'in',
       bg = 'white')


# Fig 7A: Excess Isolation Days ---------------

df.isol <- data.frame(model = c("Ten Day", "Five Day", 
                                "Simple", "Best", "Perfect"),
                      days = c(sum(as.numeric(df.indiv.traj.ten$excess_isol_days)),
                               sum(as.numeric(df.indiv.traj.five$excess_isol_days)),
                               sum(as.numeric(df.indiv.traj.simple$second_cons_neg_excess_isol_days)),
                               sum(as.numeric(df.indiv.traj.best$second_cons_neg_excess_isol_days)),
                               sum(df.indiv.traj.perfect$excess_isol_days, na.rm = TRUE)))

df.isol$model <- factor(df.isol$model, levels = c("Ten Day", "Five Day", 
                                                  "Simple", "Best", "Perfect"))


# Make the figure
fig.isol <- ggplot(df.isol, aes(x = model, y = days, fill = model)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.7) +
  geom_beeswarm(data = subset(df.indiv.traj, second_cons_neg_excess_isol_days != 0), 
                aes(x = model, 
                    y = as.numeric(second_cons_neg_excess_isol_days) * 15,
                    fill = model),
                shape = 21, stroke = 0.3) +
  geom_beeswarm(data = subset(df.indiv.traj.threshold, excess_isol_days > 0),
                aes(x = model, y = excess_isol_days * 15, fill = model),
                shape = 21, stroke = 0.3) +
  geom_beeswarm(data = subset(df.indiv.traj.perfect, excess_isol_days > 0),
                aes(x = model, y = excess_isol_days * 15, fill = model),
                shape = 21, stroke = 0.3) +
  scale_fill_manual(values = c("Best" = "#85B068",
                               "Simple" = "#0F4C5E",
                               "Ten Day" = "deeppink4",
                               "Five Day" = "orchid4")) +
  scale_y_continuous(limits = c(0, 550),
                     expand = c(0, 0), sec.axis = sec_axis(trans=~./15, 
                                                           name="Individual days\nunnecessarily isolated")) +
  labs(y = "Cumulative days\nunnecessarily isolated") +
  coord_cartesian(clip = "off") +
  theme(legend.position = "none", 
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill='white'),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.ticks.x = element_blank(),
        text = element_text(size = 11),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.2, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        plot.margin = margin(2, 2, 2, 2)); fig.isol



# Fig 7B: Infectious Days ------------------------------------------------------

df.inf <- data.frame(model = character(),
                     rebound_class = character(),
                     sum_days = character())

for (class.ii in c("Early Rebound", "Late Rebound", "None")) {
  df.ten <- subset(df.indiv.traj.ten, rebound_class == class.ii)
  df.five <- subset(df.indiv.traj.five, rebound_class == class.ii)
  df.simple <- subset(df.indiv.traj.simple, rebound_class == class.ii)
  df.best <- subset(df.indiv.traj.best, rebound_class == class.ii)
  df.perfect <- subset(df.indiv.traj.perfect, rebound_class == class.ii)
  
  next_row <- data.frame(model = c("Ten Day", 
                                   "Five Day", 
                                   "Simple", 
                                   "Best", 
                                   "Perfect"),
                         rebound_class = rep(class.ii, 5),
                         days = c(sum(df.ten$inf_days),
                                  sum(df.five$inf_days),
                                  sum(as.numeric(df.simple$second_cons_neg_inf_days)),
                                  sum(as.numeric(df.best$second_cons_neg_inf_days)),
                                  sum(as.numeric(df.perfect$inf_days), na.rm =TRUE)))
  
  df.inf <- rbind(df.inf, next_row)
  
}

df.inf$model <- factor(df.inf$model, levels = c("Ten Day", "Five Day",
                                                "Simple", "Best", "Perfect"))
df.inf$rebound_class <- factor(df.inf$rebound_class, levels = rev(c("None", "Late Rebound",
                                                                    "Early Rebound")))


# Make the figure
fig.inf <- ggplot(df.inf, aes(x = model, y = days, fill = model)) +
  geom_bar(stat = "identity", color = "black", aes(alpha = rebound_class)) +
  geom_beeswarm(data = subset(df.indiv.traj, second_cons_neg_inf_days != 0 & rebound_class != "None"), 
                aes(x = model, y = as.numeric(second_cons_neg_inf_days) * 5, 
                    fill = "Rebound", color = "Rebound"),
                shape = 21, stroke = 0.3) +
  geom_beeswarm(data = subset(df.indiv.traj.threshold, 
                              second_cons_neg_inf_days != 0 & rebound_class != "None" &
                                model == "Five Day"), 
                aes(x = model, y = as.numeric(second_cons_neg_inf_days) * 5, 
                    fill = "Rebound", color = "Rebound"),
                shape = 21, stroke = 0.3) +
  geom_beeswarm(data = subset(df.indiv.traj.perfect, 
                              inf_days != 0 & rebound_class != "None"), 
                aes(x = model, y = as.numeric(inf_days) * 5, 
                    fill = "Rebound", color = "Rebound"),
                shape = 21, stroke = 0.3) +
  geom_beeswarm(data = subset(df.indiv.traj, second_cons_neg_inf_days != 0 & rebound_class %in% c("None")), 
                aes(x = model, y = as.numeric(second_cons_neg_inf_days) * 5, fill = model),
                shape = 21, stroke = 0.3) +
  
  geom_beeswarm(data = subset(df.indiv.traj.threshold, inf_days != 0 & rebound_class %in% c("None")), 
                aes(x = model, y = inf_days * 5, fill = model),
                shape = 21, stroke = 0.3) +
  
  scale_fill_manual(values = c("Best" = "#85B068",
                               "Simple" = "#0F4C5E",
                               "Ten Day" = "deeppink4",
                               "Five Day" = "orchid4",
                               "Rebound" = zissou_pal[11]),
                    breaks = c("Rebound")) +
  scale_color_manual(values = c("Rebound" = "black")) +
  scale_alpha_manual(values = c("None" = 0.9,
                                "Early Rebound" = 0.2,
                                "Late Rebound" = 0.6),
                     breaks = c("Early Rebound", "Late Rebound", "None"),
                     labels = c("Early rebound", "Late rebound", "Not a rebound")) +
  scale_y_continuous(limits = c(0, 68),
                     expand = c(0, 0), sec.axis = sec_axis(trans=~./5, 
                                                           name="Individual days\nstill infectious")) +
  labs(y = "Cumulative days\nstill infectious") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_none(),
         color =  guide_legend(override.aes = list(fill = zissou_pal[11]))) +
  theme(legend.position = c(0.21, 0.81), 
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, 'cm'),
        legend.background = element_rect(fill='transparent'),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.spacing.y = unit(-0.03, "cm"),
        legend.margin=margin(c(0,3,3,3)),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.ticks.x = element_blank(),
        text = element_text(size = 11),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.2, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        plot.margin = margin(2, 2, 2, 2)); fig.inf



# Fig 7C: True 2nd consecutive negative ----------------------------

df.indiv.traj.true <- subset(df.indiv.traj, second_cons_neg_true_dpfp != "Not Observed")
df.indiv.traj.true$category[as.numeric(df.indiv.traj.true$second_cons_neg_pred_dpfp) == as.numeric(df.indiv.traj.true$second_cons_neg_true_dpfp)] <- "Correct"
df.indiv.traj.true$category[as.numeric(df.indiv.traj.true$second_cons_neg_pred_dpfp) > as.numeric(df.indiv.traj.true$second_cons_neg_true_dpfp)] <- "Late"
df.indiv.traj.true$category[as.numeric(df.indiv.traj.true$second_cons_neg_pred_dpfp) < as.numeric(df.indiv.traj.true$second_cons_neg_true_dpfp)] <- "Early"

# Ten Day
df.indiv.traj.ten <- subset(df.indiv.traj.true, model == "Best")
df.indiv.traj.ten$category[10 == as.numeric(df.indiv.traj.ten$second_cons_neg_true_dpfp)] <- "Correct"
df.indiv.traj.ten$category[10 > as.numeric(df.indiv.traj.ten$second_cons_neg_true_dpfp)] <- "Late"
df.indiv.traj.ten$category[10 < as.numeric(df.indiv.traj.ten$second_cons_neg_true_dpfp)] <- "Early"
df.indiv.traj.ten$model <- "Ten Day"

# Five Day
df.indiv.traj.five <- subset(df.indiv.traj.true, model == "Best")
df.indiv.traj.five$category[5 == as.numeric(df.indiv.traj.five$second_cons_neg_true_dpfp)] <- "Correct"
df.indiv.traj.five$category[5 > as.numeric(df.indiv.traj.five$second_cons_neg_true_dpfp)] <- "Late"
df.indiv.traj.five$category[5 < as.numeric(df.indiv.traj.five$second_cons_neg_true_dpfp)] <- "Early"
df.indiv.traj.five$model <- "Five Day"

# Combine
df.indiv.traj.true <- rbind(df.indiv.traj.true,
                            df.indiv.traj.ten,
                            df.indiv.traj.five)

# Set order
df.indiv.traj.true$category <- factor(df.indiv.traj.true$category, 
                                      levels = rev(c("Correct", "Early", "Late")))

# Plot
fig.time <- ggplot() +
  geom_bar(data = df.indiv.traj.true,
           aes(x = model, 
               fill = model, alpha = category), 
           color = "black", position="fill") +
  annotate("text", x = "Best", y = 0.26, label = "Correct", 
           angle = -90, vjust = -2.6, size = 3.5, fontface = 'bold.italic') +
  annotate("text", x = "Best", y = 0.685, label = "Early", 
           angle = -90, vjust = -2.6, size = 3.5, fontface = 'bold.italic') +
  annotate("text", x = "Best", y = 0.93, label = "Late", 
           angle = -90, vjust = -2.6, size = 3.5, fontface = 'bold.italic') +
  geom_text(data = subset(df.indiv.traj.true, model == "Best"),
            aes(group = category,
                label = paste0(round(..count../sum(..count..)*100, 1), "%"),
                x = model),
            stat = "count",
            colour = "black",
            position = position_fill(vjust = 0.5),
            size = 2.5) +
  geom_text(data = subset(df.indiv.traj.true, model == "Simple"),
            aes(group = category,
                label = paste0(round(..count../sum(..count..)*100, 1), "%"),
                x = model),
            stat = "count",
            colour = "black",
            position = position_fill(vjust = 0.5),
            size = 2.5) +
  geom_text(data = subset(df.indiv.traj.true, model == "Ten Day"),
            aes(group = category,
                label = paste0(round(..count../sum(..count..)*100, 1), "%"),
                x = model),
            stat = "count",
            colour = "black",
            position = position_fill(vjust = 0.5),
            size = 2.5) +
  geom_text(data = subset(df.indiv.traj.true, model == "Five Day"),
            aes(group = category,
                label = paste0(round(..count../sum(..count..)*100, 1), "%"),
                x = model),
            stat = "count",
            colour = "black",
            position = position_fill(vjust = 0.5),
            size = 2.5) +
  scale_x_discrete(limits = c("Ten Day", "Five Day", "Simple", "Best")) +
  scale_fill_manual(values = c("Best" = "#85B068",
                               "Simple" = "#0F4C5E",
                               "Ten Day" = "deeppink4",
                               "Five Day" = "orchid4")) +
  scale_alpha_manual(values = c("Correct" = 1, 
                                "Late" = 0.3,
                                "Early" = 0.6),
                     limits = c("Correct", 
                                "Late",
                                "Early")) +
  scale_y_continuous(expand = c(0, 0), label = scales::percent) +
  coord_cartesian(clip = "off") +
  labs(y = "Identification of known\n2nd negative (% of individuals)") +
  guides(fill = guide_none()) +
  theme(legend.position = "none", 
        legend.title = element_blank(),
        legend.text = element_text(size = 8, angle = 45),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill='transparent'),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 11),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.2, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        plot.margin = margin(2, 2, 2, 2)); fig.time


table(df.indiv.traj.true$category[df.indiv.traj.true$model == "Best"])


# Fig 7D: Probabilities --------------------------------------------------------

df.indiv.traj.best.probs <- subset(df.indiv.traj.best, second_cons_neg_true_dpfp != "Not Observed")
df.indiv.traj.simple.probs <- subset(df.indiv.traj.simple, second_cons_neg_true_dpfp != "Not Observed")

# Order the dataframes so the individuals are the same
df.indiv.traj.best.probs <- df.indiv.traj.best.probs[order(df.indiv.traj.best.probs$indiv_sample), ]
df.indiv.traj.simple.probs <- df.indiv.traj.simple.probs[order(df.indiv.traj.simple.probs$indiv_sample), ]
# Check that worked (needs to be TRUE only)
unique(df.indiv.traj.best.probs$indiv_sample == df.indiv.traj.simple.probs$indiv_sample)

# Make numerics
df.indiv.traj.simple.probs$second_cons_neg_true_prob <- as.numeric(df.indiv.traj.simple.probs$second_cons_neg_true_prob)
df.indiv.traj.simple.probs$first_cons_neg_true_prob <- as.numeric(df.indiv.traj.simple.probs$first_cons_neg_true_prob)
df.indiv.traj.best.probs$second_cons_neg_true_prob <- as.numeric(df.indiv.traj.best.probs$second_cons_neg_true_prob)
df.indiv.traj.best.probs$first_cons_neg_true_prob <- as.numeric(df.indiv.traj.best.probs$first_cons_neg_true_prob)


# Calculate the differences in the probabilities
# The smaller the value, the better
# So if simple > best, this will be positive (best does better)
# And if best > simple, this will be negative (simple does better)
prob_diffs_2nd <- df.indiv.traj.simple.probs$second_cons_neg_true_prob - df.indiv.traj.best.probs$second_cons_neg_true_prob 
prob_diffs_1st <- df.indiv.traj.simple.probs$first_cons_neg_true_prob  - df.indiv.traj.best.probs$first_cons_neg_true_prob 

df.diffs <- data.frame(diffs = c(prob_diffs_2nd, prob_diffs_1st))
df.diffs$diffs[]

fig.diffs <- ggplot() +
  geom_point(aes(x = "Difference", y = -0.25), size = 32, 
             color = "#85B068", alpha = 0.4,
             shape = 22, fill = "#85B068") +
  geom_point(aes(x = "Difference", y = 0.25), size = 32, 
             color = "#0F4C5E", alpha = 0.4,
             shape = 22, fill = "#0F4C5E") +
  
  geom_hline(yintercept = 0, linewidth = 0.5) +
  annotate("text", x = "Difference", y = 0.25, label = "Best better",
           angle = -90, vjust = -3.8, size = 3) +
  annotate("text", x = "Difference", y = -0.25, label = "Simple better",
           angle = -90, vjust = -3.8, size = 3) +
  geom_quasirandom(data  = df.diffs, aes(x = "Difference", y = diffs),
                   alpha = 0.8, shape = 21, fill = "grey44", size = 2) +
  scale_y_continuous(limits = c(-0.5, 0.5),
                     breaks = seq(-0.5, 0.5, 0.5), expand = c(0, 0),
                     labels = paste0(seq(-50, 50, 50), "%"),
                     position = "left") +
  labs(y = "Difference in culture\npositive chance") +
  theme(legend.position = "none", 
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill='transparent'),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 11),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.2, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        plot.margin = margin(2, 0, 2, 2)); fig.diffs

fig.probs <- ggplot() +
  geom_point(aes(x = "Simple", y = 0.75), size = 31.5, 
             color = "#FAF2C6", 
             shape = 22, fill = "#FAF2C6") +
  geom_point(aes(x = "Best", y = 0.75), size = 31.5, 
             color = "#FAF2C6",
             shape = 22, fill = "#FAF2C6") +
  geom_point(aes(x = "Simple", y = 0.25), size = 31.5, 
             color = "grey99", 
             shape = 22, fill = "grey90") +
  geom_point(aes(x = "Best", y = 0.25), size = 31.5, 
             color = "grey90", 
             shape = 22, fill = "grey90") +
  geom_hline(yintercept = 0.5, linewidth = 0.5) +
  annotate("text", x = "Best", y = 0.75, label = "Predicted Pos.",
           angle = -90, vjust = -1.7, size = 3) +
  annotate("text", x = "Best", y = 0.25, label = "Predicted Neg.",
           angle = -90, vjust = -1.7, size = 3) +
  
  geom_quasirandom(data = df.indiv.traj.best.probs,
                   aes(x = "Best", y = second_cons_neg_true_prob, fill = "Best"), 
                   alpha = 0.8, shape = 21, size = 2) +
  geom_quasirandom(data = df.indiv.traj.best.probs,
                   aes(x = "Best", y = first_cons_neg_true_prob, fill = "Best"), 
                   alpha = 0.8, shape = 21, size = 2) +
  geom_quasirandom(data = df.indiv.traj.simple.probs,
                   aes(x = "Simple", y = second_cons_neg_true_prob, fill = "Simple"), 
                   alpha = 0.8, 
                   shape = 21, size = 2) +
  geom_quasirandom(data = df.indiv.traj.simple.probs,
                   aes(x = "Simple", y = first_cons_neg_true_prob, fill = "Simple"), 
                   alpha = 0.8, 
                   shape = 21, size = 2) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.5), expand = c(0, 0),
                     labels = paste0(seq(0, 100, 50), "%"),
                     position = "right") +
  scale_x_discrete(limits = c("Simple", "Best")) +
  labs(y = "Culture positive\nchance") +
  scale_fill_manual(values = c("Best" = "#85B068",
                               "Simple" = "#0F4C5E",
                               "Ten Day" = "deeppink4",
                               "Five Day" = "orchid4")) +
  theme(legend.position = "none", 
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill='transparent'),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(), 
        text = element_text(size = 10),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.2, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        plot.margin = margin(2, 2, 2, 0)); fig.probs


fig.7D <- fig.diffs + labs(tag = "D") + fig.probs + 
  plot_layout(nrow = 1, widths = c(2, 2)); fig.7D


# Fig S19: Sampling Times -------------------------------------------------------

# Get time between all samples
df.all.diffs <- data.frame(indiv_sample = character(),
                           dpi = numeric(),
                           time_diff = numeric())
indiv_max <- c()

for (indiv.ii in unique(df.indiv.all$indiv_sample)) {
  df.sub <- subset(df.indiv.all, indiv_sample == indiv.ii & model == "Best")
  all_diffs <- c()
  indiv_max <- c(indiv_max, max(df.sub$days_post_first_pos, na.rm = TRUE))
  
  for (row_num in 2:nrow(df.sub)) {
    df.sub <- df.sub[order(as.numeric(df.sub$days_post_first_pos)), ]
    
    all_diffs <- c(all_diffs, 
                   df.sub$days_post_first_pos[row_num] - df.sub$days_post_first_pos[row_num - 1])
    
  }
  
  df.all.diffs <- rbind(df.all.diffs,
                        data.frame(indiv_sample = indiv.ii,
                                   dpi = df.sub$days_post_first_pos[2:nrow(df.sub)],
                                   time_diff = all_diffs))
  
}



p <- ggplot(df.all.diffs) +
  geom_count(aes(x = dpi, y = time_diff),
             fill = "grey66", color = "black",
             shape = 21, alpha = 1) +
  scale_x_continuous(breaks = seq(0, 30, 4)) +
  labs(x = "days since first positive test",
       y = "days since previous test") +
  theme(legend.position = c(0.14, 0.8), 
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill='white', 
                                         color = 'black',
                                         linewidth = 0.1),
        legend.key = element_blank(),
        legend.box.spacing = unit(1, "pt"),
        legend.margin = margin(-2, 3, 2, 2),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        text = element_text(size = 11),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.2, linetype = "solid"),
        
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        plot.margin = margin(2, 2, 2, 2)); p

p <- ggExtra::ggMarginal(p, type = "histogram", fill = "grey66", stroke = 1); p

ggsave('./outputs/figures/figS19-sample-delays-vs-time.tiff',
       plot = p,
       device = 'tiff',
       height = 3,
       width = 3,
       units = 'in',
       bg = 'white')


# Prep for threshold comparisons ---------------------------------------------------

df.indiv.traj.0 <- subset(df.indiv.traj.0, indiv_sample %notin% remove_indivs)
df.indiv.traj.1 <- subset(df.indiv.traj.1, indiv_sample %notin% remove_indivs)
df.indiv.traj.2 <- subset(df.indiv.traj.2, indiv_sample %notin% remove_indivs)
df.indiv.traj.3 <- subset(df.indiv.traj.3, indiv_sample %notin% remove_indivs)
df.indiv.traj.4 <- subset(df.indiv.traj.4, indiv_sample %notin% remove_indivs)

df.indiv.traj.0$unclassified <- "No"
df.indiv.traj.0$unclassified[df.indiv.traj.0$second_cons_neg_pred_dpfp == "Not Observed"] <- "Yes"
df.indiv.traj.1$unclassified <- "No"
df.indiv.traj.1$unclassified[df.indiv.traj.1$second_cons_neg_pred_dpfp == "Not Observed"] <- "Yes"
df.indiv.traj.2$unclassified <- "No"
df.indiv.traj.2$unclassified[df.indiv.traj.2$second_cons_neg_pred_dpfp == "Not Observed"] <- "Yes"
df.indiv.traj.3$unclassified <- "No"
df.indiv.traj.3$unclassified[df.indiv.traj.3$second_cons_neg_pred_dpfp == "Not Observed"] <- "Yes"
df.indiv.traj.4$unclassified <- "No"
df.indiv.traj.4$unclassified[df.indiv.traj.4$second_cons_neg_pred_dpfp == "Not Observed"] <- "Yes"




df.indiv.traj.0$second_cons_neg_pred_dpfp[df.indiv.traj.0$unclassified == "Yes"] <- 10
df.indiv.traj.0$second_cons_neg_excess_isol_days[df.indiv.traj.0$unclassified == "Yes"] <- 10 - df.indiv.traj.0$end_inf_period_dpfp[df.indiv.traj.0$unclassified == "Yes"]
df.indiv.traj.0$second_cons_neg_inf_days[df.indiv.traj.0$unclassified == "Yes"] <-  0

df.indiv.traj.1$second_cons_neg_pred_dpfp[df.indiv.traj.1$unclassified == "Yes"] <- 10
df.indiv.traj.1$second_cons_neg_excess_isol_days[df.indiv.traj.1$unclassified == "Yes"] <- 10 - df.indiv.traj.1$end_inf_period_dpfp[df.indiv.traj.1$unclassified == "Yes"]
df.indiv.traj.1$second_cons_neg_inf_days[df.indiv.traj.1$unclassified == "Yes"] <-  0

df.indiv.traj.2$second_cons_neg_pred_dpfp[df.indiv.traj.2$unclassified == "Yes"] <- 10
df.indiv.traj.2$second_cons_neg_excess_isol_days[df.indiv.traj.2$unclassified == "Yes"] <- 10 - df.indiv.traj.2$end_inf_period_dpfp[df.indiv.traj.2$unclassified == "Yes"]
df.indiv.traj.2$second_cons_neg_inf_days[df.indiv.traj.2$unclassified == "Yes"] <-  0

df.indiv.traj.3$second_cons_neg_pred_dpfp[df.indiv.traj.3$unclassified == "Yes"] <- 10
df.indiv.traj.3$second_cons_neg_excess_isol_days[df.indiv.traj.3$unclassified == "Yes"] <- 10 - df.indiv.traj.3$end_inf_period_dpfp[df.indiv.traj.3$unclassified == "Yes"]
df.indiv.traj.3$second_cons_neg_inf_days[df.indiv.traj.3$unclassified == "Yes"] <-  0

df.indiv.traj.4$second_cons_neg_pred_dpfp[df.indiv.traj.4$unclassified == "Yes"] <- 10
df.indiv.traj.4$second_cons_neg_excess_isol_days[df.indiv.traj.4$unclassified == "Yes"] <- 10 - df.indiv.traj.4$end_inf_period_dpfp[df.indiv.traj.4$unclassified == "Yes"]
df.indiv.traj.4$second_cons_neg_inf_days[df.indiv.traj.4$unclassified == "Yes"] <-  0


# Fig 7E: Threshold Excess Isolation Days --------------------------------------

df.indiv.traj.0.simple <- subset(df.indiv.traj.0, model == "Simple")
df.indiv.traj.1.simple <- subset(df.indiv.traj.1, model == "Simple")
df.indiv.traj.2.simple <- subset(df.indiv.traj.2, model == "Simple")
df.indiv.traj.3.simple <- subset(df.indiv.traj.3, model == "Simple")
df.indiv.traj.4.simple <- subset(df.indiv.traj.4, model == "Simple")

df.indiv.traj.0 <- subset(df.indiv.traj.0, model == "Best")
df.indiv.traj.1 <- subset(df.indiv.traj.1, model == "Best")
df.indiv.traj.2 <- subset(df.indiv.traj.2, model == "Best")
df.indiv.traj.3 <- subset(df.indiv.traj.3, model == "Best")
df.indiv.traj.4 <- subset(df.indiv.traj.4, model == "Best")



df.isol.thresh.5 <- data.frame(model = c("Simple", "Best"),
                               threshold = rep(0.5, 2),
                               days = c(sum(as.numeric(df.indiv.traj.simple$second_cons_neg_excess_isol_days)),
                                        sum(as.numeric(df.indiv.traj.best$second_cons_neg_excess_isol_days))),
                               percent = c(table(df.indiv.traj.simple$unclassified)[1]/nrow(df.indiv.traj.simple),
                                           table(df.indiv.traj.best$unclassified)[1]/nrow(df.indiv.traj.best)))


df.isol.thresh.4 <- data.frame(model = c("Simple", "Best"),
                               threshold = rep(0.4, 2),
                               days = c(sum(as.numeric(df.indiv.traj.4.simple$second_cons_neg_excess_isol_days)),
                                        sum(as.numeric(df.indiv.traj.4$second_cons_neg_excess_isol_days))),
                               percent = c(table(df.indiv.traj.4.simple$unclassified)[1]/nrow(df.indiv.traj.4.simple),
                                           table(df.indiv.traj.4$unclassified)[1]/nrow(df.indiv.traj.4)))


df.isol.thresh.3 <- data.frame(model = c("Simple", "Best"),
                               threshold = rep(0.3, 2),
                               days = c(sum(as.numeric(df.indiv.traj.3.simple$second_cons_neg_excess_isol_days)),
                                        sum(as.numeric(df.indiv.traj.3$second_cons_neg_excess_isol_days))),
                               percent = c(table(df.indiv.traj.3.simple$unclassified)[1]/nrow(df.indiv.traj.3.simple),
                                           table(df.indiv.traj.3$unclassified)[1]/nrow(df.indiv.traj.3)))


df.isol.thresh.2 <- data.frame(model = c("Simple", "Best"),
                               threshold = rep(0.2, 2),
                               days = c(sum(as.numeric(df.indiv.traj.2.simple$second_cons_neg_excess_isol_days)),
                                        sum(as.numeric(df.indiv.traj.2$second_cons_neg_excess_isol_days))),
                               percent = c(table(df.indiv.traj.2.simple$unclassified)[1]/nrow(df.indiv.traj.2.simple),
                                           table(df.indiv.traj.2$unclassified)[1]/nrow(df.indiv.traj.2)))

df.isol.thresh.1 <- data.frame(model = c("Simple", "Best"),
                               threshold = rep(0.1, 2),
                               days = c(sum(as.numeric(df.indiv.traj.1.simple$second_cons_neg_excess_isol_days)),
                                        sum(as.numeric(df.indiv.traj.1$second_cons_neg_excess_isol_days))),
                               percent = c(table(df.indiv.traj.1.simple$unclassified)[1]/nrow(df.indiv.traj.1.simple),
                                           table(df.indiv.traj.1$unclassified)[1]/nrow(df.indiv.traj.1)))

df.isol.thresh.0 <- data.frame(model = c("Simple", "Best"),
                               threshold = rep(0, 2),
                               days = c(sum(as.numeric(df.indiv.traj.0.simple$second_cons_neg_excess_isol_days)),
                                        sum(as.numeric(df.indiv.traj.0$second_cons_neg_excess_isol_days))),
                               percent = c(table(df.indiv.traj.0.simple$unclassified)[1]/nrow(df.indiv.traj.0.simple),
                                           table(df.indiv.traj.0$unclassified)[1]/nrow(df.indiv.traj.0)))


df.isol.thresh <- rbind(df.isol.thresh.5, df.isol.thresh.4, df.isol.thresh.3,
                        df.isol.thresh.2, df.isol.thresh.1)

df.isol.thresh$model <- factor(df.isol.thresh$model, levels = c("Simple",
                                                                "Best", "Perfect"))

df.indiv.traj.0$second_cons_neg_excess_isol_days_adj <- df.indiv.traj.0$second_cons_neg_excess_isol_days
df.indiv.traj.0$second_cons_neg_excess_isol_days_adj[df.indiv.traj.0$second_cons_neg_excess_isol_days == 1] <- rep(c(0.5, 1, 1.5), 25)

df.indiv.traj.1$second_cons_neg_excess_isol_days_adj <- df.indiv.traj.1$second_cons_neg_excess_isol_days
df.indiv.traj.1$second_cons_neg_excess_isol_days_adj[df.indiv.traj.1$second_cons_neg_excess_isol_days == 1] <- rep(c(0.5, 1, 1.5), 25)

df.indiv.traj.2$second_cons_neg_excess_isol_days_adj <- df.indiv.traj.2$second_cons_neg_excess_isol_days
df.indiv.traj.2$second_cons_neg_excess_isol_days_adj[df.indiv.traj.2$second_cons_neg_excess_isol_days == 1] <- rep(c(0.5, 1, 1.5), 25)

df.indiv.traj.3$second_cons_neg_excess_isol_days_adj <- df.indiv.traj.3$second_cons_neg_excess_isol_days
df.indiv.traj.3$second_cons_neg_excess_isol_days_adj[df.indiv.traj.3$second_cons_neg_excess_isol_days == 1] <- rep(c(0.5, 1, 1.5), 25)

df.indiv.traj.4$second_cons_neg_excess_isol_days_adj <- df.indiv.traj.4$second_cons_neg_excess_isol_days
df.indiv.traj.4$second_cons_neg_excess_isol_days_adj[df.indiv.traj.4$second_cons_neg_excess_isol_days == 1] <- rep(c(0.5, 1, 1.5), 25)

df.isol.thresh$threshold <- paste0(df.isol.thresh$threshold * 100, "%")
df.isol.thresh$threshold <- factor(df.isol.thresh$threshold, 
                                   levels = c("10%", "20%", "30%", "40%", "50%"))

fig.7E <- ggplot(df.isol.thresh, aes(x = threshold, y = days, color = model)) +

  geom_bar(data = subset(df.isol.thresh, model == "Best"), 
           aes(fill = model),
           stat = "identity",
           color = "black") +
  geom_hline(data = subset(df.isol, model == "Ten Day"),
             aes(yintercept = days, color = model),
             linewidth = 1) +
  geom_hline(data = subset(df.isol, model == "Five Day"),
             aes(yintercept = days, color = model),
             linewidth = 1) +
  geom_beeswarm(data = subset(df.indiv.traj.1, as.numeric(second_cons_neg_excess_isol_days_adj) != 0), 
                aes(x = "10%", 
                    y = as.numeric(second_cons_neg_excess_isol_days_adj) * 15,
                    fill = model),
                color = "black",
                shape = 21, stroke = 0.3) +
  geom_beeswarm(data = subset(df.indiv.traj.2, as.numeric(second_cons_neg_excess_isol_days_adj) != 0), 
                aes(x = "20%", 
                    y = as.numeric(second_cons_neg_excess_isol_days_adj) * 15,
                    fill = model),
                color = "black",
                shape = 21, stroke = 0.3) +
  geom_beeswarm(data = subset(df.indiv.traj.3, as.numeric(second_cons_neg_excess_isol_days_adj) != 0), 
                aes(x = "30%", 
                    y = as.numeric(second_cons_neg_excess_isol_days_adj) * 15,
                    fill = model),
                color = "black",
                shape = 21, stroke = 0.3) +
  geom_beeswarm(data = subset(df.indiv.traj.4, as.numeric(second_cons_neg_excess_isol_days_adj) != 0), 
                aes(x = "40%", 
                    y = as.numeric(second_cons_neg_excess_isol_days_adj) * 15,
                    fill = model),
                color = "black",
                shape = 21, stroke = 0.3) +
  geom_beeswarm(data = subset(df.indiv.traj.best, as.numeric(second_cons_neg_excess_isol_days) != 0), 
                aes(x = "50%", 
                    y = as.numeric(second_cons_neg_excess_isol_days) * 15,
                    fill = model),
                color = "black",
                shape = 21, stroke = 0.3) +
  geom_line(data = subset(df.isol.thresh, model == "Simple"),
            aes(x = threshold, y = days, group = model),
            alpha = 0.8,
            linewidth = 1, color = "#0F4C5E") +
  geom_point(data = subset(df.isol.thresh, model == "Simple"),
             aes(x = threshold, y = days, fill = model), 
             alpha = 0.8,
             shape = 21, size = 2.5, color = "black") +
  
  scale_color_manual(values = c("Best" = "#85B068",
                               "Simple" = "#0F4C5E",
                               "Ten Day" = "deeppink4",
                               "Five Day" = "orchid4")) +
  scale_fill_manual(values = c("Best" = "#85B068",
                               "Simple" = "#0F4C5E",
                               "Ten Day" = "deeppink4",
                               "Five Day" = "orchid4")) +
  scale_y_continuous(limits = c(0, 560),
                     expand = c(0, 0), sec.axis = sec_axis(trans=~./15, 
                                                           name="Individual days\nunnecessarily isolated")) +
  labs(y = "Cumulative days\nunnecessarily isolated", 
       x = "Threshold probability of\nbeing culture positive") +
  coord_cartesian(clip = "off") +
  guides(color = "none",
         fill = "none") +
  theme(legend.position = "none", 
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, 'cm'),
        legend.background = element_rect(fill='white'),
        legend.key=element_blank(),
        legend.margin = margin(t=-10),
        legend.box.spacing = unit(1, "pt"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        axis.title.x = element_text(margin = margin(t = -15)),
        axis.text.x = element_text(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 11),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.2, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        plot.margin = margin(2, 2, 2, 2)); fig.7E


# Fig 7F: Threshold Infectious Days --------------------------------------------

df.inf.thresh.5 <- data.frame(model = c("Simple", "Best"),
                               threshold = rep(0.5, 2),
                               days = c(sum(as.numeric(df.indiv.traj.simple$second_cons_neg_inf_days)),
                                        sum(as.numeric(df.indiv.traj.best$second_cons_neg_inf_days))),
                              percent = c(table(df.indiv.traj.simple$unclassified)[1]/nrow(df.indiv.traj.simple),
                                          table(df.indiv.traj.best$unclassified)[1]/nrow(df.indiv.traj.best)))


df.inf.thresh.4 <- data.frame(model = c("Simple", "Best"),
                               threshold = rep(0.4, 2),
                               days = c(sum(as.numeric(df.indiv.traj.4.simple$second_cons_neg_inf_days)),
                                        sum(as.numeric(df.indiv.traj.4$second_cons_neg_inf_days))),
                               percent = c(table(df.indiv.traj.4.simple$unclassified)[1]/nrow(df.indiv.traj.4.simple),
                                           table(df.indiv.traj.4$unclassified)[1]/nrow(df.indiv.traj.4)))


df.inf.thresh.3 <- data.frame(model = c("Simple", "Best"),
                               threshold = rep(0.3, 2),
                               days = c(sum(as.numeric(df.indiv.traj.3.simple$second_cons_neg_inf_days)),
                                        sum(as.numeric(df.indiv.traj.3$second_cons_neg_inf_days))),
                               percent = c(table(df.indiv.traj.3.simple$unclassified)[1]/nrow(df.indiv.traj.3.simple),
                                          table(df.indiv.traj.3$unclassified)[1]/nrow(df.indiv.traj.3)))


df.inf.thresh.2 <- data.frame(model = c("Simple", "Best"),
                               threshold = rep(0.2, 2),
                               days = c(sum(as.numeric(df.indiv.traj.2.simple$second_cons_neg_inf_days)),
                                        sum(as.numeric(df.indiv.traj.2$second_cons_neg_inf))),
                              percent = c(table(df.indiv.traj.2.simple$unclassified)[1]/nrow(df.indiv.traj.2.simple),
                                          table(df.indiv.traj.2$unclassified)[1]/nrow(df.indiv.traj.2)))

df.inf.thresh.1 <- data.frame(model = c("Simple",  "Best"),
                               threshold = rep(0.1, 2),
                               days = c(sum(as.numeric(df.indiv.traj.1.simple$second_cons_neg_inf_days)),
                                        sum(as.numeric(df.indiv.traj.1$second_cons_neg_inf_days))),
                               percent = c(table(df.indiv.traj.1.simple$unclassified)[1]/nrow(df.indiv.traj.1.simple),
                                           table(df.indiv.traj.1$unclassified)[1]/nrow(df.indiv.traj.1)))

df.inf.thresh.0 <- data.frame(model = c("Simple", "Best"),
                              threshold = rep(0, 2),
                              days = c(sum(as.numeric(df.indiv.traj.0.simple$second_cons_neg_inf_days)),
                                       sum(as.numeric(df.indiv.traj.0$second_cons_neg_inf_days))),
                              percent = c(table(df.indiv.traj.0.simple$unclassified)[1]/nrow(df.indiv.traj.0.simple),
                                          table(df.indiv.traj.0$unclassified)[1]/nrow(df.indiv.traj.0)))


df.inf.thresh <- rbind(df.inf.thresh.5, df.inf.thresh.4, df.inf.thresh.3,
                        df.inf.thresh.2, df.inf.thresh.1)


df.inf.thresh$model <- factor(df.inf.thresh$model, levels = c("Simple", "Best", "Perfect"))

df.inf.thresh$threshold <- paste0(df.inf.thresh$threshold * 100, "%")



fig.7F <- ggplot(df.inf.thresh, aes(x = threshold, 
                                    y = days, color = model)) +
  
  #geom_point() +
  #geom_line() +
  geom_bar(data = subset(df.inf.thresh, model == "Best"), 
           aes(fill = model),
           stat = "identity",
           color = "black") +
  geom_hline(aes(yintercept = sum(df.inf$days[df.inf$model == "Ten Day"]), 
                 color = "Ten Day"),
             linewidth = 1) +
  geom_hline(aes(yintercept = sum(df.inf$days[df.inf$model == "Five Day"]), 
                 color = "Five Day"),
             linewidth = 1) +
  
  geom_beeswarm(data = subset(df.indiv.traj.1, as.numeric(second_cons_neg_inf_days) != 0), 
                aes(x = "10%", 
                    y = as.numeric(second_cons_neg_inf_days) * 5,
                    fill = rebound_class), color = "black",
                shape = 21, stroke = 0.3) +
  geom_beeswarm(data = subset(df.indiv.traj.2, as.numeric(second_cons_neg_inf_days) != 0), 
                aes(x = "20%", 
                    y = as.numeric(second_cons_neg_inf_days) * 5,
                    fill = rebound_class),
                color = "black",
                shape = 21, stroke = 0.3) +
  geom_beeswarm(data = subset(df.indiv.traj.3, as.numeric(second_cons_neg_inf_days) != 0), 
                aes(x = "30%", 
                    y = as.numeric(second_cons_neg_inf_days) * 5,
                    fill = rebound_class), color = "black",
                shape = 21, stroke = 0.3) +
  geom_beeswarm(data = subset(df.indiv.traj.4, as.numeric(second_cons_neg_inf_days) != 0), 
                aes(x = "40%", 
                    y = as.numeric(second_cons_neg_inf_days) * 5,
                    fill = rebound_class), color = "black",
                shape = 21, stroke = 0.3) +
  geom_beeswarm(data = subset(df.indiv.traj.best, as.numeric(second_cons_neg_inf_days) != 0), 
                aes(x = "50%", 
                    y = as.numeric(second_cons_neg_inf_days) * 5,
                    fill = rebound_class), color = "black",
                shape = 21, stroke = 0.3) +
  geom_line(data = subset(df.inf.thresh, model == "Simple"),
            aes(x = threshold, y = days, group = model),
            linewidth = 1, color = "#0F4C5E", alpha = 0.8) +
  geom_point(data = subset(df.inf.thresh, model == "Simple"),
             aes(x = threshold, y = days, fill = model), 
             color = "black", alpha = 0.8,
             shape = 21, size = 2.5) +
  scale_color_manual(values = c("Best" = "#85B068",
                              "Simple" = "#0F4C5E",
                              "Ten Day" = "deeppink4",
                              "Five Day" = "orchid4")) +
  scale_fill_manual(values = c("Best" = "#85B068",
                               "Simple" = "#0F4C5E",
                               "Ten Day" = "deeppink4",
                               "Five Day" = "orchid4",
                               "Early Rebound" = zissou_pal[11],
                               "Late Rebound" = zissou_pal[11],
                               "None" = "#85B068")) +
  scale_y_continuous(limits = c(0, 70),
                     expand = c(0, 0), sec.axis = sec_axis(trans=~./5, 
                                                           name="Individual days\nstill infectious")) +
  labs(y = "Cumulative days\nstill infectious", 
       x = "Threshold probability of\nbeing culture positive") +
  coord_cartesian(clip = "off") +
  theme(legend.position = "none", 
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill='white'),
        legend.key=element_blank(),
        legend.box.spacing = unit(1, "pt"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        axis.title.x = element_text(margin = margin(t = -15)),
        axis.text.x = element_text(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 11),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.2, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.05, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        plot.margin = margin(2, 2, 2, 2)); fig.7F



# Combine and save ------


fig7 <- fig.isol + labs(tag = "A") + 
  fig.inf + labs(tag = "B") +
  fig.time + labs(tag = "C") +
  fig.7D + labs(tag = "D") + 
  fig.7E + labs(tag = "E") + 
  fig.7F +  labs(tage = "F") + 
  plot_layout(ncol = 3, byrow = FALSE,
              widths = c(1, 0.9, 1)); fig7



ggsave('./outputs/figures/fig7-clinical-metrics.tiff',
       plot = fig7,
       device = 'tiff',
       height = 5,
       width = 10,
       units = 'in',
       bg = 'white')


