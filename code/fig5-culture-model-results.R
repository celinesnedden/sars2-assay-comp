# Generates figure 5 from the manuscript, including raw data & model fits,
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


# Add jitter to < LOD values for visualization ---------------------------------

dat.culture$pos_inf_jitter <- NA
for (row_num in 1:nrow(dat.culture)){
  if (dat.culture$pos_inf[row_num] == 1){
    # culture POS
    dat.culture$pos_inf_jitter[row_num] <- runif(1, 0.52, 0.98)
  }
  else if (dat.culture$pos_inf[row_num] == 0) {
    # culture NEG
    dat.culture$pos_inf_jitter[row_num] <- runif(1, 0.02, 0.48)
  }
}


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


# 4A: Raw data -----------------------------------------------------------------

dat.culture$val_total[dat.culture$cens_total == 1] <- runif(length(dat.culture$val_total[dat.culture$cens_total == 1]),
                                                            -1.9, -0.1)

fig.data <- ggplot() + 
  geom_rect(aes(xmin = -2, xmax = 0, ymin = 0, ymax = 1), fill = "grey56", color = "black") +
  geom_hline(yintercept = 0.5, linewidth = 0.6, color = "black") +
  geom_vline(xintercept = 0, linewidth = 1, color = "black") +
  geom_point(data = dat.culture,
             aes(x = as.numeric(val_total), 
                 y = as.numeric(pos_inf_jitter),
                 fill = as.character(pos_inf)),
             alpha = 0.6, size = 1.5, shape = 21) +
  scale_fill_manual(values = c("0" = brewer.pal(11, "RdBu")[10], "1" = brewer.pal(11, "RdBu")[2])) +
  labs(y = "",
       x = "log10 total RNA copies / sample",
       color = "log10 \n dose \n (pfu)",
       tag = "A") +
  guides(fill = "none",
         color = guide_colorbar(frame.colour = "black",
                                #title.hjust = -1,
                                title.vjust = 0)) +
  theme(legend.position = "none",
        axis.title = element_text(size = 8.5),
        axis.text = element_text(size = 7),
        axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1, size = 8.5),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(colour = c("transparent", rep("black", 7))),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        plot.margin = margin(2, 2, 2, 2)) + 
  scale_y_continuous(breaks = c(0.25, 0.75), limits = c(0, 1),
                     labels = c("Culture\nNegative", "Culture\nPositive"), 
                     expand = c(0, 0)) +
  scale_x_continuous(breaks = c(-1, 0, 2, 4, 6, 8, 10, 12), limits = c(-2, 12.5),
                     labels = c("< LOD", 0, 2, 4, 6, 8, 10, 12),
                     expand = c(0, 0)) +
  geom_hline(yintercept = 1, size = 0.6, color = "black") +
  geom_hline(yintercept = 0, size = 0.6, color = "black"); fig.data


# 4B: Model Probabilities ------------------------------------------------------

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
pos_pred <- prob_mean[dat.culture$pos_inf == 1]

# Probabilities of true negatives
neg_pred <- prob_mean[dat.culture$pos_inf == 0]



# Place in a data frame
prediction.df <- as.data.frame(matrix(data = NA, 
                                      nrow = length(c(pos_pred, neg_pred)),
                                      ncol = 3))
prediction.df[, 1] <- "Simple"
prediction.df[, 2] <- c(pos_pred, neg_pred)
prediction.df[, 3] <- c(rep("Positive", length(pos_pred)),
                        rep("Negative", length(neg_pred)))



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
pos_pred <- prob_mean[dat.culture$pos_inf == 1]

# Probabilities of true negatives
neg_pred <- prob_mean[dat.culture$pos_inf == 0]

# Percent of positives that are correct:
sum(pos_pred >= 0.5)/(length(pos_pred) + nrow(subset(dat, pos_inf == 1 & cens_total == 1))) * 100

# Percent of negatives that are correct:
(sum(neg_pred < 0.5)+ nrow(subset(dat, pos_inf == 0 & cens_total == 1)))/(length(neg_pred) + nrow(subset(dat, pos_inf == 0 & cens_total == 1))) * 100


# Place in joint data frame
prediction.df <- rbind(prediction.df, 
                       cbind(rep("Best", length(c(pos_pred, neg_pred))),
                             c(pos_pred, neg_pred),
                             c(rep("Positive", length(pos_pred)), 
                               rep("Negative", length(neg_pred)))))
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
  facet_grid(.~ "Observed Culture Negative") +
  scale_fill_manual(values = c("Best" = "#85B068",
                               "Simple" = "#0F4C5E")) +
  scale_color_manual(values = c("Best" = "#56783F",
                                "Simple" = "#0F4C5E")) +
  guides(color = "none",
         fill = guide_legend(override.aes = list(alpha = 0.6,
                                                 color = c( "#0F4C5E", "#56783F")))) +
  labs(x = "Chance of Culture Neg. (%)", y = "Density",
       fill = "Model", tag = "B") + 
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0),
                     breaks = seq(0, 1, 0.25), 
                     labels = c("0","25", "50", "75", "100")) +
  scale_y_continuous(limits = c(0, 2.6), breaks = seq(0, 2, 1),
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
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(2, 6, 2, 2)); fig3B.1


# Plot all observed culture positives
fig3B.2 <- ggplot(subset(prediction.df, posneg == "Positive")) +
  geom_vline(xintercept = threshold, size = 0.7, color = zissou_pal[1], alpha = 0.8,
             linetype = "dashed") +
  geom_density(aes(x = (as.numeric(prob_mean)),
                   fill = factor(model, levels = c("Simple", "Best")),
                   color = factor(model, levels = c("Simple", "Best"))),
               alpha = 0.4) +
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
  scale_y_continuous(limits = c(0, 2.6), expand = c(0, 0),
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
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(2, 2, 2, 0)); fig3B.2


fig.perf <- fig3B.1 + fig3B.2 + plot_layout(ncol = 2); fig.perf


# 4AB: Combine into Top Panel --------------------------------------------------

fig.top <- fig.data + fig.perf + plot_layout(nrow = 1, widths = c(1, 0.7)); fig.top


# 4C: Heat Map -----------------------------------------------------------------

x_vals <- seq(from = 0, to = 12, by = 1)



# 4C: Simple predictions -------------------------------------------------------

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
y_vals <- exp(gamma_median + psi_median * x_vals)/(1 + exp(gamma_median + psi_median * x_vals))
lines.simple.df <- as.data.frame(cbind(x_vals, y_vals))
lines.simple.df$predictor <- ""


## Determine in what range simple model has different outcome ------------------

range.df <- data.frame(dose = NA, dpi = NA, sp = NA, age = NA, 
                       cell = NA, assay = NA, tg = NA,
                       x_val = NA, median_prob = NA, quantile_probs = NA,
                       simple_prob = NA)

for (dose.option in c(4, 5.5, 7)) {
    for (assay.option in c(0, 1)){
      for (tg.option in c(1, 2, 3)) {
        
        if (tg.option == 1){tg.sample <-  param.best.df.wide$'psiTG[1]'}
        else if(tg.option == 2){tg.sample <- param.best.df.wide$'psiTG[2]'}
        else if(tg.option == 3){tg.sample <- param.best.df.wide$'psiTG[3]'}
        
        for (cell.option in c(1, 2, 3)) {
          if (cell.option == 1){cell.sample <-  param.best.df.wide$'psiCELL[1]'}
          else if(cell.option == 2){cell.sample <- param.best.df.wide$'psiCELL[2]'}
          else if(cell.option == 3){cell.sample <- param.best.df.wide$'psiCELL[3]'}
          
          for (age.option in c(1, 2, 3)){
            if (age.option == 1){age.sample <-  param.best.df.wide$'psiAGE[1]'}
            else if(age.option == 2){age.sample <- param.best.df.wide$'psiAGE[2]'}
            else if(age.option == 3){age.sample <- param.best.df.wide$'psiAGE[3]'}
            
            for (sp.option in c(1, 2, 3)){
              if (sp.option == 1){sp.sample <-  param.best.df.wide$'psiSP[1]'}
              else if(sp.option == 2){sp.sample <- param.best.df.wide$'psiSP[2]'}
              else if(sp.option == 3){sp.sample <- param.best.df.wide$'psiSP[3]'}
              
              for (dpi.option in c(1, 2, 3)) {
                if (dpi.option == 1){dpi.sample <-  param.best.df.wide$'psiDPI[1]'}
                else if(dpi.option == 2){dpi.sample <- param.best.df.wide$'psiDPI[2]'}
                else if(dpi.option == 3){dpi.sample <- param.best.df.wide$'psiDPI[3]'}
                
                for (x_val in x_vals) {
                  y <- param.best.df.wide$gamma +
                    param.best.df.wide$psiT * x_val +
                    dpi.sample +
                    age.sample +
                    sp.sample +
                    param.best.df.wide$psiASSAY * assay.option +
                    tg.sample +
                    param.best.df.wide$psiDOSE * dose.option +
                    cell.sample
                  y_prob <- round(exp(y)/(1+exp(y)), 2) * 100
                  y_med <- median(y_prob) 
                  y_quant <- paste0(quantile(y_prob, c(0.05, 0.95)), collapse = ", ")
                  
                  range.df <- rbind(range.df,
                                    data.frame(
                                      dose = dose.option,
                                      dpi = dpi.option,
                                      sp = sp.option,
                                      age = age.option,
                                      cell = cell.option, 
                                      assay = assay.option, 
                                      tg = tg.option,
                                      x_val = x_val, 
                                      median_prob = y_med, 
                                      quantile_probs = y_quant,
                                      simple_prob = lines.simple.df$y_vals[lines.simple.df$x_vals == x_val]))
                  
                }
              }
            }
          }
        }
      }
    }
}

View(subset(range.df, (median_prob > 0.5 & simple_prob <= 0.5) |
              median_prob <= 0.5 & simple_prob > 0.5))

sort(unique(range.df$x_val[(range.df$median_prob > 0.5 & range.df$simple_prob <= 0.5) |
                             range.df$median_prob  <= 0.5 & range.df$simple_prob > 0.5]))

range.df <- range.df[-1, ]


## Make the Supplemental Table for 90% CrI -------------------------------------

tblS7 <- subset(range.df, select = -c(median_prob, simple_prob))

tblS7 <- tblS7 %>%
  pivot_wider(names_from = x_val, values_from = c(quantile_probs))


tblS7.doses <- subset(tblS7, dose %in% c(4, 5.5, 7) & sp == 1 & tg == 1 &
                        assay == 1 & age == 2 & cell == 1 & dpi == 2) %>%
  rename(predictor = dose) %>%
  select(-c(sp, tg, assay, cell, dpi, age)); tblS7.doses

tblS7.species <- subset(tblS7, dose == 5.5 & tg == 1 &
                          assay == 1 & age == 2 & cell == 1 & dpi == 2) %>%
  rename(predictor = sp) %>%
  select(-c(dose, tg, assay, cell, dpi, age)); tblS7.species

tblS7.tg <- subset(tblS7, dose == 5.5 & sp == 1 & 
                     assay == 1 & age == 2 & cell == 1 & dpi == 2) %>%
  rename(predictor = tg) %>%
  select(-c(dose, sp, assay, cell, age, dpi)); tblS7.tg

tblS7.dpi <- subset(tblS7, dose == 5.5 & sp == 1 & tg == 1 & 
                     assay == 1 & age == 2 & cell == 1) %>%
  rename(predictor = dpi) %>%
  select(-c(dose, sp, assay, cell, age, tg)); tblS7.dpi

tblS7.cell <- subset(tblS7, dose == 5.5 & sp == 1 & dpi == 2 & tg == 1 & 
                      assay == 1 & age == 2) %>%
  rename(predictor = cell) %>%
  select(-c(dose, sp, assay, age, tg, dpi)); tblS7.cell

tblS7.assay<- subset(tblS7, dose == 5.5 & sp == 1 & dpi == 2 & tg == 1 & 
                       age == 2 & cell == 1) %>%
  rename(predictor = assay) %>%
  select(-c(dose, sp, cell, age, tg, dpi)); tblS7.assay

tblS7.age<- subset(tblS7, dose == 5.5 & sp == 1 & dpi == 2 & tg == 1 & 
                       assay == 1 & cell == 1) %>%
  rename(predictor = age) %>%
  select(-c(dose, sp, cell, assay, tg, dpi)); tblS7.age



tblS7.combined <- rbind(tblS7.doses, tblS7.dpi, tblS7.species, tblS7.age,
                        tblS7.cell, tblS7.assay, tblS7.tg); tblS7.combined
tblS7.combined$group <- c(rep("Dose", 3),
                          rep("DPI", 3),
                          rep("Species", 3),
                          rep("Age", 3),
                          rep("Cell Line", 3),
                          rep("Assay", 2),
                          rep("Target Gene", 3))
tblS7.combined <- tblS7.combined %>%
  select(c(group, predictor, everything())); tblS7.combined


write.csv(tblS7.combined, file = "./outputs/tables/tblS8-culture-90CrI.csv",
          row.names = FALSE)



## Get quantiles for reported stats --------------------------------------------



## Dose
subset(range.df, dose %in% c(4, 7) & sp == 1 & tg == 1 & dpi == 2 & 
         age == 2 & cell == 1 & assay == 1 & x_val == 7)

## DPI
subset(range.df, dose %in% c(5.5) & sp == 1 & tg == 1 & dpi %in% c(1, 2, 3) & 
         age == 2 & cell == 1 & assay == 1 & x_val == 7)

## AGE
subset(range.df, dose %in% c(5.5) & sp == 1 & tg == 1 & dpi %in% c(2) & 
         age %in% c(1, 2, 3) & cell == 1 & assay == 1 & x_val == 7)

## SPECIES
subset(range.df, dose %in% c(5.5) & sp %in% c(1, 2, 3) & tg == 1 & dpi %in% c(2) & 
         age == 2 & cell == 1 & assay == 1 & x_val == 7)

## CELL
subset(range.df, dose %in% c(5.5) & sp == 1 & tg == 1 & dpi %in% c(2) & 
         age == 2 & cell %in% c(1, 2, 3) & assay == 1 & x_val == 7)

## ASSAY
subset(range.df, dose %in% c(5.5) & sp == 1 & tg == 1 & dpi %in% c(2) & 
         age == 2 & cell == 1 & assay %in% c(0 , 1) & x_val == 7)

## TG
subset(range.df, dose %in% c(5.5) & sp == 1 & tg %in% c(1, 2, 3) & dpi %in% c(2) & 
         age == 2 & cell == 1 & assay == 1 & x_val == 7)

# 4C: SIMPLE -----------------------------------------------------------------


fig.simple <- ggplot(as.data.frame(lines.simple.df), 
                     aes(x = as.numeric(x_vals), 
                         y = predictor, 
                         fill = as.numeric(y_vals))) +
  geom_tile() +
  geom_rect(aes(xmin = -0.5, xmax = 12.5, ymin = 0.5, ymax = 1.5), 
            fill = "transparent", color = "grey33", size = 1.5) +
  geom_text(aes(label = (round(as.numeric(y_vals), 2)) * 100), color = "black", size = 2.5) +
  scale_fill_distiller(palette = "RdBu") +
  coord_cartesian(xlim = c(-0.5, 12.5), clip = "off") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(-0.5, 12.5, 2), 
                     labels = seq(0, 12, 2),
                     expand = c(0, 0)) +
  labs(x = "log10 total RNA copies / sample",
       y = "Simple") +
  theme(legend.position = "none",
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 8),
        axis.text.y = element_text(angle = 0, 
                                   vjust = 1, hjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(2, 2, 2, 2)); fig.simple


# 4C: Heat-map Demographics ---------------------------------------------------


### Age -----------------------------------------------------------------------

age.df <- data.frame(tRNA = NA,
                     age = NA,
                     probs = NA)


for (x_val in x_vals) {
  y <- param.best.df.wide$gamma +
    param.best.df.wide$psiT * x_val +
    param.best.df.wide$'psiDPI[2]' +
    param.best.df.wide$'psiAGE[1]' +
    param.best.df.wide$'psiSP[1]' +
    param.best.df.wide$psiASSAY * 1 +
    param.best.df.wide$'psiTG[1]'  +
    param.best.df.wide$psiDOSE *  5.5 +
    param.best.df.wide$'psiCELL[1]'
  
  y_juv <- median(exp(y)/(1+exp(y)))
  
  y <- param.best.df.wide$gamma +
    param.best.df.wide$psiT * x_val +
    param.best.df.wide$'psiDPI[2]' +
    param.best.df.wide$'psiAGE[2]' +
    param.best.df.wide$'psiSP[1]' +
    param.best.df.wide$psiASSAY * 1 +
    param.best.df.wide$'psiTG[1]'  +
    param.best.df.wide$psiDOSE *  5.5 +
    param.best.df.wide$'psiCELL[1]'
  
  y_ad <- median(exp(y)/(1+exp(y)))
  
  y <- param.best.df.wide$gamma +
    param.best.df.wide$psiT * x_val +
    param.best.df.wide$'psiDPI[2]' +
    param.best.df.wide$'psiAGE[3]' +
    param.best.df.wide$'psiSP[1]' +
    param.best.df.wide$psiASSAY * 1 +
    param.best.df.wide$'psiTG[1]'  +
    param.best.df.wide$psiDOSE *  5.5 +
    param.best.df.wide$'psiCELL[1]'
  
  y_ger <- median(exp(y)/(1+exp(y)))
  
  age.df <- rbind(age.df, 
                  data.frame(tRNA = x_val, 
                             age = "Juvenile", 
                             probs = y_juv),
                  data.frame(tRNA = x_val, 
                             age = "Adult",
                             probs = y_ad),
                  data.frame(tRNA = x_val, 
                             age = "Geriatric",
                             probs = y_ger))
}

age.df <- age.df[-1, ]

fig.age <- ggplot(age.df, aes(x = as.numeric(tRNA), 
                              y = factor(age, levels = rev(c("Juvenile", "Adult", "Geriatric"))),
                              fill = as.numeric(probs))) +
  geom_segment(x = -2, xend = -2, y = -4, yend = 4, linewidth = 0.2) +
  geom_tile() + 
  geom_rect(aes(xmin = 2.5, xmax = 7.5, ymin = 0.5, ymax = 3.5), 
            fill = "transparent", color = "grey33", size = 1.5) +
  geom_text(aes(label = (round(as.numeric(probs), 2)) * 100), color = "black", size = 2.5) +
  scale_fill_gradientn(colours = rev(brewer_pal), limits = c(0, 1)) +
  labs(x = "log10 total RNA copies / sample", fill = "P(Culture Positive)",
       y = "Age Class") +
  coord_cartesian(xlim = c(-0.5, 12.5), ylim = c(0.5, 3.5), clip = "off") +
  scale_x_continuous(breaks = seq(-0.5, 12.5, 2), 
                     labels = seq(0, 12, 2),
                     expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(legend.position = "none",
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 8),
        axis.text.y = element_text(angle = 45, size = 6.5,
                                   face = c("plain", "bold", "plain"), 
                                   vjust = 1, hjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(2, 2, 2, 2)); fig.age


### Species ---------------------------------------------------------------------

sp.df <- data.frame(tRNA = NA,
                     sp = NA,
                     probs = NA)


for (x_val in x_vals) {
  y <- param.best.df.wide$gamma +
    param.best.df.wide$psiT * x_val +
    param.best.df.wide$'psiDPI[2]' +
    param.best.df.wide$'psiAGE[2]' +
    param.best.df.wide$'psiSP[1]' +
    param.best.df.wide$psiASSAY * 1 +
    param.best.df.wide$'psiTG[1]'  +
    param.best.df.wide$psiDOSE *  5.5 +
    param.best.df.wide$'psiCELL[1]'
  
  y_rm <- median(exp(y)/(1+exp(y)))
  
  y <- param.best.df.wide$gamma +
    param.best.df.wide$psiT * x_val +
    param.best.df.wide$'psiDPI[2]' +
    param.best.df.wide$'psiAGE[2]' +
    param.best.df.wide$'psiSP[2]' +
    param.best.df.wide$psiASSAY * 1 +
    param.best.df.wide$'psiTG[1]'  +
    param.best.df.wide$psiDOSE *  5.5  +
    param.best.df.wide$'psiCELL[1]'
  
  y_cm <- median(exp(y)/(1+exp(y)))
  
  y <- param.best.df.wide$gamma +
    param.best.df.wide$psiT * x_val +
    param.best.df.wide$'psiDPI[2]' +
    param.best.df.wide$'psiAGE[2]' +
    param.best.df.wide$'psiSP[3]' +
    param.best.df.wide$psiASSAY * 1 +
    param.best.df.wide$'psiTG[1]'  +
    param.best.df.wide$psiDOSE *  5.5  +
    param.best.df.wide$'psiCELL[1]'
  
  y_agm <- median(exp(y)/(1+exp(y)))
  
  sp.df <- rbind(sp.df, 
                  data.frame(tRNA = x_val, 
                             sp = "RM", 
                             probs = y_rm),
                  data.frame(tRNA = x_val, 
                             sp = "CM",
                             probs = y_cm),
                  data.frame(tRNA = x_val, 
                             sp = "AGM",
                             probs = y_agm))
}

sp.df <- sp.df[-1, ]

fig.sp <- ggplot(sp.df, aes(x = as.numeric(tRNA), 
                              y = factor(sp, levels = rev(c("RM", "CM", "AGM"))),
                              fill = as.numeric(probs))) +
  geom_tile() + 
  geom_segment(x = -2, xend = -2, y = -4, yend = 3.2, linewidth = 0.2) +
  geom_rect(aes(xmin = 7.5, xmax = 8.5, ymin = 0.5, ymax = 3.5), 
            fill = "transparent", color = "grey33", size = 1.5) +
  geom_text(aes(label = (round(as.numeric(probs), 2)) * 100), color = "black", size = 2.5) +
  scale_fill_gradientn(colours = rev(brewer_pal), limits = c(0, 1)) +
  coord_cartesian(xlim = c(-0.5, 12.5), ylim = c(0.5, 3.5), clip = "off") +
  scale_x_continuous(breaks = seq(-0.5, 12.5, 2), 
                     labels = seq(0, 12, 2),
                     expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "log10 total RNA copies / sample", fill = "P(Culture Positive)",
       y = "Species") +
  theme(legend.position = "none",
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 8),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(angle = 45, size = 6.5,
                                   face = c("plain", "plain", "bold"), 
                                   vjust = 1, hjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(2, 2, 2, 2)); fig.sp



# 4C: Heat-map Exposure Conditions ---------------------------------------------

## Dose -----------------------------------------------------------------------

dose.df <- data.frame(tRNA = NA,
                      dose = NA,
                      probs = NA)

for (x_val in x_vals) {
  y <- param.best.df.wide$gamma +
    param.best.df.wide$psiT * x_val +
    param.best.df.wide$'psiDPI[2]' +
    param.best.df.wide$'psiAGE[2]' +
    param.best.df.wide$'psiSP[1]' +
    param.best.df.wide$psiASSAY * 1 +
    param.best.df.wide$'psiTG[1]'  +
    param.best.df.wide$psiDOSE *  4 +
    param.best.df.wide$'psiCELL[1]'
  
  y_small <- median(exp(y)/(1+exp(y)))
  
  y <- param.best.df.wide$gamma +
    param.best.df.wide$psiT * x_val +
    param.best.df.wide$'psiDPI[2]' +
    param.best.df.wide$'psiAGE[2]' +
    param.best.df.wide$'psiSP[1]' +
    param.best.df.wide$psiASSAY * 1 +
    param.best.df.wide$'psiTG[1]'  +
    param.best.df.wide$psiDOSE *  5.5 +
    param.best.df.wide$'psiCELL[1]'
  
  y_mid <- median(exp(y)/(1+exp(y)))
  
  y <- param.best.df.wide$gamma +
    param.best.df.wide$psiT * x_val +
    param.best.df.wide$'psiDPI[2]' +
    param.best.df.wide$'psiAGE[2]' +
    param.best.df.wide$'psiSP[1]' +
    param.best.df.wide$psiASSAY * 1 +
    param.best.df.wide$'psiTG[1]'  +
    param.best.df.wide$psiDOSE *  7 +
    param.best.df.wide$'psiCELL[1]'
  
  y_big <- median(exp(y)/(1+exp(y)))
  
  dose.df <- rbind(dose.df, 
                    data.frame(tRNA = x_val, 
                               dose = "4", 
                               probs = y_small),
                    data.frame(tRNA = x_val, 
                               dose = "5.5",
                               probs = y_mid),
                    data.frame(tRNA = x_val, 
                              dose = "7",
                              probs = y_big))
}

dose.df <- dose.df[-1, ]

fig.dose <- ggplot(dose.df, aes(x = as.numeric(tRNA), 
                                y = factor(dose, levels = c("7", "5.5", "4")), 
                                fill = as.numeric(probs))) +
  geom_tile() + 
  geom_segment(x = -2, xend = -2, y = 0.5, yend = 3.2, linewidth = 0.2) +
  geom_text(aes(label = (round(as.numeric(probs), 2)) * 100), color = "black", size = 2.5) +
  geom_rect(aes(xmin = 5.5, xmax = 8.5, ymin = 0.5, ymax = 3.5), 
            fill = "transparent", color = "grey33", size = 1.5) +
  scale_fill_gradientn(colours = rev(brewer_pal), limits = c(0, 1)) +
  coord_cartesian(xlim = c(-0.5, 12.5), ylim = c(0.5, 3.5), clip = "off") +
  scale_x_continuous(breaks = seq(-0.5, 12.5, 2), 
                     labels = seq(0, 12, 2),
                     expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "log10 total RNA copies / sample", fill = "P(Culture Positive)",
       y = "log10 Dose") +
  theme(legend.position = "none",
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 8),
        axis.text.y = element_text(angle = 45, size = 6.5, face = c("plain", "bold", "plain"), 
                                   vjust = 1, hjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(2, 2, 2, 2)); fig.dose


## DPI -----------------------------------------------------------------------

dpi.df <- data.frame(tRNA = NA,
                      dpi = NA,
                      probs = NA)

for (x_val in x_vals) {
  y <- param.best.df.wide$gamma +
    param.best.df.wide$psiT * x_val +
    param.best.df.wide$'psiDPI[1]' +
    param.best.df.wide$'psiAGE[2]' +
    param.best.df.wide$'psiSP[1]' +
    param.best.df.wide$psiASSAY * 1 +
    param.best.df.wide$'psiTG[1]'  +
    param.best.df.wide$psiDOSE * 5.5 +
    param.best.df.wide$'psiCELL[1]'
  
  y_dpi1 <- median(exp(y)/(1+exp(y)))
  
  y <- param.best.df.wide$gamma +
    param.best.df.wide$psiT * x_val +
    param.best.df.wide$'psiDPI[2]' +
    param.best.df.wide$'psiAGE[2]' +
    param.best.df.wide$'psiSP[1]' +
    param.best.df.wide$psiASSAY * 1 +
    param.best.df.wide$'psiTG[1]'  +
    param.best.df.wide$psiDOSE * 5.5 +
    param.best.df.wide$'psiCELL[1]'
  
  y_dpi2 <- median(exp(y)/(1+exp(y)))
  
  y <- param.best.df.wide$gamma +
    param.best.df.wide$psiT * x_val +
    param.best.df.wide$'psiDPI[3]' +
    param.best.df.wide$'psiAGE[2]' +
    param.best.df.wide$'psiSP[1]' +
    param.best.df.wide$psiASSAY * 1 +
    param.best.df.wide$'psiTG[1]'  +
    param.best.df.wide$psiDOSE *  5.5 +
    param.best.df.wide$'psiCELL[1]'
  
  y_dpi3 <- median(exp(y)/(1+exp(y)))
  
  dpi.df <- rbind(dpi.df, 
                   data.frame(tRNA = x_val, 
                              dpi = "I, 1", 
                              probs = y_dpi1),
                   data.frame(tRNA = x_val, 
                              dpi = "I, 2+",
                              probs = y_dpi2),
                   data.frame(tRNA = x_val, 
                              dpi = "NI, 1+",
                              probs = y_dpi3))
}

dpi.df <- dpi.df[-1, ]

fig.dpi <- ggplot(dpi.df, aes(x = as.numeric(tRNA), 
                              y = factor(dpi, levels = rev(c("I, 1", "I, 2+", "NI, 1+"))), 
                              fill = as.numeric(probs))) +
  geom_tile() + 
  geom_segment(x = -2, xend = -2, y = 0.7, yend = 6, linewidth = 0.2) +
  geom_text(aes(label = (round(as.numeric(probs), 2)) * 100), color = "black", size = 2.5) +
  geom_rect(aes(xmin = 5.5, xmax = 7.5, ymin = 0.5, ymax = 3.5), 
            fill = "transparent", color = "grey33", size = 1.5) +
  scale_fill_gradientn(colours = rev(brewer_pal), limits = c(0, 1)) +
  coord_cartesian(xlim = c(-0.5, 12.5), ylim = c(0.5, 3.5), clip = "off") +
  scale_x_continuous(breaks = seq(-0.5, 12.5, 2), 
                     labels = seq(0, 12, 2),
                     expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "log10 total RNA copies / sample", fill = "P(Culture Positive)",
       y = "DPI") +
  theme(legend.position = "none",
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 8),
        axis.text.y = element_text(angle = 45, size = 6.5, face = c("plain", "bold", "plain"), 
                                   vjust = 1, hjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(2, 2, 2, 2)); fig.dpi

fig.cond <- fig.dose + fig.dpi + plot_layout(ncol = 1); fig.cond


# 4C: Heat-map Assay Stuff -----------------------------------------------------

## Cell ------------------------------------------------------------------------

cell.df <- data.frame(tRNA = NA,
                      cell = NA,
                      probs = NA)

for (x_val in x_vals) {
  y <- param.best.df.wide$gamma +
    param.best.df.wide$psiT * x_val +
    param.best.df.wide$'psiDPI[2]' +
    param.best.df.wide$'psiAGE[2]' +
    param.best.df.wide$'psiSP[1]' +
    param.best.df.wide$psiASSAY * 1 +
    param.best.df.wide$'psiTG[1]'  +
    param.best.df.wide$psiDOSE *  5.5 +
    param.best.df.wide$'psiCELL[1]'
  
  y1 <- median(exp(y)/(1+exp(y)))
  
  y <- param.best.df.wide$gamma +
    param.best.df.wide$psiT * x_val +
    param.best.df.wide$'psiDPI[2]' +
    param.best.df.wide$'psiAGE[2]' +
    param.best.df.wide$'psiSP[1]' +
    param.best.df.wide$psiASSAY * 1 +
    param.best.df.wide$'psiTG[1]'  +
    param.best.df.wide$psiDOSE * 5.5 +
    param.best.df.wide$'psiCELL[2]'
  
  y2 <- median(exp(y)/(1+exp(y)))
  
  y <- param.best.df.wide$gamma +
    param.best.df.wide$psiT * x_val +
    param.best.df.wide$'psiDPI[2]' +
    param.best.df.wide$'psiAGE[2]' +
    param.best.df.wide$'psiSP[1]' +
    param.best.df.wide$psiASSAY * 1 +
    param.best.df.wide$'psiTG[1]'  +
    param.best.df.wide$psiDOSE *  5.5 +
    param.best.df.wide$'psiCELL[3]'
  
  y3 <- median(exp(y)/(1+exp(y)))
  
  cell.df <- rbind(cell.df, 
                  data.frame(tRNA = x_val, 
                             cell = "E6", 
                             probs = y1),
                  data.frame(tRNA = x_val, 
                             cell = "76",
                             probs = y2),
                  data.frame(tRNA = x_val, 
                             cell = "E6-SS2",
                             probs = y3))
}

cell.df <- cell.df[-1, ]

fig.cell <- ggplot(cell.df, aes(x = as.numeric(tRNA), 
                                y = factor(cell, levels = rev(c("76", "E6", "E6-SS2"))), 
                                fill = as.numeric(probs))) +
  geom_tile() + 
  geom_segment(x = -2, xend = -2, y = 0.5, yend = 3.2, linewidth = 0.2) +
  geom_text(aes(label = (round(as.numeric(probs), 2)) * 100), color = "black", size = 2.5) +
  geom_rect(aes(xmin = 5.5, xmax = 8.5, ymin = 0.5, ymax = 3.5), 
            fill = "transparent", color = "grey33", size = 1.5) +
  scale_fill_gradientn(colours = rev(brewer_pal), limits = c(0, 1)) +
  coord_cartesian(xlim = c(-0.5, 12.5), ylim = c(0.5, 3.5), clip = "off") +
  scale_x_continuous(breaks = seq(-0.5, 12.5, 2), 
                     labels = seq(0, 12, 2),
                     expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "log10 total RNA copies / sample", fill = "P(Culture Positive)",
       y = "Cell Line") +
  theme(legend.position = "none",
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.text.y = element_text(angle = 45, size = 6.5, face = c("plain", "bold", "plain"), 
                                   vjust = 1, hjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(2, 2, 2, 2)); fig.cell


## Assay -----------------------------------------------------------------------

assay.df <- data.frame(tRNA = NA,
                       assay = NA,
                       probs = NA)

for (x_val in x_vals) {
  y <- param.best.df.wide$gamma +
    param.best.df.wide$psiT * x_val +
    param.best.df.wide$'psiDPI[2]' +
    param.best.df.wide$'psiAGE[2]' +
    param.best.df.wide$'psiSP[1]' +
    param.best.df.wide$psiASSAY * 0 +
    param.best.df.wide$'psiTG[1]'  +
    param.best.df.wide$psiDOSE *  5.5 +
    param.best.df.wide$'psiCELL[1]'
  
  y1 <- median(exp(y)/(1+exp(y)))
  
  y <- param.best.df.wide$gamma +
    param.best.df.wide$psiT * x_val +
    param.best.df.wide$'psiDPI[2]' +
    param.best.df.wide$'psiAGE[2]' +
    param.best.df.wide$'psiSP[1]' +
    param.best.df.wide$psiASSAY * 1 +
    param.best.df.wide$'psiTG[1]'  +
    param.best.df.wide$psiDOSE * 5.5 +
    param.best.df.wide$'psiCELL[1]'
  
  y2 <- median(exp(y)/(1+exp(y)))
  
  
  assay.df <- rbind(assay.df, 
                    data.frame(tRNA = x_val, 
                               assay = "TCID50", 
                               probs = y1),
                    data.frame(tRNA = x_val, 
                               assay = "Plaque",
                               probs = y2))
}

assay.df <- assay.df[-1, ]

fig.assay <- ggplot(assay.df, aes(x = as.numeric(tRNA), y = assay, fill = as.numeric(probs))) +
  geom_tile() + 
  geom_segment(x = -2, xend = -2, y = -3, yend = 4, linewidth = 0.2) +
  geom_text(aes(label = (round(as.numeric(probs), 2)) * 100), color = "black", size = 2.5) +
  geom_rect(aes(xmin = 4.5, xmax = 7.5, ymin = 0.5, ymax = 2.5), 
            fill = "transparent", color = "grey33", size = 1.5) +
  scale_fill_gradientn(colours = rev(brewer_pal), limits = c(0, 1)) +
  coord_cartesian(xlim = c(-0.5, 12.5), ylim = c(0.5, 2.5), clip = "off") +
  scale_x_continuous(breaks = seq(-0.5, 12.5, 2), 
                     labels = seq(0, 12, 2),
                     expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "log10 total RNA copies / sample", fill = "P(Culture Positive)",
       y = "Assay") +
  theme(legend.position = "none",
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 8),
        axis.text.y = element_text(angle = 45, size = 6.5, face = c("bold", "plain"), 
                                   vjust = 1, hjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(2, 2, 2, 2)); fig.assay



## TG -----------------------------------------------------------------------

tg.df <- data.frame(tRNA = NA,
                      tg = NA,
                      probs = NA)

for (x_val in x_vals) {
  y <- param.best.df.wide$gamma +
    param.best.df.wide$psiT * x_val +
    param.best.df.wide$'psiDPI[2]' +
    param.best.df.wide$'psiAGE[2]' +
    param.best.df.wide$'psiSP[1]' +
    param.best.df.wide$psiASSAY * 1 +
    param.best.df.wide$'psiTG[1]'  +
    param.best.df.wide$psiDOSE *  5.5 +
    param.best.df.wide$'psiCELL[1]'
  
  y1 <- median(exp(y)/(1+exp(y)))
  
  y <- param.best.df.wide$gamma +
    param.best.df.wide$psiT * x_val +
    param.best.df.wide$'psiDPI[2]' +
    param.best.df.wide$'psiAGE[2]' +
    param.best.df.wide$'psiSP[1]' +
    param.best.df.wide$psiASSAY * 1 +
    param.best.df.wide$'psiTG[2]'  +
    param.best.df.wide$psiDOSE * 5.5 +
    param.best.df.wide$'psiCELL[1]'
  
  y2 <- median(exp(y)/(1+exp(y)))
  
  y <- param.best.df.wide$gamma +
    param.best.df.wide$psiT * x_val +
    param.best.df.wide$'psiDPI[2]' +
    param.best.df.wide$'psiAGE[2]' +
    param.best.df.wide$'psiSP[1]' +
    param.best.df.wide$psiASSAY * 1 +
    param.best.df.wide$'psiTG[3]'  +
    param.best.df.wide$psiDOSE *  5.5 +
    param.best.df.wide$'psiCELL[1]'
  
  y3 <- median(exp(y)/(1+exp(y)))
  
  tg.df <- rbind(tg.df, 
                   data.frame(tRNA = x_val, 
                              tg = "N", 
                              probs = y1),
                   data.frame(tRNA = x_val, 
                              tg = "E",
                              probs = y2),
                   data.frame(tRNA = x_val, 
                              tg = "S",
                              probs = y3))
}

tg.df <- tg.df[-1, ]

fig.tg <- ggplot(tg.df, aes(x = as.numeric(tRNA), 
                            y = factor(tg, levels = rev(c("N", "E", "S"))), 
                            fill = as.numeric(probs))) +
  geom_tile() + 
  geom_segment(x = -2, xend = -2, y = 0.8, yend = 6, linewidth = 0.2) +
  geom_text(aes(label = (round(as.numeric(probs), 2)) * 100), color = "black", size = 2.5) +
  geom_rect(aes(xmin = 6.5, xmax = 8.5, ymin = 0.5, ymax = 3.5), 
            fill = "transparent", color = "grey33", size = 1.5) +
  scale_fill_gradientn(colours = rev(brewer_pal), limits = c(0, 1),
                       breaks = c(0, 0.5, 1), labels = c(0, "50", "100")) +
  coord_cartesian(xlim = c(-0.5, 12.5), ylim = c(0.5, 3.5), clip = "off") +
  scale_x_continuous(breaks = seq(-0.5, 12.5, 2), 
                     labels = seq(0, 12, 2),
                     expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  guides(fill = guide_colorbar(frame.colour = "black",
                                title.vjust = 0.7)) +
  labs(x = "log10 total RNA copies / sample", fill = "Chance of Culture (%)",
       y = "Target Gene") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 7.5, vjust = 0.7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        legend.margin = margin(t = -1),
        axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(angle = 45, size = 6.5, face = c("plain", "plain", "bold"), 
                                   vjust = 1, hjust = 0.5),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(2, 2, 2, 2)); fig.tg


# 4C: Combine into a Heat-map --------------------------------------------------

fig4C <- fig.simple + labs(tag = "C") + fig.dose + fig.dpi + 
          fig.sp + fig.age + 
          fig.cell + fig.assay + fig.tg +
          plot_layout(ncol = 1, heights = c(1, 3, 3, 3, 3, 3, 2, 3)); fig4C


# 4D: Fit lines ----------------------------------------------------------------

## Set number of lines 
n_lines <- 300

## Get simple model results ----------------------------------------------------

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

# Extract the mean for the reference line
psi_mean <- mean(param.simple.df.wide$psiT)
gamma_mean <- mean(param.simple.df.wide$gamma)

# Generate values & a data frame for the reference line
x_vals <- seq(from = 0, to = 12.5, by = 0.05)
y_vals <- exp(gamma_mean + psi_mean * x_vals)/(1 + exp(gamma_mean + psi_mean * x_vals))
lines.simple.df <- cbind(x_vals, y_vals)



## Dose-panel trajectories -----------------------------------------------------

assay <- 1
sex <- 1
dose <- 5.5


### Best model: get sample lines -----------------------------------------------

# Extract lines for different dose values, holding all else constant
set.seed(1)
dose.vary <- sample(runif(n_lines, min(dat$log10_dose_pfu), max(dat$log10_dose_pfu)))

lines.best.df <- as.data.frame(matrix(nrow = 0, ncol = 4))

for (ii in 1:length(dose.vary)) {
  row <- sample(1:nrow(param.best.df.wide), size = 1)
  y <- param.best.df.wide$gamma[row] +
    param.best.df.wide$psiT[row] * x_vals +
    param.best.df.wide$'psiDPI[2]'[row] +
    param.best.df.wide$'psiAGE[2]'[row] +
    param.best.df.wide$'psiSP[1]'[row] +
    param.best.df.wide$psiASSAY[row] * assay +
    param.best.df.wide$'psiTG[1]'[row]  +
    param.best.df.wide$psiDOSE[row] * dose.vary[ii] +
    param.best.df.wide$'psiCELL[1]'[row]  
  
  y_vals <- exp(y)/(1 + exp(y))
  line_num <- rep(ii, length(x_vals))
  log10_dose <- rep(dose.vary[ii], length(x_vals))
  
  lines.best.df <- rbind(lines.best.df,
                         cbind(line_num, log10_dose, x_vals, y_vals))
}


### Plot ------------------------------------------------------------------------

fig.dose.fits <- ggplot() + 
  geom_line(data = lines.best.df,
            aes(x = x_vals, y = y_vals, 
                group = line_num,
                color = log10_dose),
            size = 0.8, alpha = 0.7) +
  geom_line(data = as.data.frame(lines.simple.df),
            aes(x = x_vals, y = y_vals),
            color = "#0F4C5E",
            size = 1) +
  scale_color_gradientn(colours = zissou_pal) +
  scale_fill_gradientn(colours = zissou_pal) +
  labs(y = "Chance of Culture Positive (%)", 
       x = "log10 total RNA copies / sample",
       color = "log10\nDose") +
  guides(fill = "none",
         color = guide_colorbar(frame.colour = "black",
                                title.vjust = 0)) +
  theme(legend.position = c(0.11, 0.665),
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(2, 2, 0, 2)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1),
                     labels = c("", "20", "40", "60", "80", 100), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12), limits = c(-0.01, 12.5),
                     expand = c(0, 0)) +
  geom_hline(yintercept = 1, size = 0.6, color = "black") +
  geom_hline(yintercept = 0, size = 0.6, color = "black"); fig.dose.fits



## DPI-dependent trajectories --------------------------------------------------

### Prep & plot DPI results ----------------------------------------------------

set.seed(12)
dpi.vary <- sample(rep(1:3, each = n_lines / 3))

lines.dpi.df <- as.data.frame(matrix(nrow = 0, ncol = 4))

for (ii in 1:length(dpi.vary)){
  row <- sample(1:nrow(param.best.df.wide), size = 1)
  
  if (dpi.vary[ii] == 1){
    dpi_sample <- param.best.df.wide$'psiDPI[1]'[row]
  }
  else if (dpi.vary[ii] == 2){
    dpi_sample <- param.best.df.wide$'psiDPI[2]'[row]
  }
  else if (dpi.vary[ii] == 3){
    dpi_sample <- param.best.df.wide$'psiDPI[3]'[row]
  }
  
  y <- param.best.df.wide$gamma[row] +
    param.best.df.wide$psiT[row] * x_vals +
    dpi_sample +
    param.best.df.wide$'psiAGE[2]'[row] +
    param.best.df.wide$'psiSP[1]'[row] +
    param.best.df.wide$psiASSAY[row] * assay +
    param.best.df.wide$'psiTG[1]'[row]  +
    param.best.df.wide$psiDOSE[row] * dose +
    param.best.df.wide$'psiCELL[1]'[row] 
  
  y_vals <- exp(y)/(1 + exp(y))
  line_num <- rep(ii, length(x_vals))
  dpi_idx <- rep(dpi.vary[ii], length(x_vals))
  
  lines.dpi.df <- rbind(lines.dpi.df,
                        cbind(line_num, dpi_idx, x_vals, y_vals))
}

lines.dpi.df$dpi_idx[lines.dpi.df$dpi_idx == 1] <- "I, 1"
lines.dpi.df$dpi_idx[lines.dpi.df$dpi_idx == 2] <- "I, 2+"
lines.dpi.df$dpi_idx[lines.dpi.df$dpi_idx == 3] <- "NI, 1+"


fig.dpi.fits <- ggplot(data = lines.dpi.df) + 
  geom_line(aes(x = x_vals, y = y_vals, 
                group = line_num,
                color = dpi_idx),
            size = 0.8, alpha = 0.3) +
  geom_line(data = as.data.frame(lines.simple.df),
            aes(x = x_vals, y = y_vals),
            color = "#002C3D",
            size = 1) +
  scale_color_manual(values = c("I, 1" = "#37123C",
                                "I, 2+" = "#A83869",
                                "NI, 1+" = "#D986B0")) +
  scale_fill_manual(values = c("1" = "#37123C",
                               "2" = "#A83869",
                               "3" = "#D986B0")) +
  labs(y = "", x = "",
       color = "DPI") +
  guides(color = guide_legend(override.aes = list(size = 2.7, alpha = 1, 
                                                  linewidth = 2), 
                              keyheight = 0.6),
         fill = "none") + #,
  theme(legend.position = c(0.15, 0.81),
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 7),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(2, 2, 0, 2)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1),
                     labels = c(0, 20, 40, 60, 80, 100), expand = c(0, 0),
                     position = "right") +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12), limits = c(-0.01, 12.5),
                     expand = c(0, 0)); fig.dpi.fits



## Age-dependent trajectories --------------------------------------------------

### Prep & plot AGE results ----------------------------------------------------

set.seed(121)
age.vary <- sample(rep(1:3, each = n_lines / 3))

lines.age.df <- as.data.frame(matrix(nrow = 0, ncol = 4))

for (ii in 1:length(age.vary)){
  row <- sample(1:nrow(param.best.df.wide), size = 1)
  
  if (age.vary[ii] == 1){
    age_sample <- param.best.df.wide$'psiAGE[1]'[row]
  }
  else if (age.vary[ii] == 2){
    age_sample <- param.best.df.wide$'psiAGE[2]'[row]
  }
  else if (age.vary[ii] == 3){
    age_sample <- param.best.df.wide$'psiAGE[3]'[row]
  }
  
  y <- param.best.df.wide$gamma[row] +
    param.best.df.wide$psiT[row] * x_vals +
    param.best.df.wide$'psiDPI[2]'[row] +
    age_sample +
    param.best.df.wide$'psiSP[1]'[row] +
    param.best.df.wide$psiASSAY[row] * assay +
    param.best.df.wide$'psiTG[1]'[row]  +
    param.best.df.wide$psiDOSE[row] * dose +
    param.best.df.wide$'psiCELL[1]'[row]  
  
  y_vals <- exp(y)/(1 + exp(y))
  line_num <- rep(ii, length(x_vals))
  age_idx <- rep(age.vary[ii], length(x_vals))
  
  lines.age.df <- rbind(lines.age.df,
                        cbind(line_num, age_idx, x_vals, y_vals))
}

lines.age.df$age_idx[lines.age.df$age_idx == 1] <- "Juvenile"
lines.age.df$age_idx[lines.age.df$age_idx == 2] <- "Adult"
lines.age.df$age_idx[lines.age.df$age_idx == 3] <- "Geriatric"


fig.age.fits <- ggplot(data = lines.age.df) + 
  annotate("rect", xmin = 0.04, xmax = 4.3, ymin = 0.6, ymax = 0.99, 
           fill = "white") +
  geom_line(aes(x = x_vals, y = y_vals, 
                group = line_num,
                color = factor(age_idx, levels = c("Juvenile", "Adult", "Geriatric"))),
            size = 0.8, alpha = 0.3) +
  geom_line(data = as.data.frame(lines.simple.df),
            aes(x = x_vals, y = y_vals),
            color = "#002C3D",
            size = 1) +
  scale_color_manual(values = c("Juvenile" = "#D1320F",
                                "Adult" =  "#980B0B", 
                                "Geriatric" = "#B3ADA2")) +
  labs(y = "", x = "",
       color = "Age Class") +
  guides(color = guide_legend(override.aes = list(size = 2.7, alpha = 1,
                                                  linewidth = 2), 
                              keyheight = 0.6),
         fill = "none") + #,
  theme(legend.position = c(0.181, 0.81),
        legend.title = element_text(size = 8.5),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(2, 2, 0, 2)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1),
                     labels = c(0, 20, 40, 60, 80, 100), expand = c(0, 0),
                     position = "right") +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10), limits = c(-0.01, 12.5),
                     expand = c(0, 0)); fig.age.fits



## Species-dependent trajectories ----------------------------------------------

### Prep & plot DPI results ----------------------------------------------------

set.seed(1212)
sp.vary <- sample(rep(1:3, each = n_lines / 3))

lines.sp.df <- as.data.frame(matrix(nrow = 0, ncol = 4))

for (ii in 1:length(sp.vary)){
  row <- sample(1:nrow(param.best.df.wide), size = 1)
  
  if (sp.vary[ii] == 1){
    sp_sample <- param.best.df.wide$'psiSP[1]'[row]
  }
  else if (sp.vary[ii] == 2){
    sp_sample <- param.best.df.wide$'psiSP[2]'[row]
  }
  else if (sp.vary[ii] == 3){
    sp_sample <- param.best.df.wide$'psiSP[3]'[row]
  }
  
  y <- param.best.df.wide$gamma[row] +
    param.best.df.wide$psiT[row] * x_vals +
    param.best.df.wide$'psiDPI[2]'[row] +
    param.best.df.wide$'psiAGE[2]'[row] +
    sp_sample +
    param.best.df.wide$psiASSAY[row] * assay +
    param.best.df.wide$'psiTG[1]'[row]  +
    param.best.df.wide$psiDOSE[row] * dose +
    param.best.df.wide$'psiCELL[1]'[row] 
  
  y_vals <- exp(y)/(1 + exp(y))
  line_num <- rep(ii, length(x_vals))
  sp_idx <- rep(sp.vary[ii], length(x_vals))
  
  lines.sp.df <- rbind(lines.sp.df,
                       cbind(line_num, sp_idx, x_vals, y_vals))
}

lines.sp.df$sp_idx[lines.sp.df$sp_idx == 1] <- "RM"
lines.sp.df$sp_idx[lines.sp.df$sp_idx == 2] <- "CM"
lines.sp.df$sp_idx[lines.sp.df$sp_idx == 3] <- "AGM"


fig.sp.fits <- ggplot(data = lines.sp.df) + 
  geom_line(aes(x = x_vals, y = y_vals, 
                group = line_num,
                color = factor(sp_idx, levels = c("RM", "CM", "AGM"))),
            size = 0.8, alpha = 0.3) +
  geom_line(data = as.data.frame(lines.simple.df),
            aes(x = x_vals, y = y_vals),
            color = "#002C3D",
            size = 1) +
  scale_color_manual(values = c("RM" = "#0B6884",
                                "CM" = "#744468",
                                "AGM" = "#D56540")) +
  labs(y = "", x = "log10 total RNA copies / sample", color = "Species") +
  guides(color = guide_legend(override.aes = list(size = 2.7, alpha = 1,
                                                  linewidth = 2), 
                              keyheight = 0.6),
         fill = "none") + #,
  theme(legend.position = c(0.145, 0.81),
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title.y = element_blank() ,
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(2, 2, 0, 2)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1),
                     labels = c("", "20", "40", "60", "80", 100), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10), limits = c(-0.01, 12.5),
                     expand = c(0, 0)); fig.sp.fits


## Assay-dependent trajectories -----------------------------------------------

set.seed(12121)
assay.vary <- sample(rep(0:1, each = n_lines / 2))

lines.assay.df <- as.data.frame(matrix(nrow = 0, ncol = 4))

for (ii in 1:length(assay.vary)){
  row <- sample(1:nrow(param.best.df.wide), size = 1)
  
  y <- param.best.df.wide$gamma[row] +
    param.best.df.wide$psiT[row] * x_vals +
    param.best.df.wide$'psiDPI[2]'[row] +
    param.best.df.wide$'psiAGE[2]'[row] +
    param.best.df.wide$'psiSP[1]'[row] +
    param.best.df.wide$psiASSAY[row] * assay.vary[ii] +
    param.best.df.wide$'psiTG[1]'[row]  +
    param.best.df.wide$psiDOSE[row] * dose +
    param.best.df.wide$'psiCELL[1]'[row]  
  
  y_vals <- exp(y)/(1 + exp(y))
  line_num <- rep(ii, length(x_vals))
  assay_idx <- rep(assay.vary[ii], length(x_vals))
  
  lines.assay.df <- rbind(lines.assay.df,
                          cbind(line_num, assay_idx, x_vals, y_vals))
}

lines.assay.df$assay_idx[lines.assay.df$assay_idx == 0] <- "TCID50"
lines.assay.df$assay_idx[lines.assay.df$assay_idx == 1] <- "Plaque"

### Plot -----------------------------------------------------------------------

fig.assay.fits <- ggplot(data = lines.assay.df) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = 1.2, 
           alpha = 0.2) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -0.2, ymax = 0, 
           alpha = 0.2) +
  geom_line(aes(x = x_vals, y = y_vals, 
                group = line_num,
                color = factor(assay_idx, levels = c("TCID50", "Plaque"))),
            size = 0.8, alpha = 0.3) +
  geom_line(data = as.data.frame(lines.simple.df),
            aes(x = x_vals, y = y_vals),
            color = "#002C3D",
            size = 1) +
  scale_color_manual(values = c("TCID50" = "#E38900",
                                "Plaque" = zissou_pal[5])) + 
  labs(y = "", color = "Assay") +
  guides(color = guide_legend(override.aes = list(size = 2.7, alpha = 1,
                                                  linewidth = 2), 
                              keyheight = 0.6),
         fill = "none") + 
  theme(legend.position = c(0.176, 0.84),
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 7),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(2, 2, 0, 2)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1),
                     labels = c(0, "20", "40", "60", "80", 100), expand = c(0, 0),
                     position = "right") +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12), limits = c(-0.01, 12.5),
                     expand = c(0, 0)); fig.assay.fits



## TG-dependent trajectories ---------------------------------------------------

### Prep results ---------------------------------------------------------------

set.seed(121212)
tg.vary <- sample(rep(1:3, each = n_lines / 3))

lines.tg.df <- as.data.frame(matrix(nrow = 0, ncol = 4))

for (ii in 1:length(tg.vary)){
  row <- sample(1:nrow(param.best.df.wide), size = 1)
  
  if (tg.vary[ii] == 1){
    tg_sample <- param.best.df.wide$'psiTG[1]'[row]
  }
  else if (tg.vary[ii] == 2){
    tg_sample <- param.best.df.wide$'psiTG[2]'[row]
  }
  else if (tg.vary[ii] == 3){
    tg_sample <- param.best.df.wide$'psiTG[3]'[row]
  }
  
  y <- param.best.df.wide$gamma[row] +
    param.best.df.wide$psiT[row] * x_vals +
    param.best.df.wide$'psiDPI[2]'[row] +
    param.best.df.wide$'psiAGE[2]'[row] +
    param.best.df.wide$'psiSP[1]'[row] +
    param.best.df.wide$psiASSAY[row] * assay +
    tg_sample  +
    param.best.df.wide$psiDOSE[row] * dose +
    param.best.df.wide$'psiCELL[1]'[row]  
  
  y_vals <- exp(y)/(1 + exp(y))
  line_num <- rep(ii, length(x_vals))
  tg_idx <- rep(tg.vary[ii], length(x_vals))
  
  lines.tg.df <- rbind(lines.tg.df,
                       cbind(line_num, tg_idx, x_vals, y_vals))
}

lines.tg.df$tg_idx[lines.tg.df$tg_idx == 1] <- "N"
lines.tg.df$tg_idx[lines.tg.df$tg_idx == 2] <- "E"
lines.tg.df$tg_idx[lines.tg.df$tg_idx == 3] <- "S"


### Plot -----------------------------------------------------------------------

fig.tg.fits <- ggplot(data = lines.tg.df) + 
  geom_line(aes(x = x_vals, y = y_vals, 
                group = line_num,
                color = factor(tg_idx, levels = c("N", "E", "S"))),
            size = 0.8, alpha = 0.3) +
  geom_line(data = as.data.frame(lines.simple.df),
            aes(x = x_vals, y = y_vals),
            color = "#002C3D",
            size = 1) +
  scale_color_manual(values = c("N" =  "#167358",  #"#0a9396",   #"#005f73",
                                "E" = "#27AEB0", #   "#089396",
                                "S" = "#DBBB5A")) +
  labs(y = "P(culture positive)", x = "log10 total RNA copies / sample", 
       color = "Target\nGene") +
  guides(color = guide_legend(override.aes = list(size = 2.7, alpha = 1,
                                                  linewidth = 2), 
                              keyheight = 0.6),
         fill = "none") + #,
  theme(legend.position = c(0.128, 0.77),
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 7),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(2, 2, 0, 2)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1),
                     labels = seq(0, 100, 20), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12), limits = c(-0.01, 12.5),
                     expand = c(0, 0)); fig.tg.fits



## Cell-dependent trajectories -------------------------------------------------

### Prep results ---------------------------------------------------------------

set.seed(1212121)
cell.vary <- sample(rep(1:3, each = n_lines / 3))

lines.cell.df <- as.data.frame(matrix(nrow = 0, ncol = 4))

for (ii in 1:length(cell.vary)){
  row <- sample(1:nrow(param.best.df.wide), size = 1)
  
  if (cell.vary[ii] == 1){
    cell_sample <- param.best.df.wide$'psiCELL[1]'[row]
  }
  else if (cell.vary[ii] == 2){
    cell_sample <- param.best.df.wide$'psiCELL[2]'[row]
  }
  else if (cell.vary[ii] == 3){
    cell_sample <- param.best.df.wide$'psiCELL[3]'[row]
  }
  
  y <- param.best.df.wide$gamma[row] +
    param.best.df.wide$psiT[row] * x_vals +
    param.best.df.wide$'psiDPI[2]'[row] +
    param.best.df.wide$'psiAGE[2]'[row] +
    param.best.df.wide$'psiSP[1]'[row] +
    param.best.df.wide$psiASSAY[row] * assay +
    param.best.df.wide$'psiTG[1]'[row]  +
    param.best.df.wide$psiDOSE[row] * dose +
    cell_sample  
  
  y_vals <- exp(y)/(1 + exp(y))
  line_num <- rep(ii, length(x_vals))
  cell_idx <- rep(cell.vary[ii], length(x_vals))
  
  lines.cell.df <- rbind(lines.cell.df,
                         cbind(line_num, cell_idx, x_vals, y_vals))
}

lines.cell.df$cell_idx[lines.cell.df$cell_idx == 1] <- "E6"
lines.cell.df$cell_idx[lines.cell.df$cell_idx == 2] <- "76"
lines.cell.df$cell_idx[lines.cell.df$cell_idx == 3] <- "E6-SS2"


### Plot -----------------------------------------------------------------------

fig.cell.fits <- ggplot(data = lines.cell.df) + 
  geom_line(aes(x = x_vals, y = y_vals, 
                group = line_num,
                color = cell_idx),
            size = 0.8, alpha = 0.3) +
  geom_line(data = as.data.frame(lines.simple.df),
            aes(x = x_vals, y = y_vals),
            color = "#002C3D",
            size = 1) +
  scale_color_manual(values = c("76" = "#bc4749",   
                                "E6" = "#386641", 
                                "E6-SS2" = "#a7c957")) +  
  labs(y = "", x = "log10 total RNA copies / sample", color = "Cell Line") +
  guides(color = guide_legend(override.aes = list(size = 2.7, alpha = 1,
                                                  linewidth = 2), 
                              keyheight = 0.6),
         fill = "none") + 
  theme(legend.position = c(0.18, 0.81),
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title = element_blank() ,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 7),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(2, 2, 0, 2)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1),
                     labels = seq(0, 100, 20), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12), limits = c(-0.01, 12.5),
                     expand = c(0, 0)); fig.cell.fits


## Plot together  -------------------------------------------------------------

# Generate text-only figure for joint labels
y_title <- ggplot() + 
  annotate(geom = "text", x = 1, y = 1, size = 3.2,
           label = "Chance of Culture Positive (%)", angle = 90) +
  theme(plot.margin = c(0, 0, 0, 0)) +
  theme_void() +
  labs(tag = "D")

x_title <- ggplot() + 
  annotate(geom = "text", x = 1, y = 1.7, size = 2.8,
           label = "log10 total RNA copies / sample") +
  theme(plot.margin = c(-10, -10, -10, -10)) +
  scale_y_continuous(limits = c(0, 2)) +
  theme_void()


layout_matrix <- matrix(c(1, 1, 2, 2, 3, 3, 4, 4, 8, 5, 5, 6, 6, 7, 7, 8), nrow = 8)

fig4D <- grid.arrange(fig.dose.fits, fig.sp.fits, fig.cell.fits, fig.tg.fits,
             fig.dpi.fits, fig.age.fits, fig.assay.fits, layout_matrix = layout_matrix)



fig4D <- ggarrange(y_title, fig4D, ncol = 2, widths = c(0.045, 1)); fig4D
fig4D <- ggarrange(fig4D, x_title, nrow = 2, heights = c(1, 0.025)); fig4D

  

# Combine bottom panel  --------------------------------------------------------

fig.bottom <- ggarrange(fig4C, fig4D, ncol = 2, widths = c(1, 1.3)); fig.bottom


# Combine all of it  -----------------------------------------------------------

fig4 <- ggarrange(fig.top, fig.bottom, nrow = 2, heights = c(1.3, 4.5)); fig4

ggsave('./outputs/figures/fig5-culture-model-results.tiff',
       plot = fig4,
       device = 'tiff',
       height = 10.5,
       width = 9,
       units = 'in',
       bg = 'white')

