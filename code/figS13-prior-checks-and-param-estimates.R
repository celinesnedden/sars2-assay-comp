# Generates prior predictive simulations for both components of the hurdle model
# Includes the code for Supplemental Figure S13

library(extraDistr)
library(patchwork)
library(ggplot2)
library(wesanderson)
library(dplyr)
library(tidyr)
zissou_pal <- wes_palette("Zissou1", 12, type = "continuous")
set.seed(159)

# Prior predictive lines -------------------------------------------------------

## PCR Model --------------------------------------------------------------------

### Logistic model ---------------------------------------------------------------

# Generate simulated data
N_samples <- 200
t_sim <- seq(from = 0, to = 10, length.out = N_samples)
st_sim <- sample(x = c(1, 2), size = N_samples, replace = TRUE)
dose_sim <- seq(from = 4, to = 7.5, length.out = N_samples)
age_sim <- sample(x = 1:3, size = N_samples, replace = TRUE)
sex_sim <- sample(x = 1:2, size = N_samples, replace = TRUE)
sp_sim <- sample(x = 1:3, size = N_samples, replace = TRUE)
tg_sim <- sample(x = 1:4, size = N_samples, replace = TRUE)
dpi_sim <- sample(x = 1:3, size = N_samples, replace = TRUE)

# Generate indicator cols for TG types
tg_sim1 <- tg_sim
tg_sim1[tg_sim == 1] <- 1
tg_sim1[tg_sim != 1] <- 0

tg_sim2 <- tg_sim
tg_sim2[tg_sim == 2] <- 1
tg_sim2[tg_sim != 2] <- 0

tg_sim3 <- tg_sim
tg_sim3[tg_sim == 3] <- 1
tg_sim3[tg_sim != 3] <- 0

tg_sim4 <- tg_sim
tg_sim4[tg_sim == 4] <- 1
tg_sim4[tg_sim != 4] <- 0

# Generate indicator cols for DPI types
dpi_sim1 <- dpi_sim
dpi_sim1[dpi_sim == 1] <- 1
dpi_sim1[dpi_sim != 1] <- 0

dpi_sim2 <- dpi_sim
dpi_sim2[dpi_sim == 2] <- 1
dpi_sim2[dpi_sim != 2] <- 0

dpi_sim3 <- dpi_sim
dpi_sim3[dpi_sim == 3] <- 1
dpi_sim3[dpi_sim != 3] <- 0

# Generate 
prior_df <- as.data.frame(cbind(t_sim, st_sim, dose_sim, age_sim, sex_sim, sp_sim,
                                tg_sim, tg_sim1, tg_sim2, tg_sim3, tg_sim4,
                                dpi_sim, dpi_sim1, dpi_sim2, dpi_sim3))



#### Simulate coefficients -------------------------------------------------------

##### Informative ----------------------------------------------------------------

gamma_inf = rnorm(n = N_samples, mean = -1, sd = 2)
deltaT_inf = rnorm(n = N_samples, mean = 2, sd = 1)
deltaST_inf = rnorm(n = N_samples, mean = 0, sd = 1)
deltaDOSE_inf = rnorm(n = N_samples, mean = -1, sd = 0.5)
deltaAGE_inf = rnorm(n = N_samples, mean = 0, sd = 1)
deltaSEX_inf = rnorm(n = N_samples, mean = 0, sd = 1)
deltaSP_inf = rnorm(n = N_samples, mean = 0, sd = 1)
deltaTG1_inf = rnorm(n = N_samples, mean = 1, sd = 1)
deltaTG2_inf = rnorm(n = N_samples, mean = 1, sd = 1)
deltaTG3_inf = rnorm(n = N_samples, mean = -1, sd = 1)
deltaTG4_inf = rnorm(n = N_samples, mean = -1, sd = 1)
deltaDPI1_inf = rnorm(n = N_samples, mean = -1, sd = 1)
deltaDPI2_inf = rnorm(n = N_samples, mean = 1, sd = 1)
deltaDPI3_inf = rnorm(n = N_samples, mean = 0, sd = 1)


##### Non-Informative ------------------------------------------------------------

gamma_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
deltaT_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
deltaST_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
deltaDOSE_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
deltaAGE_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
deltaSEX_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
deltaSP_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
deltaTG1_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
deltaTG2_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
deltaTG3_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
deltaTG4_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
deltaDPI1_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
deltaDPI2_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
deltaDPI3_noninf = rnorm(n = N_samples, mean = 0, sd = 1)


#### Generate predictions -------------------------------------------------------

pred_df_inf <- as.data.frame(t_sim)
pred_df_noninf <- as.data.frame(t_sim)

for (row_num in 1:nrow(prior_df)) {
  y_inf <- gamma_inf[row_num] + deltaT_inf[row_num] * prior_df$t_sim + 
           deltaST_inf[row_num] * prior_df$st_sim[row_num] + 
           deltaDOSE_inf[row_num] * prior_df$dose_sim[row_num] + 
           deltaAGE_inf[row_num] * prior_df$age_sim[row_num] + 
           deltaSEX_inf[row_num] * prior_df$sex_sim[row_num] + 
           deltaSP_inf[row_num] * prior_df$sp_sim[row_num] +
           deltaTG1_inf[row_num] * prior_df$tg_sim1[row_num] +
           deltaTG2_inf[row_num] * prior_df$tg_sim2[row_num] +
           deltaTG3_inf[row_num] * prior_df$tg_sim3[row_num] +
           deltaTG4_inf[row_num] * prior_df$tg_sim4[row_num] +
           deltaDPI1_inf[row_num] * prior_df$dpi_sim1[row_num] +
           deltaDPI2_inf[row_num] * prior_df$dpi_sim2[row_num] +
           deltaDPI3_inf[row_num] * prior_df$dpi_sim3[row_num]
  y_noninf <- gamma_noninf[row_num] + deltaT_noninf[row_num] * prior_df$t_sim + 
              deltaST_noninf[row_num] * prior_df$st_sim[row_num] + 
              deltaDOSE_noninf[row_num] * prior_df$dose_sim[row_num] + 
              deltaAGE_noninf[row_num] * prior_df$age_sim[row_num] + 
              deltaSEX_noninf[row_num] * prior_df$sex_sim[row_num] + 
              deltaSP_noninf[row_num] * prior_df$sp_sim[row_num] +
              deltaTG1_noninf[row_num] * prior_df$tg_sim1[row_num] +
              deltaTG2_noninf[row_num] * prior_df$tg_sim2[row_num] +
              deltaTG3_noninf[row_num] * prior_df$tg_sim3[row_num] +
              deltaTG4_noninf[row_num] * prior_df$tg_sim4[row_num] +
              deltaDPI1_noninf[row_num] * prior_df$dpi_sim1[row_num] +
              deltaDPI2_noninf[row_num] * prior_df$dpi_sim2[row_num] +
              deltaDPI3_noninf[row_num] * prior_df$dpi_sim3[row_num]
  
  y_inf <- exp(y_inf)/(1 + exp(y_inf))
  y_noninf <- exp(y_noninf)/(1 + exp(y_noninf))
  pred_df_inf <- cbind(pred_df_inf, y_inf)
  pred_df_noninf <- cbind(pred_df_noninf, y_noninf)
}

# Convert colnames and pivot longer to maintain groups & plot
colnames(pred_df_inf) <- c("t_val", paste0("y", 1:N_samples))
pred_df_inf <- pred_df_inf %>% pivot_longer(cols = !c("t_val"),
                                            names_to = c("sample"))
pred_df_inf$prior <- "Informative"

colnames(pred_df_noninf) <- c("t_val", paste0("y", 1:N_samples))
pred_df_noninf <- pred_df_noninf %>% pivot_longer(cols = !c("t_val"),
                                                  names_to = c("sample"))
pred_df_noninf$prior <- "Non-informative"


# Combine into one dataframe for facetting
pred_df <- rbind(pred_df_inf, pred_df_noninf)


#### Make the plot  -------------------------------------------------------------

pred_df$y_facet <- "PCR Logistic"
figSXA <- ggplot(data = pred_df) +
  geom_line(aes(x = t_val, y = value, group = sample, color = prior), 
            alpha = 0.8, lwd = 0.5) +
  labs(y = "Chance of sgRNA Positive (%)", x = "log10 total RNA copies / sample", tag = "A") +
  facet_grid(y_facet ~ prior) +
  scale_color_manual(values = c("Informative" = zissou_pal[1],
                                "Non-informative" = zissou_pal[12])) +
  scale_x_continuous(limits = c(0, 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0),
                     breaks = c(0, 0.5, 1), 
                     labels = c("0", "50", "100")) +
  theme(legend.position = "none", 
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title.y = element_text(size = 9),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.background = element_rect(fill='white'),
        text = element_text(size = 12),
        strip.text = element_text(size = 9),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(size = 0.1, linetype = 'solid',
                                          colour = "light grey"),
        plot.margin = margin(2, 2, 2, 2)); figSXA
  


### Linear model ---------------------------------------------------------------

# Generate simulated data
N_samples <- 200
t_sim <- seq(from = 0, to = 10, length.out = N_samples)
st_sim <- sample(x = c(0, 1), size = N_samples, replace = TRUE)
dose_sim <- seq(from = 4, to = 7.5, length.out = N_samples)
age_sim <- sample(x = 1:3, size = N_samples, replace = TRUE)
sex_sim <- sample(x = 1:2, size = N_samples, replace = TRUE)
sp_sim <- sample(x = 1:3, size = N_samples, replace = TRUE)
tg_sim <- sample(x = 1:4, size = N_samples, replace = TRUE)
dpi_sim <- sample(x = 1:3, size = N_samples, replace = TRUE)

# Generate indicator cols for TG types
tg_sim1 <- tg_sim
tg_sim1[tg_sim == 1] <- 1
tg_sim1[tg_sim != 1] <- 0

tg_sim2 <- tg_sim
tg_sim2[tg_sim == 2] <- 1
tg_sim2[tg_sim != 2] <- 0

tg_sim3 <- tg_sim
tg_sim3[tg_sim == 3] <- 1
tg_sim3[tg_sim != 3] <- 0

tg_sim4 <- tg_sim
tg_sim4[tg_sim == 4] <- 1
tg_sim4[tg_sim != 4] <- 0

# Generate indicator cols for DPI types
dpi_sim1 <- dpi_sim
dpi_sim1[dpi_sim == 1] <- 1
dpi_sim1[dpi_sim != 1] <- 0

dpi_sim2 <- dpi_sim
dpi_sim2[dpi_sim == 2] <- 1
dpi_sim2[dpi_sim != 2] <- 0

dpi_sim3 <- dpi_sim
dpi_sim3[dpi_sim == 3] <- 1
dpi_sim3[dpi_sim != 3] <- 0

# Generate 
prior_df <- as.data.frame(cbind(t_sim, st_sim, dose_sim, age_sim, sex_sim, sp_sim,
                                tg_sim, tg_sim1, tg_sim2, tg_sim3,
                                dpi_sim, dpi_sim1, dpi_sim2, dpi_sim3))


#### Simulate coefficients -------------------------------------------------------

##### Informative ----------------------------------------------------------------

alpha_inf = rnorm(n = N_samples, mean = -3, sd = 1)

betaT_inf = rgamma(n = N_samples, shape = 2, scale = 0.5); #hist(betaT)
betaST_inf = rnorm(n = N_samples, mean = 0, sd = 1)
betaDOSE_inf = rnorm(n = N_samples, mean = -0.25, sd = 0.25); #hist(betaDOSE)
betaAGE_inf = rnorm(n = N_samples, mean = 0, sd = 0.1)
betaSEX_inf = rnorm(n = N_samples, mean = 0, sd = 0.1)
betaSP_inf = rnorm(n = N_samples, mean = 0, sd = 0.1)
betaTG1_inf = rnorm(n = N_samples, mean = 0.5, 1)
betaTG2_inf = rnorm(n = N_samples, mean = 0.5, 1)
betaTG3_inf = rnorm(n = N_samples, mean = -0.5, 1)
betaTG4_inf = rnorm(n = N_samples, mean = -0.5, 1)
betaDPI1_inf = rnorm(n = N_samples, mean = -0.5, 1)
betaDPI2_inf = rnorm(n = N_samples, mean = 0.5, 1)
betaDPI3_inf = rnorm(n = N_samples, mean = 0, 1)

##### Non-Informative ------------------------------------------------------------

#alpha_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
#betaT_noninf = rnorm(n = N_samples, mean = 0, sd = 1); #hist(betaT)
alpha_noninf = rnorm(n = N_samples, mean = -2, sd = 1)
betaT_noninf = rgamma(n = N_samples, shape = 2, scale = 0.5); #hist(betaT)
betaST_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
betaDOSE_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
betaAGE_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
betaSEX_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
betaSP_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
betaTG1_noninf = rnorm(n = N_samples, mean = 0, 1)
betaTG2_noninf = rnorm(n = N_samples, mean = 0, 1)

betaTG3_noninf = rnorm(n = N_samples, mean = 0, 1)
betaDPI1_noninf = rnorm(n = N_samples, mean = 0, 1)
betaDPI2_noninf = rnorm(n = N_samples, mean = 0, 1)
betaDPI3_noninf = rnorm(n = N_samples, mean = 0, 1)


##### Generate predictions -------------------------------------------------------

pred_df_inf <- as.data.frame(t_sim)
pred_df_noninf <- as.data.frame(t_sim)

for (row_num in 1:nrow(prior_df)) {
  y_inf <-  alpha_inf[row_num] + betaT_inf[row_num] * prior_df$t_sim + 
            betaST_inf[row_num] * prior_df$st_sim[row_num] + 
            betaDOSE_inf[row_num] * prior_df$dose_sim[row_num] + 
            betaAGE_inf[row_num] * prior_df$age_sim[row_num] + 
            betaSEX_inf[row_num] * prior_df$sex_sim[row_num] + 
            betaSP_inf[row_num] * prior_df$sp_sim[row_num] +
            betaTG1_inf[row_num] * prior_df$tg_sim1[row_num] +
            betaTG2_inf[row_num] * prior_df$tg_sim2[row_num] +
            betaTG3_inf[row_num] * prior_df$tg_sim3[row_num] +
            betaDPI1_inf[row_num] * prior_df$dpi_sim1[row_num] +
            betaDPI2_inf[row_num] * prior_df$dpi_sim2[row_num] +
            betaDPI3_inf[row_num] * prior_df$dpi_sim3[row_num]
  y_noninf <-  alpha_noninf[row_num] + betaT_noninf[row_num] * prior_df$t_sim + 
               betaST_noninf[row_num] * prior_df$st_sim[row_num] + 
               betaDOSE_noninf[row_num] * prior_df$dose_sim[row_num] + 
               betaAGE_noninf[row_num] * prior_df$age_sim[row_num] + 
               betaSEX_noninf[row_num] * prior_df$sex_sim[row_num] + 
               betaSP_noninf[row_num] * prior_df$sp_sim[row_num] +
               betaTG1_noninf[row_num] * prior_df$tg_sim1[row_num] +
               betaTG2_noninf[row_num] * prior_df$tg_sim2[row_num] +
               betaTG3_noninf[row_num] * prior_df$tg_sim3[row_num] +
               betaDPI1_noninf[row_num] * prior_df$dpi_sim1[row_num] +
               betaDPI2_noninf[row_num] * prior_df$dpi_sim2[row_num] +
               betaDPI3_noninf[row_num] * prior_df$dpi_sim3[row_num]
  pred_df_inf <- cbind(pred_df_inf, y_inf)
  pred_df_noninf <- cbind(pred_df_noninf, y_noninf)
}

# Convert colnames and pivot longer to  maintain groups & plot
colnames(pred_df_inf) <- c("t_val", as.character(1:N_samples))
pred_df_inf <- pred_df_inf %>% pivot_longer(cols = !c("t_val"),
                                    names_to = c("sample"))
pred_df_inf$prior <- "Informative"

colnames(pred_df_noninf) <- c("t_val", as.character(1:N_samples))
pred_df_noninf <- pred_df_noninf %>% pivot_longer(cols = !c("t_val"),
                                            names_to = c("sample"))
pred_df_noninf$prior <- "Non-informative"


# Combine together
pred_df <- rbind(pred_df_inf, pred_df_noninf)


##### Make the plot  -------------------------------------------------------------

pred_df$y_facet <- "PCR Linear"
figSXB <- ggplot(data = pred_df) +
  geom_line(aes(x = t_val, y = value, group = sample, color = prior), 
            alpha = 0.8, lwd = 0.5) +
  labs(y = "log10 sgRNA copies / sample", x = element_blank()) +
  facet_grid(y_facet ~ prior) +
  scale_color_manual(values = c("Informative" = zissou_pal[1],
                                "Non-informative" = zissou_pal[12])) +
  scale_x_continuous(limits = c(0, 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-2, 10), expand = c(0, 0),
                     breaks = seq(0, 10, 2), labels = seq(0, 10, 2)) +
  theme(legend.position = "none", 
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title.y = element_text(size = 9),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.background = element_rect(fill='white'),
        text = element_text(size = 12),
        strip.text.x = element_blank(),
        strip.text = element_text(size = 9),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(size = 0.1, linetype = 'solid',
                                          colour = "light grey"),
        plot.margin = margin(2, 2, 2, 2)); figSXB


## Culture Model --------------------------------------------------------------------

### Logistic model ---------------------------------------------------------------

# Generate simulated data
N_samples <- 200
t_sim <- seq(from = 0, to = 10, length.out = N_samples)
st_sim <- sample(x = c(1, 2), size = N_samples, replace = TRUE)
dose_sim <- seq(from = 4, to = 7.5, length.out = N_samples)
age_sim <- sample(x = 1:3, size = N_samples, replace = TRUE)
sex_sim <- sample(x = 1:2, size = N_samples, replace = TRUE)
sp_sim <- sample(x = 1:3, size = N_samples, replace = TRUE)
tg_sim <- sample(x = 1:3, size = N_samples, replace = TRUE)
dpi_sim <- sample(x = 1:3, size = N_samples, replace = TRUE)
cell_sim <- sample(x = 1:3, size = N_samples, replace = TRUE)
assay_sim <- sample(x = 1:2, size = N_samples, replace = TRUE)

# Generate indicator cols for TG types
tg_sim1 <- tg_sim
tg_sim1[tg_sim == 1] <- 1
tg_sim1[tg_sim != 1] <- 0

tg_sim2 <- tg_sim
tg_sim2[tg_sim == 2] <- 1
tg_sim2[tg_sim != 2] <- 0

tg_sim3 <- tg_sim
tg_sim3[tg_sim == 3] <- 1
tg_sim3[tg_sim != 3] <- 0

# Generate indicator cols for DPI types
dpi_sim1 <- dpi_sim
dpi_sim1[dpi_sim == 1] <- 1
dpi_sim1[dpi_sim != 1] <- 0

dpi_sim2 <- dpi_sim
dpi_sim2[dpi_sim == 2] <- 1
dpi_sim2[dpi_sim != 2] <- 0

dpi_sim3 <- dpi_sim
dpi_sim3[dpi_sim == 3] <- 1
dpi_sim3[dpi_sim != 3] <- 0

# Generate indicator cols for cell types
cell_sim1 <- cell_sim
cell_sim1[cell_sim == 1] <- 1
cell_sim1[cell_sim != 1] <- 0

cell_sim2 <- cell_sim
cell_sim2[cell_sim == 2] <- 1
cell_sim2[cell_sim != 2] <- 0

cell_sim3 <- cell_sim
cell_sim3[cell_sim == 3] <- 1
cell_sim3[cell_sim != 3] <- 0


# Generate 
prior_df <- as.data.frame(cbind(t_sim, st_sim, dose_sim, age_sim, sex_sim, sp_sim,
                                cell_sim, cell_sim1, cell_sim2, cell_sim3,
                                assay_sim,
                                tg_sim, tg_sim1, tg_sim2, tg_sim3,
                                dpi_sim, dpi_sim1, dpi_sim2, dpi_sim3))



#### Simulate coefficients ------------------------------------------------------

##### Informative ---------------------------------------------------------------

gamma_inf = rnorm(n = N_samples, mean = -1, sd = 1)
psiT_inf = rnorm(n = N_samples, mean = 1, sd = 0.5)
psiST_inf = rnorm(n = N_samples, mean = 0, sd = 1)
psiDOSE_inf = rnorm(n = N_samples, mean = 0, sd = 1)
psiAGE_inf = rnorm(n = N_samples, mean = 0, sd = 1)
psiSEX_inf = rnorm(n = N_samples, mean = 0, sd = 1)
psiSP_inf = rnorm(n = N_samples, mean = 0, sd = 1)
psiTG1_inf = rnorm(n = N_samples, mean = -0.5, sd = 1)
psiTG2_inf = rnorm(n = N_samples, mean = 0, sd = 1)
psiTG3_inf = rnorm(n = N_samples, mean = 0.5, sd = 1)
psiDPI1_inf = rnorm(n = N_samples, mean = 0, sd = 1)
psiDPI2_inf = rnorm(n = N_samples, mean = 0, sd = 1)
psiDPI3_inf = rnorm(n = N_samples, mean = 0, sd = 1)
psiCELL1_inf = rnorm(n = N_samples, mean = 0, sd = 1)
psiCELL2_inf = rnorm(n = N_samples, mean = 0, sd = 1)
psiCELL3_inf = rnorm(n = N_samples, mean = 0.5, sd = 1)
psiASSAY_inf = rnorm(n = N_samples, mean = -0.5, sd = 1)


##### Non-Informative ------------------------------------------------------------

gamma_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
psiT_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
psiST_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
psiDOSE_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
psiAGE_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
psiSEX_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
psiSP_noninf = rnorm(n = N_samples, mean = 0, sd = 1)
psiTG1_noninf = rnorm(n = N_samples, mean = 0, 1)
psiTG2_noninf = rnorm(n = N_samples, mean = 0, 1)
psiTG3_noninf = rnorm(n = N_samples, mean = 0, 1)
psiDPI1_noninf = rnorm(n = N_samples, mean = 0, 1)
psiDPI2_noninf = rnorm(n = N_samples, mean = 0, 1)
psiDPI3_noninf = rnorm(n = N_samples, mean = 0, 1)
psiCELL1_noninf = rnorm(n = N_samples, mean = 0, 1)
psiCELL2_noninf = rnorm(n = N_samples, mean = 0, 1)
psiCELL3_noninf = rnorm(n = N_samples, mean = 0, 1)
psiASSAY_noninf = rnorm(n = N_samples, mean = 0, 1)


##### Generate predictions -------------------------------------------------------

pred_df_inf <- as.data.frame(t_sim)
pred_df_noninf <- as.data.frame(t_sim)

for (row_num in 1:nrow(prior_df)) {
  y_inf <- gamma_inf[row_num] + psiT_inf[row_num] * prior_df$t_sim + 
    psiST_inf[row_num] * prior_df$st_sim[row_num] + 
    psiDOSE_inf[row_num] * prior_df$dose_sim[row_num] + 
    psiAGE_inf[row_num] * prior_df$age_sim[row_num] + 
    psiSEX_inf[row_num] * prior_df$sex_sim[row_num] + 
    psiSP_inf[row_num] * prior_df$sp_sim[row_num] +
    psiTG1_inf[row_num] * prior_df$tg_sim1[row_num] +
    psiTG2_inf[row_num] * prior_df$tg_sim2[row_num] +
    psiTG3_inf[row_num] * prior_df$tg_sim3[row_num] +
    psiDPI1_inf[row_num] * prior_df$dpi_sim1[row_num] +
    psiDPI2_inf[row_num] * prior_df$dpi_sim2[row_num] +
    psiDPI3_inf[row_num] * prior_df$dpi_sim3[row_num] +
    psiCELL1_inf[row_num] * prior_df$cell_sim1[row_num] +
    psiCELL2_inf[row_num] * prior_df$cell_sim2[row_num] +
    psiCELL3_inf[row_num] * prior_df$cell_sim3[row_num] +
    psiASSAY_inf[row_num] * prior_df$assay_sim[row_num]
  y_noninf <- gamma_noninf[row_num] + psiT_noninf[row_num] * prior_df$t_sim + 
    psiST_noninf[row_num] * prior_df$st_sim[row_num] + 
    psiDOSE_noninf[row_num] * prior_df$dose_sim[row_num] + 
    psiAGE_noninf[row_num] * prior_df$age_sim[row_num] + 
    psiSEX_noninf[row_num] * prior_df$sex_sim[row_num] + 
    psiSP_noninf[row_num] * prior_df$sp_sim[row_num] +
    psiTG1_noninf[row_num] * prior_df$tg_sim1[row_num] +
    psiTG2_noninf[row_num] * prior_df$tg_sim2[row_num] +
    psiTG3_noninf[row_num] * prior_df$tg_sim3[row_num] +
    psiDPI1_noninf[row_num] * prior_df$dpi_sim1[row_num] +
    psiDPI2_noninf[row_num] * prior_df$dpi_sim2[row_num] +
    psiDPI3_noninf[row_num] * prior_df$dpi_sim3[row_num] +
    psiCELL1_noninf[row_num] * prior_df$cell_sim1[row_num] +
    psiCELL2_noninf[row_num] * prior_df$cell_sim2[row_num] +
    psiCELL3_noninf[row_num] * prior_df$cell_sim3[row_num] +
    psiASSAY_noninf[row_num] * prior_df$assay_sim[row_num]
  
  y_inf <- exp(y_inf)/(1 + exp(y_inf))
  y_noninf <- exp(y_noninf)/(1 + exp(y_noninf))
  pred_df_inf <- cbind(pred_df_inf, y_inf)
  pred_df_noninf <- cbind(pred_df_noninf, y_noninf)
}

# Convert colnames and pivot longer to maintain groups & plot
colnames(pred_df_inf) <- c("t_val", paste0("y", 1:N_samples))
pred_df_inf <- pred_df_inf %>% pivot_longer(cols = !c("t_val"),
                                            names_to = c("sample"))
pred_df_inf$prior <- "Informative"

colnames(pred_df_noninf) <- c("t_val", paste0("y", 1:N_samples))
pred_df_noninf <- pred_df_noninf %>% pivot_longer(cols = !c("t_val"),
                                                  names_to = c("sample"))
pred_df_noninf$prior <- "Non-informative"


# Combine into one dataframe for facetting
pred_df <- rbind(pred_df_inf, pred_df_noninf)


##### Make the plot  -------------------------------------------------------------

pred_df$y_facet <- "Culture Logistic"
figSXC <- ggplot(data = pred_df) +
  geom_line(aes(x = t_val, y = value, group = sample, color = prior), 
            alpha = 0.8, lwd = 0.5) +
  labs(y = "Chance of culture positive (%)", x = "log10 total RNA copies / sample") +
  #facet_wrap(.~ prior) +
  facet_grid(y_facet ~ prior) +
  scale_color_manual(values = c("Informative" = zissou_pal[1],
                                "Non-informative" = zissou_pal[12])) +
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, 2),
                     expand = c(0, 0), labels = c("", seq(2, 10, 2))) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme(legend.position = "none", 
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key = element_blank(),
        axis.title.y = element_text(size = 9),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.background = element_rect(fill='white'),
        text = element_text(size = 12),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        strip.text.x = element_blank(),
        strip.text = element_text(size = 9),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(size = 0.1, linetype = 'solid',
                                          colour = "light grey"),
        plot.margin = margin(2, 2, 2, 2)); figSXC



## Combine them into 1 plot  ----------------------------------------------------

figSXA <- figSXA + figSXB +  figSXC + plot_layout(nrow = 3); figSXA





# Parameter estimates ----------------------------------------------------------

## Load model fits -------------------------------------------------------------

fit.best.inf <- readRDS(file = "./outputs/fits/fit-sgRNA-best-inf-priors.RDS")
fit.best.noninf <- readRDS(file = "./outputs/fits/fit-sgRNA-best-noninf-priors.RDS")

fit.culture.inf <- readRDS(file = "./outputs/fits/fit-culture-best-inf-priors.RDS")
fit.culture.noninf <- readRDS(file = "./outputs/fits/fit-culture-best-noninf-priors.RDS")


## PCR Logistic ----------------------------------------------------------------

### Extract samples ------------------------------------------------------------

# Extract parameter estimates from best model with informative priors
pars_best <- c("gamma", "deltaT", 
               "deltaDOSE", 
               "deltaTG[1]", "deltaTG[2]", "deltaTG[3]", "deltaTG[4]",
               "deltaSP[1]", "deltaSP[2]", "deltaSP[3]")
param.best.df.wide <- fit.best.inf$draws(variables = pars_best,
                                     inc_warmup = FALSE,
                                     format = "df")
param.best.df <- param.best.df.wide %>% pivot_longer(cols = everything(),
                                                     names_to = "param",
                                                     values_to = "value")
param.best.df <- as.data.frame(param.best.df)
param.best.df <- subset(param.best.df, param %notin% c(".chain", ".iteration", ".draw"))


# Extract parameter estimates from best model with non-informative priors
param.best.noninf.df.wide <- fit.best.noninf$draws(variables = pars_best,
                                                   inc_warmup = FALSE,
                                                   format = "df")
param.best.noninf.df <- param.best.noninf.df.wide %>% pivot_longer(cols = everything(),
                                                     names_to = "param",
                                                     values_to = "value")
param.best.noninf.df <- as.data.frame(param.best.noninf.df)
param.best.noninf.df <- subset(param.best.noninf.df, param %notin% c(".chain", ".iteration", ".draw"))

# Combine
param.best.df$prior <- "Informative"
param.best.noninf.df$prior <- "Non-informative"
param.log <- rbind(param.best.df, param.best.noninf.df)

# Set more informative names 
param.log$param_short <- str_remove_all(param.log$param, "delta")
param.log$param_short[param.log$param_short == "TG[1]"] <- "TG [T↑ SG↑]"
param.log$param_short[param.log$param_short == "TG[2]"] <- "TG [T↓ SG↑]"
param.log$param_short[param.log$param_short == "TG[3]"] <- "TG [T↑ SG↓]"
param.log$param_short[param.log$param_short == "TG[4]"] <- "TG [T↓ SG↓]"
param.log$param_short[param.log$param_short == "SP[1]"] <- "SP [RM]"
param.log$param_short[param.log$param_short == "SP[2]"] <- "SP [CM]"
param.log$param_short[param.log$param_short == "SP[3]"] <- "SP [AGM]"

param.log$param_short_factor <- factor(param.log$param_short, 
                                       levels = rev(c("gamma", "T",
                                                      "TG [T↑ SG↑]",
                                                      "TG [T↓ SG↑]",
                                                      "TG [T↑ SG↓]",
                                                      "TG [T↓ SG↓]",
                                                      "SP [RM]",
                                                      "SP [CM]",
                                                      "SP [AGM]",
                                                      "DPI [I, 1]",
                                                      "DPI [I, 2+]",
                                                      "DPI [NI, 1+]",
                                                      "DOSE")))


### Plot -----------------------------------------------------------------------

figSXB1 <- ggplot(param.log, aes(x = value, y = param_short_factor, fill = prior)) + 
  geom_vline(xintercept = 0, size = 0.6, color = zissou_pal[10], alpha = 0.3) +
  geom_density_ridges(scale = 1, rel_min_height = 0.01, alpha = 0.7) + 
  labs(y = "", x = "Parameter value", tag = "B") +
  facet_grid(.~"PCR Logistic" ) +
  scale_fill_manual(values = c("Informative" = zissou_pal[1],
                               "Non-informative" = zissou_pal[12])) +
  scale_x_continuous(breaks = seq(from = -6, to = 6, by = 2), 
                     limits = c(-7.5, 7.5), expand = c(0, 0)) +
  theme(legend.position = "none",
        text = element_text(size = 12),
        #axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 8),
        #axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 9),
        legend.background = element_rect(fill='white'),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(0, 0, 0, 0)); figSXB1



## PCR Linear ------------------------------------------------------------------

### Extract samples ------------------------------------------------------------

# Extract parameter estimates from best model with informative priors
pars_best <- c("alpha", "betaT", 
               #"betaTG[1]", # excluded because no quantitative sgRNA vals for this category
               "betaTG[2]", "betaTG[3]", "betaTG[4]",
               "betaSP[1]", "betaSP[2]", "betaSP[3]",
               "betaDPI[1]", "betaDPI[2]", "betaDPI[3]",
               "betaDOSE")
param.best.df.wide <- fit.best.inf$draws(variables = pars_best,
                                     inc_warmup = FALSE,
                                     format = "df")
param.best.df.wide <- fit.best.inf$draws(variables = pars_best,
                                         inc_warmup = FALSE,
                                         format = "df")
param.best.df <- param.best.df.wide %>% pivot_longer(cols = everything(),
                                                     names_to = "param",
                                                     values_to = "value")
param.best.df <- as.data.frame(param.best.df)
param.best.df <- subset(param.best.df, param %notin% c(".chain", ".iteration", ".draw"))


# Extract parameter estimates from best model with non-informative priors
param.best.noninf.df.wide <- fit.best.noninf$draws(variables = pars_best,
                                                   inc_warmup = FALSE,
                                                   format = "df")
param.best.noninf.df <- param.best.noninf.df.wide %>% pivot_longer(cols = everything(),
                                                                   names_to = "param",
                                                                   values_to = "value")
param.best.noninf.df <- as.data.frame(param.best.noninf.df)
param.best.noninf.df <- subset(param.best.noninf.df, param %notin% c(".chain", ".iteration", ".draw"))


# Combine
param.best.df$prior <- "Informative"
param.best.noninf.df$prior <- "Non-informative"
param.lin <- rbind(param.best.df, param.best.noninf.df)

# Set more informative names 
param.lin$param_short <- str_remove_all(param.lin$param, "beta")
param.lin$param_short[param.lin$param_short == "DPI[1]"] <- "DPI [I, 1]"
param.lin$param_short[param.lin$param_short == "DPI[2]"] <- "DPI [I, 2+]"
param.lin$param_short[param.lin$param_short == "DPI[3]"] <- "DPI [NI, 1+]"
param.lin$param_short[param.lin$param_short == "TG[1]"] <- "TG [T↑ SG↑]"
param.lin$param_short[param.lin$param_short == "TG[2]"] <- "TG [T↓ SG↑]"
param.lin$param_short[param.lin$param_short == "TG[3]"] <- "TG [T↑ SG↓]"
param.lin$param_short[param.lin$param_short == "TG[4]"] <- "TG [T↓ SG↓]"
param.lin$param_short[param.lin$param_short == "SP[1]"] <- "SP [RM]"
param.lin$param_short[param.lin$param_short == "SP[2]"] <- "SP [CM]"
param.lin$param_short[param.lin$param_short == "SP[3]"] <- "SP [AGM]"

param.lin$param_short_factor <- factor(param.lin$param_short, 
                                       levels = rev(c("alpha", "T", 
                                                      "TG [T↓ SG↑]",
                                                      "TG [T↑ SG↓]",
                                                      "TG [T↓ SG↓]",
                                                      "SP [RM]",
                                                      "SP [CM]",
                                                      "SP [AGM]",
                                                      "DPI [I, 1]",
                                                      "DPI [I, 2+]",
                                                      "DPI [NI, 1+]",
                                                      "DOSE")))




### Plot -----------------------------------------------------------------------

figSXB2 <- ggplot(param.lin, aes(x = value, y = param_short_factor, fill = prior)) + 
  geom_vline(xintercept = 0, size = 0.6, color = zissou_pal[10], alpha = 0.3) +
  geom_density_ridges(scale = 1, rel_min_height = 0.01, alpha = 0.7) + 
  labs(y = "", x = "Parameter value") +
  facet_grid(.~ "PCR Linear") +
  scale_fill_manual(values = c("Informative" = zissou_pal[1],
                               "Non-informative" = zissou_pal[12])) +
  scale_x_continuous(breaks = seq(from = -4, to = 4, by = 2), 
                     limits = c(-4.5, 4.5), expand = c(0, 0)) +
  theme(legend.position = "none",
        #axis.text.x = element_blank(),
        #axis.title.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        text = element_text(size = 10),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        strip.text = element_text(size = 9),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(0, 0, 0, 0)); figSXB2



## Culture Logistic ------------------------------------------------------------

### Extract samples ------------------------------------------------------------

# Extract parameter estimates from best model with informative priors
pars_culture <- c("gamma", "psiT", 
               "psiDPI[1]",  "psiDPI[2]", "psiDPI[3]",
               "psiAGE[1]", "psiAGE[2]", "psiAGE[3]",
               "psiSP[1]", "psiSP[2]", "psiSP[3]",
               "psiASSAY",
               "psiTG[1]", "psiTG[2]", "psiTG[3]", 
               "psiDOSE", 
               "psiCELL[1]", "psiCELL[2]", "psiCELL[3]")

param.culture.df.wide <- fit.culture.inf$draws(variables = pars_culture,
                                            inc_warmup = FALSE,
                                            format = "df")
param.culture.df <- param.culture.df.wide %>% pivot_longer(cols = everything(),
                                                     names_to = "param",
                                                     values_to = "value")
param.culture.df <- as.data.frame(param.culture.df)
param.culture.df <- subset(param.culture.df, param %notin% c(".chain", ".iteration", ".draw"))


# Extract parameter estimates from culture model with non-informative priors
param.culture.noninf.df.wide <- fit.culture.noninf$draws(variables = pars_culture,
                                                         inc_warmup = FALSE,
                                                         format = "df")
param.culture.noninf.df <- param.culture.noninf.df.wide %>% pivot_longer(cols = everything(),
                                                                         names_to = "param",
                                                                         values_to = "value")
param.culture.noninf.df <- as.data.frame(param.culture.noninf.df)
param.culture.noninf.df <- subset(param.culture.noninf.df, param %notin% c(".chain", ".iteration", ".draw"))

# Combine
param.culture.df$prior <- "Informative"
param.culture.noninf.df$prior <- "Non-informative"
param.cul <- rbind(param.culture.df, param.culture.noninf.df)

# Set more informative names 
param.cul$param_short <- str_remove_all(param.cul$param, "psi")
param.cul$param_short[param.cul$param_short == "DPI[1]"] <- "DPI [I, 1]"
param.cul$param_short[param.cul$param_short == "DPI[2]"] <- "DPI [I, 2+]"
param.cul$param_short[param.cul$param_short == "DPI[3]"] <- "DPI [NI, 1+]"
param.cul$param_short[param.cul$param_short == "TG[1]"] <- "TG [N]"
param.cul$param_short[param.cul$param_short == "TG[2]"] <- "TG [E]"
param.cul$param_short[param.cul$param_short == "TG[3]"] <- "TG [S]"
param.cul$param_short[param.cul$param_short == "SP[1]"] <- "SP [RM]"
param.cul$param_short[param.cul$param_short == "SP[2]"] <- "SP [CM]"
param.cul$param_short[param.cul$param_short == "SP[3]"] <- "SP [AGM]"
param.cul$param_short[param.cul$param_short == "AGE[1]"] <- "AGE [Juvenile]"
param.cul$param_short[param.cul$param_short == "AGE[2]"] <- "AGE [Adult]"
param.cul$param_short[param.cul$param_short == "AGE[3]"] <- "AGE [Geriatric]"
param.cul$param_short[param.cul$param_short == "CELL[1]"] <- "CELL [76]"
param.cul$param_short[param.cul$param_short == "CELL[2]"] <- "CELL [E6]"
param.cul$param_short[param.cul$param_short == "CELL[3]"] <- "CELL [E6-SS2]"

param.cul$param_short_factor <- factor(param.cul$param_short, 
                                       levels = rev(c("gamma", "T", 
                                                      "TG [N]",
                                                      "TG [E]",
                                                      "TG [S]",
                                                      "SP [RM]",
                                                      "SP [CM]",
                                                      "SP [AGM]",
                                                      "AGE [Juvenile]",
                                                      "AGE [Adult]",
                                                      "AGE [Geriatric]",
                                                      "DPI [I, 1]",
                                                      "DPI [I, 2+]",
                                                      "DPI [NI, 1+]",
                                                      "DOSE",
                                                      "CELL [76]",
                                                      "CELL [E6]",
                                                      "CELL [E6-SS2]",
                                                      "ASSAY")))




### Plot -----------------------------------------------------------------------

figSXB3 <- ggplot(param.cul, aes(x = value, y = param_short_factor, fill = prior)) + 
  geom_vline(xintercept = 0, size = 0.6, color = zissou_pal[10], alpha = 0.3) +
  geom_density_ridges(scale = 1, rel_min_height = 0.01, alpha = 0.7) + 
  labs(y = "", x = "Parameter value") +
  facet_grid(.~ "Culture Logistic" ) +
  scale_fill_manual(values = c("Informative" = zissou_pal[1],
                               "Non-informative" = zissou_pal[12])) +
  scale_x_continuous(breaks = seq(from = -4, to = 4, by = 2), 
                     limits = c(-4.5, 4.5), expand = c(0, 0)) +
  theme(legend.position = "none",
        text = element_text(size = 12),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 9),
        legend.background = element_rect(fill='white'),
        strip.background  = element_rect(colour = "black", fill = "grey92"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "light grey"),
        plot.margin = margin(0, 0, 0, 0)); figSXB3



## Combine them into 1 plot  ----------------------------------------------------

figSXB <- figSXB1 + figSXB2 + figSXB3 + 
  plot_layout(ncol = 3, widths = c(0.8, 0.8, 0.8)); figSXB


# Combine into FigSX  ----------------------------------------------------------

figSX <- ggarrange(figSXA, figSXB, ncol = 2, widths = c(2, 3)); figSX

ggsave('./outputs/figures/figS13-prior-checks-and-params.tiff',
       plot = figSX,
       device = 'png',
       height = 7,
       width = 12,
       units = 'in',
       bg = 'white')



