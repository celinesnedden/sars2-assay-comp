# Creates table S7 with the credible intervals for all parameter values from
#     our best model

# Load model fits -------------------------------------------------------------

fit.best.inf <- readRDS(file = "./outputs/fits/fit-sgRNA-best-inf-priors.RDS")
fit.culture.inf <- readRDS(file = "./outputs/fits/fit-culture-best-inf-priors.RDS")


# Prep dataframe ---------------------------------------------------------------


## PCR Logistic ----------------------------------------------------------------

pars_best <- c("gamma", "deltaT", 
               "deltaDOSE", 
               "deltaSP[1]", "deltaSP[2]", "deltaSP[3]",
               "deltaTG[1]", "deltaTG[2]", "deltaTG[3]", "deltaTG[4]"
               )


tbl.log <- as.data.frame(fit.best.inf$summary(pars_best))
tbl.log$model <- "PCR Logistic"

## PCR Linear ------------------------------------------------------------------

pars_best <- c("alpha", "betaT", 
               "betaDOSE",
               "betaDPI[1]", "betaDPI[2]", "betaDPI[3]",
               "betaSP[1]", "betaSP[2]", "betaSP[3]",
               "betaTG[2]", "betaTG[3]", "betaTG[4]"
               )

tbl.lin <- as.data.frame(fit.best.inf$summary(pars_best))
tbl.lin$model <- "PCR Linear"


## Culture Logistic ------------------------------------------------------------

pars_culture <- c("gamma", "psiT", 
                  "psiDOSE",
                  "psiDPI[1]",  "psiDPI[2]", "psiDPI[3]",
                  "psiSP[1]", "psiSP[2]", "psiSP[3]",
                  "psiAGE[1]", "psiAGE[2]", "psiAGE[3]",
                  "psiTG[1]", "psiTG[2]", "psiTG[3]", 
                  "psiASSAY",
                  "psiCELL[1]", "psiCELL[2]", "psiCELL[3]")

tbl.cul <- as.data.frame(fit.culture.inf$summary(pars_culture))
tbl.cul$model <- "Culture"


## Combine & order ------------------------------------------------------------

tbl.full <- rbind(tbl.log, tbl.lin, tbl.cul)
tbl.full <- subset(tbl.full, select = -c(ess_bulk, ess_tail, rhat, mean, mad))

tbl.full <- subset(tbl.full, select = c(model, variable,
                                        q5, median, q95, sd))

tbl.full$variable[str_detect(tbl.full$variable, "aTG\\[1")] <- "TG [T↑ SG↑]"
tbl.full$variable[str_detect(tbl.full$variable, "aTG\\[2")] <- "TG [T↓ SG↑]"
tbl.full$variable[str_detect(tbl.full$variable, "aTG\\[3")] <- "TG [T↑ SG↓]"
tbl.full$variable[str_detect(tbl.full$variable, "aTG\\[4")] <- "TG [T↓ SG↓]"
tbl.full$variable[str_detect(tbl.full$variable, "psiTG\\[1")] <- "TG [N]"
tbl.full$variable[str_detect(tbl.full$variable, "psiTG\\[2")] <- "TG [E]"
tbl.full$variable[str_detect(tbl.full$variable, "psiTG\\[3")] <- "TG [S]"

tbl.full$variable[str_detect(tbl.full$variable, "SP\\[1")] <- "SP [RM]"
tbl.full$variable[str_detect(tbl.full$variable, "SP\\[2")] <- "SP [CM]"
tbl.full$variable[str_detect(tbl.full$variable, "SP\\[3")] <- "SP [AGM]"

tbl.full$variable[str_detect(tbl.full$variable, "DPI\\[1")] <- "DPI [I, 1]"
tbl.full$variable[str_detect(tbl.full$variable, "DPI\\[2")] <- "DPI [I, 2+]"
tbl.full$variable[str_detect(tbl.full$variable, "DPI\\[3")] <- "DPI [NI, 1+]"

tbl.full$variable[str_detect(tbl.full$variable, "AGE\\[1")] <- "AGE [Juvenile]"
tbl.full$variable[str_detect(tbl.full$variable, "AGE\\[2")] <- "AGE [Adult]"
tbl.full$variable[str_detect(tbl.full$variable, "AGE\\[3")] <- "AGE [Geriatric]"

tbl.full$variable[str_detect(tbl.full$variable, "CELL\\[1")] <- "CELL [76]"
tbl.full$variable[str_detect(tbl.full$variable, "CELL\\[2")] <- "CELL [E6]"
tbl.full$variable[str_detect(tbl.full$variable, "CELL\\[3")] <- "CELL [E6-SS2]"

tbl.full$variable <- str_remove_all(tbl.full$variable, "beta|delta|psi")

tbl.full$variable[tbl.full$variable %in% c("gamma", "alpha")] <- "intercept"

tbl.full$q5 <- round(tbl.full$q5, 2)
tbl.full$q95 <- round(tbl.full$q95, 2)
tbl.full$median <- round(tbl.full$median, 2)
tbl.full$sd <- round(tbl.full$sd, 2)


# Save the table ----------------------------------------------------------

write.csv(tbl.full,
          file = './outputs/tables/tblS7-parameter-values.csv',
          row.names = FALSE)



