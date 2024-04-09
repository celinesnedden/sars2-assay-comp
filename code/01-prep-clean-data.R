# Prepares clean data from the raw data
# Run this file first, if verifying the generation of the clean data
# Otherwise, skip to subsequent files to use the cleaned data directly

# Prep environment -------------------------------------------------------------

# Install required packages 
req_pkgs <- c("stringr", "tidyverse")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependenices = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))

# Function for convenience
`%notin%` <- Negate(`%in%`) 


# Load raw data ----------------------------------------------------------------

dat <- read.csv("./data/raw-data.csv")


# Update meta-data  ------------------------------------------------------------

## Consolidate target genes ----------------------------------------------------

dat$tg_total[dat$tg_total == "N (N2)"] <- "N"
dat$tg_total[dat$tg_total == "S (RBD)"] <- "S"

 
## Convert simple categorical predictors  --------------------------------------

# Sample types
dat$st_idx[dat$sample_type == "Non-invasive"] <- 0
dat$st_idx[dat$sample_type == "Invasive"] <- 1

# Species
dat$sp_idx[dat$animal_species == "rhesus macaque"] <- 1
dat$sp_idx[dat$animal_species == "cynomolgus macaque"] <- 2
dat$sp_idx[dat$animal_species == "african green monkey"] <- 3

# Target gene classification based on high/low totRNA and sgRNA quantities
dat$tg_idx[dat$tg_total == "N" & dat$tg_sg == "N"] <- 1 # T-High, SG-High
dat$tg_idx[dat$tg_total == "E" & dat$tg_sg == "N"] <- 2 # T-Low, SG-High
dat$tg_idx[dat$tg_total == "N" & dat$tg_sg %in% c("E", "ORF7")] <- 3 # T-High, SG-Low
dat$tg_idx[dat$tg_total == "E" & dat$tg_sg %in% c("E", "ORF7")] <- 4 # T-Low, SG-Low
# Note: samples without sgRNA values will return NA 

# Sex
dat$sex_idx[dat$sex == "Male"] <- 0
dat$sex_idx[dat$sex == "Female"] <- 1
dat$sex_idx[dat$sex == "Unknown"] <- -9

# Age
dat$age_idx[dat$age_class_stand == "Juvenile"] <- 1
dat$age_idx[dat$age_class_stand == "Adult"] <- 2
dat$age_idx[dat$age_class_stand == "Geriatric"] <- 3
dat$age_idx[dat$age_class_stand == "Unknown"] <- -9

# Cell line
# Note: samples without culture data will be NA 
dat$cell_idx[dat$cell_line == "Vero E6"] <- 1
dat$cell_idx[dat$cell_line == "Vero 76"] <- 2
dat$cell_idx[dat$cell_line == "Vero E6/TMPRSS2"] <- 3

# Culture assay
# Note: samples without culture data will be NA 
dat$assay_idx[dat$culture_assay == "plaque"] <- 1
dat$assay_idx[dat$culture_assay == "TCID50"] <- 0



## Classify organ systems ------------------------------------------------------

dat$organ_system[dat$sample_grp %in% c("Nasal Sample", 
                                       "Nasopharyngeal Sample", 
                                       "Oropharynx",
                                       "Oropharyngeal Sample",
                                       "Nasal tissue (Cavity)", 
                                       "Nasal tissue (Mucosa)")] <- "URT"

dat$organ_system[dat$sample_grp %in% c("BAL", "Lung (Right)", "Lung (Left)" , 
                                       "Lung (Right upper)",
                                       "Lung (Right middle)", 
                                       "Lung (Right lower)", 
                                       "Lung (Left upper)", 
                                       "Lung (Left middle)",
                                       "Lung (Left lower)", 
                                       "Lung (Lower)", 
                                       "Trachea", 
                                       "Lung (Upper)",
                                       "Lung",
                                       "Lung (Azygos / accessory)",
                                       "Lung (Average of 7 lobes)", 
                                       "Bronchus (Right)", 
                                       "Bronchus (Left)",
                                       "Bronchial Sample", 
                                       "Tracheal Sample")] <- "LRT"

dat$organ_system[dat$sample_grp %in% c("Anal / Rectal Sample", 
                                       "Colon (Cecum)",
                                       "Colon", 
                                       "Colon (Descending)",
                                       "Rectum" ,
                                       "Small intestine", 
                                       "Stomach",
                                       "Small intestine (Duodenum)",
                                       "Small intestine (Jejunum)",
                                       "Small intestine (Ileum)")] <- "GI"

dat$organ_system[dat$organ_system %notin% c("URT", "LRT", "GI")] <- "Other"


## Generate DPI predictor  ------------------------------------------------------

# Based on day post infection, but stratified by whether the tissue was inoculated

in_inoc <- c("Nasal Sample")
it_inoc <- c()
oc_inoc <- c()
ig_inoc <- c("Stomach")
ae_inoc <- c("Oropharyngeal Sample", 
             "Nasopharyngeal Sample")
itin_inoc <- c("Nasal Sample", 
               "Nasopharyngeal Sample", 
               "Oropharyngeal Sample", 
               "Trachea", 
               "Nasal tissue (Cavity)", 
               "BAL")
itinoc_inoc <- c("Ocular Sample", 
                 "Nasal Sample",  
                 "Trachea",
                 "Oropharyngeal Sample")
multi_inoc <- c("Tonsil", 
                "Eye (Conjunctiva)", 
                "Nasal Sample", 
                "Oropharyngeal Sample", 
                "BAL", 
                "Nasal tissue (Mucosa)", 
                "Oropharynx", 
                "Trachea",
                "Ocular Sample", 
                "Oral Sample", 
                "Tracheal Sample", 
                "Bronchial Sample",
                "Bronchus (Right)", 
                "Bronchus (Left)")

dat <- dat %>% 
  mutate(dpi_idx = 
           case_when(inoc_route == "IN" & 
                       sample_grp %in% in_inoc &
                       dpi == 1 ~ 1,
                     
                     inoc_route == "IN" & 
                       sample_grp %in% in_inoc &
                       dpi != 1 ~ 2,
                     
                     inoc_route == "IG" & 
                       sample_grp %in% ig_inoc &
                       dpi == 1 ~ 1,
                     
                     inoc_route == "IG" & 
                       sample_grp %in% ig_inoc &
                       dpi != 1 ~ 2,
                     
                     inoc_route == "AE" & 
                       sample_grp %in% ae_inoc &
                       dpi == 1 ~ 1,
                     
                     inoc_route == "AE" & 
                       sample_grp %in% ae_inoc &
                       dpi != 1 ~ 2,
                      
                     inoc_route == "IT, IN" & 
                       sample_grp %in% itin_inoc &
                       dpi == 1 ~ 1,
                     
                     inoc_route == "IT, IN" & 
                       sample_grp %in% itin_inoc &
                       dpi != 1 ~ 2,
                     
                     inoc_route == "IT, IN, OC" & 
                       sample_grp %in% itinoc_inoc &
                       dpi == 1 ~ 1,
                     
                     inoc_route == "IT, IN, OC" & 
                       sample_grp %in% itinoc_inoc &
                       dpi != 1 ~ 2,
                     
                     inoc_route == "IT, IN, OR, OC" &
                       sample_grp %in% multi_inoc &
                       dpi == 1 ~ 1,
                     
                     inoc_route == "IT, IN, OR, OC" & 
                       sample_grp %in% multi_inoc &
                       dpi != 1 ~ 2)) %>%
  mutate(dpi_idx = replace(dpi_idx, is.na(dpi_idx), 3))




## Add log10_dose column -------------------------------------------------------

dat$log10_dose_pfu <- log10(as.numeric(dat$inoc_dose_total_pfu))



# Handle censored data ---------------------------------------------------------

dat$cens_total[dat$val_total %in% c("< LOD", "< LOQ")] <- 1
dat$cens_total[!is.na(dat$val_total) & dat$val_total %in% c("Pos (LOD < x < LOQ)")] <- 3 # positive but not quantified
dat$cens_total[!is.na(as.numeric(dat$val_total))] <- 0

dat$cens_sg[!is.na(dat$val_sg) & dat$val_sg %in% c("< LOD", "< LOQ", "Neg")] <- 1
dat$cens_sg[!is.na(dat$val_sg) & dat$val_sg %in% c("Pos", "Pos (LOD < x < LOQ)")] <- 3 # positive but not quantified
dat$cens_sg[!is.na(as.numeric(dat$val_sg))] <- 0

dat$pos_total[dat$val_total %in% c("< LOD", "< LOQ")] <- 0
dat$pos_total[!is.na(dat$val_total) & dat$val_total %in% c("Pos (LOD < x < LOQ)")] <- 1
dat$pos_total[!is.na(as.numeric(dat$val_total))] <- 1

dat$pos_sg[!is.na(dat$val_sg) & dat$val_sg %in% c("< LOD", "< LOQ", "Neg")] <- 0
dat$pos_sg[!is.na(dat$val_sg) & dat$val_sg %in% c("Pos", "Pos (LOD < x < LOQ)")] <- 1
dat$pos_sg[!is.na(as.numeric(dat$val_sg))] <- 1


## Remove total RNA positive censored values -----------------------------------

# We rely on positive totRNA samples having quantitative values
#   So we remove positive totRNA samples w/o quantitative information
dat <- subset(dat, cens_total %in% c(0, 1, NA))


## Add identifier for culture results ------------------------------------------

dat$pos_inf[!is.na(dat$val_inf) & dat$val_inf == "Pos"] <- 1
dat$pos_inf[!is.na(as.numeric(dat$val_inf))] <- 1
dat$pos_inf[!is.na(dat$val_inf) & dat$val_inf == "Neg"] <- 0


# Combine LOD/LOQ into 1 column ------------------------------------------------

dat$lodq_sg <- NA
for (row_num in 1:nrow(dat)) {
  if (!is.na(dat$pos_sg[row_num])) {
    if (dat$lod_sg[row_num] == "Unknown" & dat$loq_sg[row_num] == "Unknown") {
      dat$lodq_sg[row_num] <- -9 # Set -9 for Unknown
    }
    else if (dat$lod_sg[row_num] != "Unknown" & dat$loq_sg[row_num] == "Unknown") {
      dat$lodq_sg[row_num] <- as.numeric(dat$lod_sg[row_num])
    }
    else if (dat$lod_sg[row_num] == "Unknown" & dat$loq_sg[row_num] != "Unknown") {
      dat$lodq_sg[row_num] <- as.numeric(dat$loq_sg[row_num])
    }
    else if (dat$lod_sg[row_num] != "Unknown" & dat$loq_sg[row_num] != "Unknown") {
      dat$lodq_sg[row_num] <- as.numeric(dat$lod_sg[row_num])
    }
  }
}

dat$lodq_total <- NA
for (row_num in 1:nrow(dat)) {
  if (!is.na(dat$pos_total[row_num])) {
    if (dat$lod_total[row_num] == "Unknown" & dat$loq_total[row_num] == "Unknown") {
      dat$lodq_total[row_num] <- -9 # Set -9 for Unknown
    }
    else if (dat$lod_total[row_num] != "Unknown" & dat$loq_total[row_num] == "Unknown") {
      dat$lodq_total[row_num] <- as.numeric(dat$lod_total[row_num])
    }
    else if (dat$lod_total[row_num] == "Unknown" & dat$loq_total[row_num] != "Unknown") {
      dat$lodq_total[row_num] <- as.numeric(dat$loq_total[row_num])
    }
    else if (dat$lod_total[row_num] != "Unknown" & dat$loq_total[row_num] != "Unknown") {
      dat$lodq_total[row_num] <- as.numeric(dat$lod_total[row_num])
    }
  }
}


# Prepare for Stan -------------------------------------------------------------

# Stan does not allow character or NA entries, so these must be changed to
#     values, which won't be called during fitting nor will they affect inference

dat$val_sg[!is.na(dat$val_sg) & dat$cens_sg == 1] <- -9
dat$val_total[!is.na(dat$val_total) & dat$cens_total == 1] <- -9
dat$val_inf[!is.na(dat$val_inf) & dat$pos_inf == 0] <- -9


# Reorder columns for better visibility ----------------------------------------

dat <- dat %>%
  select(article, indiv, dpi, dpi_idx,
         sample_type, st_idx, sample_rep,
         val_total, pos_total, cens_total, 
         val_sg, pos_sg, cens_sg, 
         val_inf, pos_inf,
         lodq_sg, lodq_total,
         tg_total, tg_sg, tg_idx, 
         cell_line, cell_idx, 
         culture_assay, assay_idx,
         animal_species, sp_idx,
         age_class_stand, age_years_rep, age_class_rep, age_idx,
         sex, sex_idx, 
         units_total, units_sg, units_inf,
         lod_total, loq_total, lod_sg, loq_sg,
         everything()) 


# Export -----------------------------------------------------------------------

write.csv(dat, file = "./data/clean-data.csv")


