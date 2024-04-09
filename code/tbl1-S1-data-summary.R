# Generates Table 1 (simplified data summary) and table S1 (full data summary)
#   Note that the manuscript versions includes some minor additional post-processing
#   formatting for aesthetic purposes only


# Prep environment -------------------------------------------------------------

## Install necessary packages  -------------------------------------------------

req_pkgs <- c("stringr", "stringi", "dplyr")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependences = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))

`%notin%` <- Negate(`%in%`)  # for convenience


# Table 1: Combined version ----------------------------------------------------

## Load data -------------------------------------------------------------------

dat <- read.csv("./data/clean-data.csv")

## Change some index vars to be more informative for the table 
dat$tg_total[is.na(dat$tg_total) & !is.na(dat$tg_sg)] <- dat$tg_sg[is.na(dat$tg_total) & !is.na(dat$tg_sg)]

dat$dpi_idx[dat$dpi_idx == 1] <- "1: Inoc, 1 dpi"
dat$dpi_idx[dat$dpi_idx == 2] <- "2: Inoc, 2+ dpi"
dat$dpi_idx[dat$dpi_idx == 3] <- "3: Non-Inoc, 1+ dpi"

dat$animal_species[dat$sp_idx == 1] <- "Rhesus macaque"
dat$animal_species[dat$sp_idx == 2] <- "Cynomolgus macaque"
dat$animal_species[dat$sp_idx == 3] <- "African green monkey"

dat$culture_assay[dat$culture_assay == "plaque"] <- "Plaque"

## Group inoculation routes 
dat$inoc_group <- NA
dat$inoc_group[!str_detect(dat$inoc_route, ",")] <- "Single"
dat$inoc_group[str_detect(dat$inoc_route, ",")] <- "Multi"


## Simple stats for text -------------------------------------------------------

# These are all statistics reported in the first results paragraph

cat("The number of included articles is: ", length(unique(dat$article)), "\n")
cat("The number of available samples is: ", nrow(dat), "\n")
cat("The number of individuals sampled is: ", length(unique(dat$indiv)))

data.bothRNA <- subset(dat, !is.na(pos_total) & !is.na(pos_sg))
cat("The number of articles with totRNA and sgRNA is: ", length(unique(data.bothRNA$article)))
cat("The number of individuals with totRNA and sgRNA is: ", length(unique(data.bothRNA$indiv)))
cat("The number of samples with totRNA and sgRNA is: ", nrow(data.bothRNA))

data.culture <- subset(dat, !is.na(pos_inf))
cat("The number of articles with culture is: ", length(unique(data.culture$article)))
cat("The number of individuals with culture is: ", length(unique(data.culture$indiv)))
cat("The number of samples with culture is: ", nrow(data.culture))

data.all.three <- subset(dat, !is.na(pos_total) & !is.na(pos_sg) & !is.na(pos_inf))
cat("The number of articles with all three sample types is: ", length(unique(data.all.three$article)))
cat("The number of individuals with all three sample types is: ", length(unique(data.all.three$indiv)))
cat("The number of samples with all three sample types is: ", nrow(data.all.three))




## Prep table ------------------------------------------------------------------

cols <- c(" ", "sgRNA", "Culture", "Total")
rows <- c("Species", "Rhesus macaque", "Cynomolgus macaque", "African green monkey",
          "Age class", "Juvenile", "Adult", "Geriatric", "Unknown",
          "Sex", "Female", "Male", "Unknown",
          "Exposure dose", "10^4-<10^6", ">=10^6",
          "Exposure route", "Single", "Multi",
          "Sample type", "Invasive", "Non-invasive",
          "Sample time", "1: Inoc, 1 dpi", "2: Inoc, 2+ dpi", "3: Non-Inoc, 1+ dpi",
          "PCR target gene", "N", "E", "S",
          "Culture assay", "TCID50", "Plaque",
          "Cell line", "Vero E6", "Vero E6/TMPRSS2", "Vero 76",
          "Total")

tbl1 <- as.data.frame(matrix(nrow = length(rows), ncol = length(cols)))
colnames(tbl1) <- cols
tbl1$" " <- rows


## Fill in the table -----------------------------------------------------------

for (col_num in 2:(dim(tbl1)[2]-1)) {
  
  # Assay-specific subsets
  if (colnames(tbl1)[col_num] == "sgRNA") {
    dat.sub <- subset(dat, !is.na(pos_sg) & !is.na(pos_total))
  }
  else if (colnames(tbl1)[col_num] == "Culture") {
    dat.sub <- subset(dat, !is.na(pos_inf))
  }
  tbl1[nrow(tbl1), col_num] <- paste(dim(dat.sub)[1],
                                     length(unique(dat.sub$indiv)), 
                                     length(unique(dat.sub$article)),
                                     sep = " / ")
  
  # Species 
  for (sp.ii in unique(dat$animal_species)) {
    dat.sub.sp <- subset(dat.sub, animal_species == sp.ii)
    tbl1[which(tbl1$" " == sp.ii)[1], col_num] <- paste(dim(dat.sub.sp)[1],
                                                        length(unique(dat.sub.sp$indiv)), 
                                                        length(unique(dat.sub.sp$article)),
                                                        sep = " / ")
    dat.sp <- subset(dat, animal_species == sp.ii)
    tbl1[which(tbl1$" " == sp.ii)[1], ncol(tbl1)] <- paste(dim(dat.sp)[1],
                                                           length(unique(dat.sp$indiv)), 
                                                           length(unique(dat.sp$article)),
                                                           sep = " / ")
  }
  
  
  # Sex 
  for (sex.ii in unique(dat$sex)) {
    dat.sub.sex <- subset(dat.sub, sex == sex.ii)
    tbl1[which(tbl1$" " == sex.ii)[1], col_num] <- paste(dim(dat.sub.sex)[1],
                                                         length(unique(dat.sub.sex$indiv)), 
                                                         length(unique(dat.sub.sex$article)),
                                                         sep = " / ")
    dat.sex <- subset(dat, sex == sex.ii)
    tbl1[which(tbl1$" " == sex.ii)[1], ncol(tbl1)] <- paste(dim(dat.sex)[1],
                                                            length(unique(dat.sex$indiv)), 
                                                            length(unique(dat.sex$article)),
                                                            sep = " / ")
  }
  
  # Age class 
  for (age.ii in unique(dat$age_class_stand)) {
    dat.sub.age <- subset(dat.sub, age_class_stand == age.ii)
    if (age.ii == "Unknown") {
      tbl1[which(tbl1$" " == age.ii)[2], col_num] <- paste(dim(dat.sub.age)[1],
                                                           length(unique(dat.sub.age$indiv)), 
                                                           length(unique(dat.sub.age$article)),
                                                           sep = " / ")
      dat.age <- subset(dat, age_class_stand == age.ii)
      tbl1[which(tbl1$" " == age.ii)[2], ncol(tbl1)] <- paste(dim(dat.age)[1],
                                                              length(unique(dat.age$indiv)), 
                                                              length(unique(dat.age$article)),
                                                              sep = " / ")
    }
    else {
      tbl1[which(tbl1$" " == age.ii)[1], col_num] <- paste(dim(dat.sub.age)[1],
                                                           length(unique(dat.sub.age$indiv)), 
                                                           length(unique(dat.sub.age$article)),
                                                           sep = " / ")
      dat.age <- subset(dat, age_class_stand == age.ii)
      tbl1[which(tbl1$" " == age.ii)[1], ncol(tbl1)] <- paste(dim(dat.age)[1],
                                                              length(unique(dat.age$indiv)), 
                                                              length(unique(dat.age$article)),
                                                              sep = " / ")
    }
    
  }
  
  # Exposure dose
  dat.sub.mid <- subset(dat.sub, log10(as.numeric(inoc_dose_total_pfu)) < 6)
  dat.sub.hi <- subset(dat.sub, log10(as.numeric(inoc_dose_total_pfu)) >= 6)
  tbl1[which(tbl1$" " == "10^4-<10^6"), col_num] <- paste(dim(dat.sub.mid)[1],
                                                          length(unique(dat.sub.mid$indiv)), 
                                                          length(unique(dat.sub.mid$article)),
                                                          sep = " / ")
  tbl1[which(tbl1$" " == ">=10^6"), col_num] <- paste(dim(dat.sub.hi)[1],
                                                      length(unique(dat.sub.hi$indiv)), 
                                                      length(unique(dat.sub.hi$article)),
                                                      sep = " / ")
  dat.mid <- subset(dat, log10(as.numeric(inoc_dose_total_pfu)) < 6)
  dat.hi <- subset(dat, log10(as.numeric(inoc_dose_total_pfu)) >= 6)
  tbl1[which(tbl1$" " == "10^4-<10^6"), ncol(tbl1)] <- paste(dim(dat.mid)[1],
                                                             length(unique(dat.mid$indiv)), 
                                                             length(unique(dat.mid$article)),
                                                             sep = " / ")
  tbl1[which(tbl1$" " == ">=10^6"), ncol(tbl1)] <- paste(dim(dat.hi)[1],
                                                         length(unique(dat.hi$indiv)), 
                                                         length(unique(dat.hi$article)),
                                                         sep = " / ")
  
  
  # Exposure route
  for (route.ii in unique(dat$inoc_group)) {
    dat.sub.rt <- subset(dat.sub, inoc_group == route.ii)
    tbl1[which(tbl1$" " == route.ii), col_num] <- paste(dim(dat.sub.rt)[1],
                                                        length(unique(dat.sub.rt$indiv)), 
                                                        length(unique(dat.sub.rt$article)),
                                                        sep = " / ")
    
    dat.rt <- subset(dat, inoc_group == route.ii)
    tbl1[which(tbl1$" " == route.ii), ncol(tbl1)] <- paste(dim(dat.rt)[1],
                                                           length(unique(dat.rt$indiv)), 
                                                           length(unique(dat.rt$article)),
                                                           sep = " / ")
    
  }
  
  # Sample type 
  for (st.ii in unique(dat$sample_type)) {
    dat.sub.st <- subset(dat.sub, sample_type == st.ii)
    tbl1[which(tbl1$" " == st.ii), col_num] <- paste(dim(dat.sub.st)[1],
                                                     length(unique(dat.sub.st$indiv)), 
                                                     length(unique(dat.sub.st$article)),
                                                     sep = " / ")
    
    dat.st <- subset(dat, sample_type == st.ii)
    tbl1[which(tbl1$" " == st.ii), ncol(tbl1)] <- paste(dim(dat.st)[1],
                                                        length(unique(dat.st$indiv)), 
                                                        length(unique(dat.st$article)),
                                                        sep = " / ")
  }
  
  # Sample time 
  for (dpi.ii in unique(dat$dpi_idx)) {
    dat.sub.dpi <- subset(dat.sub, dpi_idx == dpi.ii)
    tbl1[which(tbl1$" " == dpi.ii), col_num] <- paste(dim(dat.sub.dpi)[1],
                                                      length(unique(dat.sub.dpi$indiv)), 
                                                      length(unique(dat.sub.dpi$article)),
                                                      sep = " / ")
    
    dat.dpi <- subset(dat, dpi_idx == dpi.ii)
    tbl1[which(tbl1$" " == dpi.ii), ncol(tbl1)] <- paste(dim(dat.dpi)[1],
                                                         length(unique(dat.dpi$indiv)), 
                                                         length(unique(dat.dpi$article)),
                                                         sep = " / ")
  }
  
  # Target gene 
  for (tg.ii in unique(dat$tg_total)) {
    dat.sub.tg <- subset(dat.sub, tg_total == tg.ii)
    tbl1[which(tbl1$" " == tg.ii), col_num] <- paste(dim(dat.sub.tg)[1],
                                                     length(unique(dat.sub.tg$indiv)), 
                                                     length(unique(dat.sub.tg$article)),
                                                     sep = " / ")
    
    dat.tg <- subset(dat, tg_total == tg.ii)
    tbl1[which(tbl1$" " == tg.ii), ncol(tbl1)] <- paste(dim(dat.tg)[1],
                                                        length(unique(dat.tg$indiv)), 
                                                        length(unique(dat.tg$article)),
                                                        sep = " / ")
  }
  
  # Culture assay 
  for (assay.ii in unique(dat$culture_assay)) {
    dat.sub.assay <- subset(dat.sub, culture_assay == assay.ii)
    tbl1[which(tbl1$" " == assay.ii), col_num] <- paste(dim(dat.sub.assay)[1],
                                                        length(unique(dat.sub.assay$indiv)), 
                                                        length(unique(dat.sub.assay$article)),
                                                        sep = " / ")
    
    dat.assay <- subset(dat, culture_assay == assay.ii)
    tbl1[which(tbl1$" " == assay.ii), ncol(tbl1)] <- paste(dim(dat.assay)[1],
                                                           length(unique(dat.assay$indiv)), 
                                                           length(unique(dat.assay$article)),
                                                           sep = " / ")
  }
  
  # Cell line 
  for (cell.ii in unique(dat$cell_line)) {
    dat.sub.cell <- subset(dat.sub, cell_line == cell.ii)
    tbl1[which(tbl1$" " == cell.ii), col_num] <- paste(dim(dat.sub.cell)[1],
                                                       length(unique(dat.sub.cell$indiv)), 
                                                       length(unique(dat.sub.cell$article)),
                                                       sep = " / ")
    
    dat.cell <- subset(dat, cell_line == cell.ii)
    tbl1[which(tbl1$" " == cell.ii), ncol(tbl1)] <- paste(dim(dat.cell)[1],
                                                          length(unique(dat.cell$indiv)), 
                                                          length(unique(dat.cell$article)),
                                                          sep = " / ")
  }
}

tbl1[nrow(tbl1), ncol(tbl1)] <- paste(nrow(dat),
                                      length(unique(dat$indiv)), 
                                      length(unique(dat$article)),
                                      sep = " / ")

tbl1[is.na(tbl1)] <- ""


colnames(tbl1)[colnames(tbl1) == "sgRNA"] <- "sgRNA & total RNA"
colnames(tbl1)[colnames(tbl1) == "Culture"] <- "Culture & either RNA"
colnames(tbl1)[colnames(tbl1) == "Total"] <- "All data"


## Save as csv  ----------------------------------------------------------------

write.csv(tbl1, "./outputs/tables/tbl1-data-summary.csv", 
          row.names = FALSE)





# Table S1: All data -----------------------------------------------------------

## Load data -------------------------------------------------------------------

dat <- read.csv("./data/clean-data.csv")


## Prep the table to be filled in ----------------------------------------------

cols <- c("Article", "N", "Species", "Sex", "Age_class", 
          "Exposure_route", "Exposure_dose", "Viral_isolate",
          "Sample_type", "Sample_location", "Sample_days", "total_RNA", "sgRNA",
          "Culture", "Viral_isolate")
tblS1 <- as.data.frame(matrix(nrow = 0, ncol = length(cols)),
                       col.names = cols)


## Fill in the table -----------------------------------------------------------

# First facet by article 
for (article.ii in sort(unique(dat$article))) {
  article_subset <- subset(dat, article == article.ii)
  
  # Facet by species
  for (species.ii in unique(article_subset$animal_species)) {
    species_subset <- subset(article_subset, animal_species == species.ii)
    
    # Facet by exposure route
    for (route.ii in unique(species_subset$inoc_route)) {
      route_subset <- subset(species_subset, inoc_route == route.ii)
      
      # Facet by exposure dose
      for (dose.ii in unique(route_subset$inoc_dose_total_pfu)) {
        dose_subset <- subset(route_subset, inoc_dose_total_pfu == dose.ii)
        
        # Convert dose.ii into log10 units for uniformity and star if reported as TCID50
        dose.ii <- round(log10(as.numeric(dose.ii)), digits = 2)
        
        if (unique(dose_subset$inoc_dose_units_rep) %in% c("TCID50", "CCID50")){
          dose.ii <- paste0(dose.ii, "*") 
        }
        
        # Just in case I end up facetting further
        final_subset <- dose_subset
        
        # Dealing with cell lines for studies including inf virus
        if (length(unique(final_subset$cell_line)) >= 2) {
          cell_line <- unique(final_subset$cell_line[!is.na(final_subset$cell_line)])
        }
        else if (is.na(unique(final_subset$cell_line))){
          cell_line <- "--"
        }
        else {
          cell_line <- unique(final_subset$cell_line)
        }
        
        if (NA %notin% unique(final_subset$culture_assay)) {
          if (unique(final_subset$culture_assay) == "TCID50"){
            cell_line <- paste0(cell_line, "+")
          }
        } 

      
        # Fill in the next row
        next_row <- list(Article = article.ii,
                         N = paste0(dim(final_subset)[1], " (", 
                                   length(unique(final_subset$indiv)), ")"),
                         Species = species.ii,
                         Sex = paste(sort(unique(final_subset$sex)), 
                                     collapse = ", "),
                         Age_class = paste(sort(unique(final_subset$age_class_stand)), 
                                           collapse = ", "),
                         Exposure_route = route.ii,
                         Exposure_dose = dose.ii,
                         Viral_isolate = str_remove(unique(final_subset$viral_strain), "SARS-CoV-2/human/"),
                         Sample_type = paste(sort(unique(final_subset$sample_type)), 
                                              collapse = ", "),
                         Sample_days = paste(sort(as.numeric(unique(final_subset$dpi_idx))), 
                                             collapse = ", "),
                         Sample_location = paste(unique(final_subset$organ_system)[order(match(unique(final_subset$organ_system), c("URT", "LRT", "GI", "Other")))], 
                                                 collapse = ", "),
                         Total_RNA = paste(unique(final_subset$tg_total[!is.na(final_subset$tg_total)]), 
                                           collapse = ", "),
                         sgRNA = paste0(paste(unique(final_subset$tg_sg[!is.na(final_subset$tg_sg)]), 
                                       collapse = ", "), " (",
                                       paste(unique(final_subset$tg_idx[!is.na(final_subset$tg_sg)]),
                                             collapse = ", "), ")"),
                         Culture = cell_line)
        
        tblS1 <- rbind(tblS1, next_row)
        
      }
    }
  }
}


## Add acronyms to decrease size -----------------------------------------------

tblS1 <- tblS1 %>% 
  mutate(Species = replace(Species, Species == "rhesus macaque", "RM")) %>%
  mutate(Species = replace(Species, Species == "cynomolgus macaque", "CM")) %>%
  mutate(Species = replace(Species, Species == "african green monkey", "AGM")) %>%
  mutate(Sample_type = replace(Sample_type, Sample_type == "Non-invasive", "NI")) %>%
  mutate(Sample_type = replace(Sample_type, Sample_type == "Invasive, Non-invasive", "I, NI")) %>%
  mutate(Sample_type = replace(Sample_type, Sample_type == "Invasive", "I")) %>%
  mutate(Sex = str_replace_all(Sex, "Female", "F")) %>%
  mutate(Sex = str_replace_all(Sex, "Male", "M")) %>%
  mutate(Sex = str_replace_all(Sex, "Unknown", "U")) %>%
  mutate(Exposure_dose = replace(Exposure_dose, is.na(Exposure_dose), "U")) %>%
  mutate(Viral_isolate = str_remove_all(Viral_isolate, "(no ID)")) %>%
  mutate(Viral_isolate = str_replace_all(Viral_isolate, "[?][?]", "U")) %>%
  mutate(Viral_isolate = str_remove_all(Viral_isolate, "[(][)]")) %>%
  mutate(sgRNA = str_replace_all(sgRNA, "[(][)]", "--")) %>%
  mutate(sgRNA = str_replace_all(sgRNA, "[(]NA[)]", "")) %>%
  mutate(Age_class = str_replace_all(Age_class, "Unknown", "U")) %>%
  mutate(Age_class = str_replace_all(Age_class, "Juvenile", "J")) %>%
  mutate(Age_class = str_replace_all(Age_class, "Adult", "A")) %>%
  mutate(Age_class = str_replace_all(Age_class, "Geriatric", "G")) %>%
  mutate(Age_class = replace(Age_class, Age_class == "A, J", "J, A"))

# One article without totRNA data
tblS1$Total_RNA[tblS1$Total_RNA == ""] <- "--"

## Remove "_" from column names
colnames(tblS1) <- str_replace_all(colnames(tblS1), "_", " ")
colnames(tblS1)[colnames(tblS1) == "Total RNA"] <- "Total RNA target gene"
colnames(tblS1)[colnames(tblS1) == "sgRNA"] <- "sgRNA target gene"



## Remove repeat article names -------------------------------------------------

for (row_num in 2:nrow(tblS1)){
  if (tblS1$Article[row_num - 1] == tblS1$Article[row_num]){
    tblS1$Article[row_num] <- ""
  }
}

View(tblS1)
  
 
## Save as csv  ----------------------------------------------------------------

write.csv(tblS1, "./outputs/tables/tblS1-data-summary.csv", 
          row.names = FALSE)

