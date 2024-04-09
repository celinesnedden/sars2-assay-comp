# Generates Table S11, which includes the assignment of inoculated vs. 
#     non-inoculated tissues, for all included exposure routes


# Prep environment -------------------------------------------------------------

## Install necessary packages  -------------------------------------------------

req_pkgs <- c("readxl", "stringr", "stringi", "dplyr")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependences = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))

`%notin%` <- Negate(`%in%`)  # for convenience


## Load the data  --------------------------------------------------------------

dat <- read.csv("./data/clean-data.csv")


# Assign consistent tissue names -----------------------------------------------

dat$tissue[str_detect(dat$sample_grp, "Lung")] <- "Lung"
dat$tissue[str_detect(dat$sample_grp, "Colon")] <- "Colon"
dat$tissue[str_detect(dat$sample_grp, "Small intestine")] <- "Small intestine"
dat$tissue[str_detect(dat$sample_grp, "Brain")] <- "Brain"
dat$tissue[str_detect(dat$sample_grp, "Bronch")] <- "Bronchus"
dat$tissue[str_detect(dat$sample_grp, "Eye|Ocular")] <- "Eye"
dat$tissue[str_detect(dat$sample_grp, "Nasal")] <- "Nose/Nasopharynx"
dat$tissue[str_detect(dat$sample_grp, "Naso")] <- "Nose/Nasopharynx"
dat$tissue[str_detect(dat$sample_grp, "Oro")] <- "Oropharynx"
dat$tissue[str_detect(dat$sample_grp, "Oral")] <- "Mouth"
dat$tissue[str_detect(dat$sample_grp, "Rectal|Rectum")] <- "Anus/Rectum"
dat$tissue[str_detect(dat$sample_grp, "Trachea")] <- "Trachea"
dat$tissue[str_detect(dat$sample_grp, "LN")] <- "Lymph Node"
dat$tissue[is.na(dat$tissue)] <- dat$sample_grp[is.na(dat$tissue)]


# Generate the table -----------------------------------------------------------

tbl_cols <- c("Inoculation_Route", "Inoculated_Locations", "NonInoculated_Locations")
tblSX <- as.data.frame(matrix(nrow = 0, ncol = length(tbl_cols)))
colnames(tblSX) <- tbl_cols

# Loop over exposure routes to fill in the table
for (route in unique(dat$inoc_route)) {
  route_sub <- subset(dat, inoc_route == route)
  inoc_samples <- paste(sort(unique(route_sub$tissue[route_sub$dpi_idx %in% c(1, 2)])), collapse = ", ")
  noninoc_samples <- paste(sort(unique(route_sub$tissue[route_sub$dpi_idx == 3])), collapse = ", ")
  next_row <- c(list(Inoculation_Route = route,
                     Inoculated_Locations = inoc_samples,
                     NonInoculated_Locations = noninoc_samples))
  tblSX <- rbind(tblSX, next_row)
}


# Format the table -------------------------------------------------------------

# Arrange order manually, from single to multi-inoculated procedures
inoc_order <- c("AE", "IT", "IN", "IG", "OC", "IT, IN", "IT, IN, OC", "IT, IN, OR, OC")
tblSX <- tblSX %>% arrange(factor(Inoculation_Route, levels = inoc_order))

# Change column names to be formatted as in the paper
colnames(tblSX) <- str_replace_all(colnames(tblSX), "_", " ")
colnames(tblSX)[colnames(tblSX) == "NonInoculated Locations"] <- "Non-Inoculated Locations"


# Save -------------------------------------------------------------------------

write.csv(tblSX, file = "./outputs/tables/tblS5-inoc-tissue-assignments.csv",
          row.names = FALSE)


