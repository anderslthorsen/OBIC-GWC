# Written by Anders Lillevik Thorsen, November 2022
# Extract info to make list of OCD patients and relevant variables for vertex-wise analysis of symptom dimensions

# Load relevant packages
library(haven)
library(tidyverse)
library(openxlsx)

# Setup environment
rm(list=ls()) # Clears variables
options(scipen = 999) # Gives decimals rather than power for large number
setwd("S:/Project/OBIC/R")

# Load data
data <- read_sav(file="S:/Project/OBIC/dataset/Cortical_myelination_FINAL_20July22.sav")
data <- subset(data, Group == 1) # Select OCD patients only

extracted_columns <- c("Sex_Rel", "Agr_Check", "Contam_Clean", "Sym_Ordering", "Hoarding","Age", "Sex", "SurfaceHoles", "Site", "ASEG_ID") # Select which columns to extract

df <- data[ , extracted_columns] # Subset by specific columns
head(df)

df_complete_cases <- df[complete.cases(df), ]

# Prepare data to fit the design matrix in FSL Palm
df_complete_cases$Age <- scale(df_complete_cases$Age, scale = FALSE) # Mean-center
df_complete_cases$SurfaceHoles <- scale(df_complete_cases$SurfaceHoles, scale = FALSE) # Mean-center

head(df_complete_cases)

# Create dummy variables for each site, where site 8 is the reference. Note that no patients in site 2 had data on symptom dimensions.
df_complete_cases$Site_1 <- ifelse(df_complete_cases$Site == 1, 1, 0)
df_complete_cases$Site_3 <- ifelse(df_complete_cases$Site == 3, 1, 0)
df_complete_cases$Site_4 <- ifelse(df_complete_cases$Site == 4, 1, 0)
df_complete_cases$Site_5 <- ifelse(df_complete_cases$Site == 5, 1, 0)
df_complete_cases$Site_6 <- ifelse(df_complete_cases$Site == 6, 1, 0)
df_complete_cases$Site_7 <- ifelse(df_complete_cases$Site == 7, 1, 0)

mytable <- table(df_complete_cases$Site_1, df_complete_cases$Site_3, df_complete_cases$Site_4, df_complete_cases$Site_5, df_complete_cases$Site_6, df_complete_cases$Site_7)
ftable(mytable)

# Write list of IDs to use when merging surfaces into one 4D file 
df_complete_cases$ASEG_ID %>% write_lines("S:/Project/OBIC/R/symptom_dimensions/symptom_dim_path_to_files.txt", sep = "\r\n")

# Write design matrix
df_design_matrix <- df_complete_cases[ , c("Sex_Rel", "Agr_Check", "Contam_Clean", "Sym_Ordering", "Hoarding","Age", "Sex", "SurfaceHoles", "Site_1", "Site_3", "Site_4", "Site_5", "Site_6", "Site_7")]

write.table(df_design_matrix, file = "S:/Project/OBIC/R/symptom_dimensions/symptom_dim_design_matrix.txt", row.names = FALSE, col.names = FALSE)

# Write exchangeability block list
for (i in 1:nrow(df_complete_cases)){
 
  df_complete_cases$Site_recoded <- NA
  if(as.numeric(df_complete_cases$Site[i]) > 1){
    df_complete_cases$Site_recoded[i] <- as.numeric(df_complete_cases$Site[i])-1
    #print(as.numeric(df_complete_cases$Site[i])-1)
    print("test")
    }
}

write.xlsx(df_complete_cases, file = "S:/Project/OBIC/R/symptom_dimensions/df_complete_cases.xlsx")
