# Written by Anders Lillevik Thorsen, March 2022
# This script analyzes GWC and cortical thickness in the OBIC dataset

library(lme4)
library(skimr)
library(ggplot2)
library(haven)
library(effectsize)
library(parameters)
library(corrplot)

# Setup environment
rm(list=ls()) # Clears variables
options(scipen = 100) # Gives decimals rather than power for large number
setwd("S:/Project/OBIC/R")

# Load data
#data <- read_sav(file="Cortical_myelination_with_ICA_ROI_230322.sav")
data <- read_sav(file="S:/Project/OBIC/dataset/Cortical_myelination_FINAL_20July22.sav")


# Skim dataset to get an overview
#skim(data)

# Define names of all ROIs
ROI <- c("bankssts", "caudalanteriorcingulate", "caudalmiddlefrontal", "cuneus", "entorhinal", "fusiform", "inferiorparietal", "inferiortemporal", "isthmuscingulate", "lateraloccipital", "lateralorbitofrontal", "lingual", "medialorbitofrontal", "middletemporal", "parahippocampal", "paracentral", "parsopercularis", "parsorbitalis", "parstriangularis", "pericalcarine", "postcentral", "posteriorcingulate", "precentral", "precuneus", "rostralanteriorcingulate", "rostralmiddlefrontal", "superiorfrontal", "superiorparietal", "superiortemporal", "supramarginal", "frontalpole", "temporalpole", "transversetemporal", "insula")

# Calculate mean of two hemispheres for GWC and thickness
for (metric in c("GWC","thickness")){
  for (region in 1:length(ROI))
  {
    lh_ROI <- paste0("lh_", metric, "_", ROI[region])
    rh_ROI <- paste0("rh_", metric, "_", ROI[region])
    
    bil_ROI_name <- paste0("bil_", metric, "_", ROI[region])
    bil_ROI_vals <- (data[[lh_ROI]]+data[[rh_ROI]])/2
    
    data[bil_ROI_name] <- bil_ROI_vals
    
  }
}

# Remove left over variables
rm(metric, region, lh_ROI,rh_ROI, bil_ROI_name, bil_ROI_vals)

# Define table for storing results of linear mixed models
stats_LME_ROI<-data.frame(b_thickness=matrix(0,34,1),
                          stdB_thickness=matrix(0,34,1),
                          t_thickness=matrix(0,34,1), 
                          p_thickness=matrix(0,34,1), 
                          FDRp_thickness=matrix(0,34,1),
                          b_thickness_adjusted=matrix(0,34,1),
                          stdB_thickness_adjusted=matrix(0,34,1),
                          t_thickness_adjusted=matrix(0,34,1), 
                          p_thickness_adjusted=matrix(0,34,1), 
                          FDRp_thickness_adjusted=matrix(0,34,1),
                          b_GWC=matrix(0,34,1),
                          stdB_GWC=matrix(0,34,1),
                          t_GWC=matrix(0,34,1), 
                          p_GWC=matrix(0,34,1),
                          FDRp_GWC=matrix(0,34,1))

# Name rows by ROI
row.names(stats_LME_ROI)<-ROI

# Calculate age squared
data$Age_sq = data$Age*data$Age

# Do lmer for GWC and cortical thickness derived from Desikan-Killiany parcellation
for (metric in c("GWC", "thickness")){
  for(region in 1:length(ROI)){
    #add the ICA value 
    data$r<- data[[paste0("bil_", metric, "_", ROI[region])]]
    data$GWC_covariate <- data[[paste0("bil_GWC_", ROI[region])]]
  
    #run the LMER
    lmer.test <- lmer(r ~ Group + Sex + Age + SurfaceHoles + (1 | Site), data = data, REML = FALSE) 
    s=summary(lmer.test)
    
    #Add values to table
    stats_LME_ROI[region, paste0("t_", metric)]<-s$coefficients[2,4]
    stats_LME_ROI[region, paste0("p_", metric)]<-s$coefficients[2,5]
  
      # run the LMER for thickness wile adjusting for GWC
      if(metric=="thickness"){
      lmer.test.GWC_covariate <- lmer(r ~ Group + Sex + (1 | Site) + Age + SurfaceHoles + GWC_covariate, data = data, REML = FALSE) 
      s=summary(lmer.test.GWC_covariate)
      
      stats_LME_ROI[region, "t_thickness_adjusted"]<-s$coefficients[2,4]
      stats_LME_ROI[region, "p_thickness_adjusted"]<-s$coefficients[2,5]

      }
  }
}

# Remove left over variables
rm(s, metric, region, lmer.test, lmer.test.GWC_covariate)

# do fdr correction
stats_LME_ROI$FDRp_thickness <- p.adjust(stats_LME_ROI$p_thickness, method="fdr")
stats_LME_ROI$FDRp_GWC <- p.adjust(stats_LME_ROI$p_GWC, method="fdr")
stats_LME_ROI$FDRp_thickness_adjusted <- p.adjust(stats_LME_ROI$p_thickness_adjusted, method="fdr")

# Compare t- and p-values for adjusted and unadjusted models for cortical thickness
stats_LME_ROI$diff_t_ratio <- (stats_LME_ROI$t_thickness_adjusted/stats_LME_ROI$t_thickness)^2
stats_LME_ROI$diff_FDRp_thickness <- stats_LME_ROI$FDRp_thickness-stats_LME_ROI$FDRp_thickness_adjusted

# Save results as xlsx
write.xlsx(stats_LME_ROI, colNames = TRUE, rowNames = TRUE, 
           file = "S:/Project/OBIC/R/Results/cort_thickness_adjusted_for_GWC_results.xlsx")

# Relating components from ICA to regional GWC
# Make list of GWC variables
variables <- vector()

# Loop over ROIs and append to variable
for(i in 1:length(ROI)){
  for(hemi in c("lh", "rh")){
  variables <- append(variables, values = paste0(hemi, "_GWC_", ROI[i]))
  }
}

# Append ICA variables
ICAs <- c("ICA7_1",
          "ICA7_2",
          "ICA7_3",
          "ICA7_4",
          "ICA7_5",
          "ICA7_6",
          "ICA7_7",
          "ICA_Rerun7z_1",
          "ICA_Rerun7z_2",
          "ICA_Rerun7z_3", 
          "ICA_Rerun7z_4", 
          "ICA_Rerun7z_5", 
          "ICA_Rerun7z_6", 
          "ICA_Rerun7z_7")
variables <- append(variables, values = ICAs)

# Create new data frame containing ICA and GWC variables
df <- data[, variables]

# Correlating ICA and GWC
res <- cor(df, method = "pearson")
res <- round(res, 2)
res <- res[,69:82]
corrplot(res, type = "full", order = "original", tl.pos = "n")

# Remove left over variables
rm(variables, ICAs, i, hemi)
