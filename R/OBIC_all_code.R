# Written by Anders Lillevik Thorsen, June 2022
# This script will do all needed data manipulations, analyses and plots in the OBIC GWC project when it is completed

# Load relevant packages
library(lme4)
library(ggplot2)
library(haven)
library(effectsize)
library(parameters)
library(performance)
library(ggseg)
library(tidyverse)
library(lmerTest)
library(EMAtools)
library(openxlsx)
library(RVIpkg)
library(sjPlot)

# Setup environment
rm(list=ls()) # Clears variables
options(scipen = 999) # Gives decimals rather than power for large number
setwd("S:/Project/OBIC/R")

# Load precalculated RVI dataset
data <- read_sav(file = "S:/Project/OBIC/dataset/Cortical_myelination_FINAL_RVI_reanalysis_05June23.sav")

# Load data
#data <- read_sav(file="S:/Project/OBIC/dataset/Cortical_myelination_FINAL_20July22.sav")

  # Run linear mixed model for all 7 ICA components
  data$r <- 0
  ROI.names_LME_OBIC <- data %>%
    select(contains("*ICA"), (contains("_Rerun")), -(contains("Medicated")))
  region_OBIC_LME <- names(ROI.names_LME_OBIC)
  
  stats_LME_ICA<-data.frame(t=matrix(0,7,1),
                            p=matrix(0,7,1),
                            StdB=matrix(0,7,1),
                            CI_low=matrix(0,7,1),
                            CI_high=matrix(0,7,1),
                            Cohensd=matrix(0,7,1),
                            ICC=matrix(0,7,1),
                            N_OCD=matrix(0,7,1),
                            N_HC=matrix(0,7,1))
  row.names(stats_LME_ICA)<-region_OBIC_LME
  
  for(region in 1:length(region_OBIC_LME)){
    #add the ICA value 
    data$r <- data[[paste0("ICA_Rerun7z_", region)]]
    
    #run the LMER
    ICA.test <- lmer(r ~ Group + Sex + Age + SurfaceHoles + (1 | Site), data = data, REML = FALSE) 
    s = summary(ICA.test)
    
    stdBeta <- effectsize(ICA.test)
    Cohensd <- lme.dscore(ICA.test, data=data, type = "lme4")
    icc <- icc(ICA.test)
    
    #Add values to table
    stats_LME_ICA[region,"t"]<-s$coefficients[2,4]
    stats_LME_ICA[region,"p"]<-s$coefficients[2,5]
    stats_LME_ICA[region,"StdB"]<-stdBeta$Std_Coefficient[2]
    stats_LME_ICA[region,"CI_low"]<-stdBeta$CI_low[2]
    stats_LME_ICA[region,"CI_high"]<-stdBeta$CI_high[2]
    stats_LME_ICA[region,"Cohensd"]<-Cohensd$d[1]
    stats_LME_ICA[region,"ICC"] <- icc$ICC_adjusted
    stats_LME_ICA[region,"N_OCD"]<-nrow(subset(data, Group == 1))
    stats_LME_ICA[region,"N_HC"]<-nrow(subset(data, Group == 0))
  }
  
  # Clean up temporary variables
  rm(region, region_OBIC_LME, ROI.names_LME_OBIC, s, ICA.test, stdBeta, Cohensd)
  
  # do fdr correction
  stats_LME_ICA$p_fdr<-p.adjust(stats_LME_ICA$p, method="fdr")

  # Round results to three decimals
  stats_LME_ICA <- round(stats_LME_ICA, digits = 3)

  # Export summary statistics to Excel
  write.xlsx(stats_LME_ICA, file = "Results/ICA_lmer_results.xlsx")
  
# Perform leave-one-site out validation for ICA component 3
  Leave_out_results<-data.frame(Left_out_site=matrix(0,8,1),
                        t=matrix(0,8,1),
                        p=matrix(0,8,1),
                        StdB=matrix(0,8,1),
                        CI_low=matrix(0,8,1),
                        CI_high=matrix(0,8,1),
                        Cohensd=matrix(0,8,1),
                        ICC=matrix(0,8,1),
                        N_OCD=matrix(0,8,1),
                        N_HC=matrix(0,8,1))
                            
  for(i in c(1:(length)(unique(data$Site)))){
    data_wo_i <- filter(data, Site != i)
    ICA.test <- lmer(ICA_Rerun7z_3 ~ Group + Sex + Age + SurfaceHoles + (1 | Site), data= data_wo_i, REML = FALSE)
  
    s <- summary(ICA.test)
    stdBeta <- effectsize(ICA.test)
    Cohensd <- lme.dscore(ICA.test, data=data_wo_i, type = "lme4")
    icc <- icc(ICA.test)
    N_OCD <- nrow(subset(data_wo_i, Group == 1))
    N_HC <- nrow(subset(data_wo_i, Group == 0))
    
    #Add values to table
    Leave_out_results[i,"Left_out_site"]<-i
    Leave_out_results[i,"t"]<-s$coefficients[2,4]
    Leave_out_results[i,"p"]<-s$coefficients[2,5]
    Leave_out_results[i,"StdB"]<-stdBeta$Std_Coefficient[2]
    Leave_out_results[i,"CI_low"]<-stdBeta$CI_low[2]
    Leave_out_results[i,"CI_high"]<-stdBeta$CI_high[2]
    Leave_out_results[i,"Cohensd"]<-Cohensd$d[1]
    Leave_out_results[i,"ICC"]<- icc$ICC_adjusted
    Leave_out_results[i,"N_OCD"]<-nrow(subset(data_wo_i, Group == 1))
    Leave_out_results[i,"N_HC"]<-nrow(subset(data_wo_i, Group == 0))
  
  }

  # Clean up temporary variables
  rm(i, icc, data_wo_i, s, ICA.test, stdBeta, Cohensd, N_HC, N_OCD)
  
  # Estimate Cohen's d for GWC per ROI to plot in ggseg

ggplot_ROIs <- c("bankssts","caudalanteriorcingulate","caudalmiddlefrontal","cuneus",
                 "entorhinal","frontalpole","fusiform","inferiorparietal","inferiortemporal","insula",
                 "isthmuscingulate","lateraloccipital","lateralorbitofrontal","lingual",
                 "medialorbitofrontal","middletemporal","paracentral","parahippocampal",
                 "parsopercularis","parsorbitalis","parstriangularis","pericalcarine",
                 "postcentral","posteriorcingulate","precentral","precuneus",
                 "rostralanteriorcingulate","rostralmiddlefrontal",
                 "superiorfrontal","superiorparietal","superiortemporal",
                 "supramarginal","temporalpole","transversetemporal")

# Set up dataframe to hold results which are to be plotted by ggseg
results_Cohensd <- tibble(
  region = rep(c(ggplot_ROIs_spaced), 2),
  mean = as.double(NA),
  pval = as.double(NA),
  hemi = c(rep("left", 34), rep("right", 34))
)

# Estimate standardized mean GWC per ROI
for(hemi in c("lh","rh")){ # Loop over hemispheres
  i=0 # Sets starting point for counter
  if (hemi == "lh"){ 
    i=i+0 # Set to zero as left is first in the dataframe
  } else if (hemi == "rh"){
    i=i+34
  }
  
    for(region in ggplot_ROIs){ # Loop over regions
      
      i=i+1
      
      data$r <- data[[paste0(hemi,"_GWC_",region)]]
      lmer.test <- lmer(r ~ Group + Sex + Age + SurfaceHoles + (1 | Site), data = data, REML = FALSE) 
      Cohensd <- lme.dscore(lmer.test, data=data, type = "lme4")
      s <- summary(lmer.test)
  
      results_Cohensd[i,"pval"]<-s$coefficients[2,5]
      results_Cohensd[i,"mean"] <- Cohensd$d[1]
      
    }
}

# Export summary statistics to Excel
write.xlsx(results_Cohensd, file = "Results/GWC_ROI_CohensD_results.xlsx")

# Plot standardized mean difference between groups using ggseg
#tiff("../Figures/GWC_DK_OCD_vs_HC.tiff", compression = "zip", res = 288, width = 1440, height = 750)
svglite(filename = "../Figures/GWC_DK_OCD_vs_HC.svg")
results_Cohensd %>%
  ggplot() +
  geom_brain(atlas = dk,
             position=position_brain(hemi ~ side),
             colour ="black",
             aes(fill = mean)) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  labs(fill="Cohens's d") +
  ggtitle("Regional GWC in OCD versus HC") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())
invisible(dev.off())

    #results_Cohensd %>%
    #  ggseg(mapping=aes(fill=mean), position="stacked", colour="black")

results_Cohensd_GM <- tibble(
  region = rep(c(ggplot_ROIs_spaced), 2),
  mean = as.double(NA),
  hemi = c(rep("left", 34), rep("right", 34))
)

# Estimate standardized mean gray matter signal intensity per ROI
for(hemi in c("lh","rh")){ # Loop over hemispheres
  i=0 # Sets starting point for counter
  if (hemi == "lh"){ 
    i=i+0 # Set to zero as left is first in the dataframe
  } else if (hemi == "rh"){
    i=i+34
  }
  
  for(region in ggplot_ROIs){ # Loop over regions
    
    i=i+1
    
    data$r <- data[[paste0(hemi,"_GM_",region)]]
    lmer.test <- lmer(r ~ Group + Sex + Age + SurfaceHoles + (1 | Site), data = data, REML = FALSE) 
    Cohensd <- lme.dscore(lmer.test, data=data, type = "lme4")
    
    results_Cohensd_GM[i,"mean"] <- Cohensd$d[1]

  }
}

# Plot standardized mean difference between groups using ggseg
results_Cohensd_GM %>%
  ggplot() +
  geom_brain(atlas = dk,
             position=position_brain(hemi ~ side),
             colour ="black",
             aes(fill = mean)) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  labs(fill="Cohens's d relative to HC") +
  ggtitle("Regional gray matter signal intensity in OCD versus HC") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())

results_Cohensd_WM <- tibble(
  region = rep(c(ggplot_ROIs_spaced), 2),
  mean = as.double(NA),
  hemi = c(rep("left", 34), rep("right", 34))
)

# Estimate standardized mean white matter signal intensity per ROI
for(hemi in c("lh","rh")){ # Loop over hemispheres
  i=0 # Sets starting point for counter
  if (hemi == "lh"){ 
    i=i+0 # Set to zero as left is first in the dataframe
  } else if (hemi == "rh"){
    i=i+34
  }
  
  for(region in ggplot_ROIs){ # Loop over regions
    
    i=i+1
    
    data$r <- data[[paste0(hemi,"_WM_",region)]]
    lmer.test <- lmer(r ~ Group + Sex + Age + SurfaceHoles + (1 | Site), data = data, REML = FALSE) 
    Cohensd <- lme.dscore(lmer.test, data=data, type = "lme4")
    
    results_Cohensd_WM[i,"mean"] <- Cohensd$d[1]
    
  }
}

# Plot standardized mean difference between groups using ggseg
results_Cohensd_WM %>%
  ggplot() +
  geom_brain(atlas = dk,
             position=position_brain(hemi ~ side),
             colour ="black",
             aes(fill = mean)) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  labs(fill="Cohens's d relative to HC") +
  ggtitle("Regional white matter signal intensity in OCD versus HC") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())

# Estimate difference in standardized mean group difference (Cohen's d) between GM and WM signal intensity
results_Cohensd_GM_WM_diff <- tibble(
  region = rep(c(ggplot_ROIs_spaced), 2),
  mean = as.double(NA),
  hemi = c(rep("left", 34), rep("right", 34))
)

# 
results_Cohensd_GM_WM_diff$mean <- results_Cohensd_GM$mean-results_Cohensd_WM$mean

# Plot difference in standardized mean difference between groups using ggseg
results_Cohensd_WM %>%
  ggplot() +
  geom_brain(atlas = dk,
             position=position_brain(hemi ~ side),
             colour ="black",
             aes(fill = mean)) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  labs(fill="Cohens's d relative to HC") +
  ggtitle("Regional white matter signal intensity in OCD versus HC") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())

# Clean up temporary variables
rm(i, hemi, group, region, ggplot_ROIs, ggplot_ROIs_spaced, Cohensd, lmer.test)

#####
# Estimate regional vulnerability index (Kochunov et al., 2022, Hum Brain Mapp) and relate to GWC
#####

# Load ENIGMA group differences from RVIpkg, add OCD data from Boedhoe, 2020 Am. J. Psychiatry
ENIGMA_subcortical <- RVIpkg::EP.Subcortical
ENIGMA_cortical <- RVIpkg::EP.GM
ENIGMA_cortical[nrow(ENIGMA_cortical) + 1,1] = "temporalpole" # Add temporal pole as this is missing in RVIpkg

ENIGMA_subcortical$OCD <- c(.11, -0.05, 0.01, 0, 0.09, -0.09, -0.06, -0.03)
# Add estimates from Boedhoe 2020 Am J Psychiatry without OBIC participants (aka unbiased estimates in this context)
ENIGMA_subcortical$OCD_reanalysis <- c(0.04, # Lateral ventricle
                                       0.02,  # Thalamus
                                       0.03,  # Caudate
                                       0.09,  # Putamen
                                       0.02,  # Pallidum
                                       -0.04,  # Hippocampus
                                       -0.07, # Amygdala
                                       -0.03) # Accumbens

cor(ENIGMA_subcortical$OCD, ENIGMA_subcortical$OCD_reanalysis) # r=0.5634289 - okayish overlap

ENIGMA_cortical$OCD <- c(-0.05, # Banksts
                         0.01, # caudalanteriorcingulate
                         -0.08, # caudalmiddlefrontal
                         0.01, # cuneus
                         -0.02, # entorhinal
                         -0.08, # fusiform
                         -0.11, # inferiorparietal
                         -0.10, # inferiortemporal
                         -0.04, # isthmuscingulate
                         -0.07, # lateraloccipital
                         -0.07, # lateralorbitofrontal
                         -0.02, # lingual
                         -0.10, # medialorbitofrontal
                         -0.10, # middletemporal
                         -0.03, # parahippocampal
                         -0.01, # paracentral
                         -0.06, # parsopercularis
                         -0.06, # parsorbitalis
                         -0.06, # parstriangularis
                         0.02, # pericalcarine
                         0.00, # postcentral
                         -0.07, # posteriorcingulate
                         -0.04, # precentral
                         -0.10, # precuneus
                         0.00, # rostralanteriorcingulate
                         -0.08, # rostralmiddlefrontal
                         -0.05, # superiorfrontal
                         -0.04, # superiorparietal
                         -0.02, # superiortemporal
                         -0.08, # supramarginal
                         0.00, # frontalpole
                         -0.01, # transversetemporal
                         -0.05, # insula
                         0.01) # temporalpole

# Add new, unbiased estimates from Boedhoe 2020 Am J Psychiatry without OBIC participants
ENIGMA_cortical$OCD_reanalysis <- c(-0.05, # Banksts
                                    -0.02, # caudalanteriorcingulate
                                    -0.10, # caudalmiddlefrontal
                                    -0.03, # cuneus
                                    -0.04, # entorhinal
                                    -0.13, # fusiform
                                    -0.13, # inferiorparietal
                                    -0.12, # inferiortemporal
                                    -0.07, # isthmuscingulate
                                    -0.06, # lateraloccipital
                                    -0.11, # lateralorbitofrontal
                                    -0.06, # lingual
                                    -0.10, # medialorbitofrontal
                                    -0.13, # middletemporal
                                    -0.06, # parahippocampal
                                    0.00, # paracentral
                                    -0.09, # parsopercularis
                                    -0.07, # parsorbitalis
                                    -0.06, # parstriangularis
                                    0.06, # pericalcarine
                                    -0.02, # postcentral
                                    -0.07, # posteriorcingulate
                                    -0.06, # precentral
                                    -0.11, # precuneus
                                    -0.06, # rostralanteriorcingulate
                                    -0.09, # rostralmiddlefrontal
                                    -0.07, # superiorfrontal
                                    -0.04, # superiorparietal
                                    -0.01, # superiortemporal
                                    -0.06, # supramarginal
                                    -0.01, # frontalpole
                                    -0.03, # transversetemporal
                                    -0.09, # insula
                                    -0.01) # temporalpole

# Change order between insula (row 33) and temporal pole (row 34) so that the order is temporal pole (32), transeversetemporal (33) and insula (34)
sorted_ENIGMA_cortical <- ENIGMA_cortical
sorted_ENIGMA_cortical[32,] <- ENIGMA_cortical[34,] # should be temporal pole
sorted_ENIGMA_cortical[33,] <- ENIGMA_cortical[32,] # should be transversetemporal
sorted_ENIGMA_cortical[34,] <- ENIGMA_cortical[33,] # should be insula
#ENIGMA_cortical <- sorted_ENIGMA_cortical

cor(ENIGMA_cortical$OCD, ENIGMA_cortical$OCD_reanalysis) # r=0.87386 - good overlap

ENIGMA_cortical$SSD[34] <- -0.241 # Add cortical thickness for temporal pole for SSD from Van Erp, 2018, Biol Psychiatry

# Define names of all cortical ROIs
ROI <- c("bankssts", "caudalanteriorcingulate", "caudalmiddlefrontal", "cuneus", "entorhinal", "fusiform", "inferiorparietal", "inferiortemporal", "isthmuscingulate", "lateraloccipital", "lateralorbitofrontal", "lingual", "medialorbitofrontal", "middletemporal", "parahippocampal", "paracentral", "parsopercularis", "parsorbitalis", "parstriangularis", "pericalcarine", "postcentral", "posteriorcingulate", "precentral", "precuneus", "rostralanteriorcingulate", "rostralmiddlefrontal", "superiorfrontal", "superiorparietal", "superiortemporal", "supramarginal", "frontalpole", "temporalpole", "transversetemporal", "insula")

# Calculate mean of two hemispheres for cortical thickness
for (region in 1:length(ROI)){
  lh_ROI <- paste0("lh_thickness", "_", ROI[region])
  rh_ROI <- paste0("rh_thickness", "_", ROI[region])
  
  bil_ROI_name <- paste0("bil_thickness", "_", ROI[region])
  bil_ROI_vals <- (data[[lh_ROI]]+data[[rh_ROI]])/2
  
  data[bil_ROI_name] <- bil_ROI_vals
}

# Define names of all subcortical ROIs
data$lh_vol_lateralventricle <- data$LeftLateralVentricle
data$rh_vol_lateralventricle <- data$RightLateralVentricle

ROI <- c("lateralventricle","Thalamus","Caudate","Putamen","Pallidum","Hippocampus","Amygdala","Accumbensarea")

# Calculate mean of two hemispheres for subcortical volume
for (region in 1:length(ROI)){
  lh_ROI <- paste0("lh_vol", "_", ROI[region])
  rh_ROI <- paste0("rh_vol", "_", ROI[region])
  
  bil_ROI_name <- paste0("bil_vol", "_", ROI[region])
  bil_ROI_vals <- (data[[lh_ROI]]+data[[rh_ROI]])/2
  
  data[bil_ROI_name] <- bil_ROI_vals
}

cortical_column1 <- which(colnames(data) == "bil_thickness_bankssts")
cortical_column34 <- cortical_column1+33
colnames(data)[cortical_column1:cortical_column34] # check order

subcortical_column1 <- which(colnames(data) == "bil_vol_lateralventricle")
subcortical_column8 <- subcortical_column1+7
colnames(data)[subcortical_column1:subcortical_column8] # check order

# Run RVIpkg for cortical thickness derived from OCD
#RVI_cortical <- RVI_func(ID='BIDS', DXcontrol='Group==0', covariates=c('Age','Sex'), resp.range=c(cortical_column1:cortical_column34),EP=ENIGMA_cortical$OCD, data=data)
#data$RVI_cortical <- RVI_cortical$RVI$RVI
#data$RVI_cortical

# Run RVIpkg for cortical thickness derived from OCD - reanalysis corrected by removing OBIC participants and correction order of columns in resp.range
RVI_cortical <- RVI_func(ID='BIDS', DXcontrol='Group==0', covariates=c('Age','Sex'), resp.range=c(cortical_column1:cortical_column34),EP=sorted_ENIGMA_cortical$OCD_reanalysis, data=data)
data$RVI_cortical_reanalysis <- RVI_cortical$RVI$RVI
data$RVI_cortical_reanalysis

cor(data$RVI_cortical, data$RVI_cortical_reanalysis) # r=0.804616, relatively okay overlap even though the previous RVI was wrong (three columns switched) and biased by including OBIC in both data and reference

# Run RVIpkg for subcortical volume derived from OCD
RVI_subcortical <- RVI_func(ID='BIDS', DXcontrol='Group==0', covariates=c('Age','Sex', 'EstimatedTotalIntraCranialVol'), resp.range=c(subcortical_column1:subcortical_column8),EP=ENIGMA_subcortical$OCD_reanalysis, data=data)
data$RVI_subcortical_reanalysis <- RVI_subcortical$RVI$RVI
cor(data$RVI_subcortical, data$RVI_subcortical_reanalysis) # r=0.7226277, relatively okay overlap given bias by including OBIC in both data and reference

# Calculate mean RVI across the cortex and subcortex
data$RVI_mean_reanalysis <- (data$RVI_cortical_reanalysis+data$RVI_subcortical_reanalysis)/2
cor(data$RVI_mean, data$RVI_mean_reanalysis) # r=0.7359412, relatively okay overlap given bias by including OBIC in both data and reference

# Save dataset with RVi and mean bilateral cortical thickness/subcortical volume
write_sav(data, 'S:/Project/OBIC/dataset/Cortical_myelination_FINAL_RVI_reanalysis_05June23.sav')

# Load precalculated RVI dataset
#data <- read_sav(file = "S:/Project/OBIC/dataset/Cortical_myelination_FINAL_RVI_07Aug22.sav")

# Test Group difference in RVI
lmer_RVI_cort <-lmer(RVI_cortical_reanalysis ~ Group + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE) 
summary(lmer_RVI_cort)
plot_model(lmer_RVI_cort, title = "Cortical RVI", type = "std", show.values = TRUE)
lme.dscore(lmer_RVI_cort, data=data, type = "lme4")

lmer_RVI_subcort <-lmer(RVI_subcortical_reanalysis ~ Group + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE) 
summary(lmer_RVI_subcort)
plot_model(lmer_RVI_subcort, title = "Subcortical RVI", type = "std", show.values = TRUE)
lme.dscore(lmer_RVI_subcort, data=data, type = "lme4")

lmer_RVI_mean <-lmer(RVI_mean_reanalysis ~ Group + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE) 
summary(lmer_RVI_mean)
plot_model(lmer_RVI_mean, title = "Mean cortical and subcortical RVI", type = "std", show.values = TRUE)
lme.dscore(lmer_RVI_mean, data=data, type = "lme4")

# Estimate mean (SD) RVI scores per group
data_OCD <- subset(data, Group == 1)
summary(data_OCD$RVI_mean_reanalysis)
sd(data_OCD$RVI_mean_reanalysis)
summary(data_OCD$RVI_mean)
sd(data_OCD$RVI_mean)

data_HC <- subset(data, Group == 0)
summary(data_HC$RVI_mean_reanalysis)
sd(data_HC$RVI_mean_reanalysis)
summary(data_HC$RVI_mean)
sd(data_HC$RVI_mean)

lmer_RVI_mean <-lmer(ICA_Rerun7z_3 ~ Group + RVI_mean + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE) 
summary(lmer_RVI_mean)
plot_model(lmer_RVI_mean, title = "Mean cortical and subcortical RVI", type = "std", show.values = TRUE)
lme.dscore(lmer_RVI_mean, data=data, type = "lme4")

# Run linear mixed model for RVI mean, cortical and subcortical to compare groups
data$r<-0
ROI.names_LME_OBIC <- data %>%
  select(contains("*ICA"), (contains("_Rerun")), -(contains("Medicated")))
region_OBIC_LME<-names(ROI.names_LME_OBIC)
stats_LME_RVI<-data.frame(t=matrix(0,7,1),
                          p=matrix(0,7,1),
                          StdB=matrix(0,7,1),
                          CI_low=matrix(0,7,1),
                          CI_high=matrix(0,7,1),
                          Cohensd=matrix(0,7,1),
                          ICC=matrix(0,7,1),
                          N_OCD=matrix(0,7,1),
                          N_HC=matrix(0,7,1))
row.names(stats_LME_RVI) <- region_OBIC_LME

for(region in 1:length(region_OBIC_LME)){
  
  print(region_OBIC_LME[region])
  #add the ICA value 
  data$r <- data[[region_OBIC_LME[region]]]
  
  #run the LMER
  RVI.test <- lmer(r ~ RVI_mean_reanalysis + Group + Sex + Age + SurfaceHoles + (1 | Site), data = data, REML = FALSE) 
  s=summary(RVI.test)
  
  stdBeta <- effectsize(RVI.test)
  Cohensd <- lme.dscore(RVI.test, data=data, type = "lme4")
  icc <- icc(RVI.test)
  
  #Add values to table
  stats_LME_RVI[region,"t"]<-s$coefficients[2,4]
  stats_LME_RVI[region,"p"]<-s$coefficients[2,5]
  stats_LME_RVI[region,"StdB"]<-stdBeta$Std_Coefficient[2]
  stats_LME_RVI[region,"CI_low"]<-stdBeta$CI_low[2]
  stats_LME_RVI[region,"CI_high"]<-stdBeta$CI_high[2]
  stats_LME_RVI[region,"Cohensd"]<-Cohensd$d[1]
  stats_LME_RVI[region,"ICC"] <- icc$ICC_adjusted
  stats_LME_RVI[region,"N_OCD"]<-nrow(subset(data, Group == 1))
  stats_LME_RVI[region,"N_HC"]<-nrow(subset(data, Group == 0))
}

# do fdr correction
stats_LME_RVI$p_fdr<-p.adjust(stats_LME_RVI$p, method="fdr")

# Round results to three decimals
stats_LME_RVI <- round(stats_LME_RVI, digits = 3)

# Export summary statistics to Excel
write.xlsx(stats_LME_RVI, file = "Results/RVI_reanalysis_lmer_results.xlsx")

# Clean up temporary variables
rm(region, region_OBIC_LME, ROI.names_LME_OBIC, s, RVI.test, stdBeta, Cohensd, RVI_cortical, RVI_subcortical)

# Run linear mixed model to relate mean RVI for schizophrenia to all 7 ICA components
data$r<-0
ROI.names_LME_OBIC<-data %>%
  select(contains("*ICA"), (contains("_Rerun")), -(contains("Medicated")))
region_OBIC_LME<-names(ROI.names_LME_OBIC)

stats_LME_ICA<-data.frame(t=matrix(0,7,1),
                          p=matrix(0,7,1),
                          StdB=matrix(0,7,1),
                          CI_low=matrix(0,7,1),
                          CI_high=matrix(0,7,1),
                          Cohensd=matrix(0,7,1),
                          ICC=matrix(0,7,1),
                          N_OCD=matrix(0,7,1),
                          N_HC=matrix(0,7,1))
row.names(stats_LME_ICA)<-region_OBIC_LME

for(region in 1:length(region_OBIC_LME)){
  #add the ICA value 
  data$r <- data[[paste0("ICA_Rerun7z_", region)]]
  
  #run the LMER
  ICA.test <- lmer(r ~ Group + RVI_mean + Sex + Age + SurfaceHoles + (1 | Site), data = data, REML = FALSE) 
  s=summary(ICA.test)
  
  stdBeta <- effectsize(ICA.test)
  Cohensd <- lme.dscore(ICA.test, data=data, type = "lme4")
  icc <- icc(ICA.test)
  
  #Add values to table
  stats_LME_ICA[region,"t"]<-s$coefficients[3,4]
  stats_LME_ICA[region,"p"]<-s$coefficients[3,5]
  stats_LME_ICA[region,"StdB"]<-stdBeta$Std_Coefficient[3]
  stats_LME_ICA[region,"CI_low"]<-stdBeta$CI_low[3]
  stats_LME_ICA[region,"CI_high"]<-stdBeta$CI_high[3]
  stats_LME_ICA[region,"Cohensd"]<-Cohensd$d[2]
  stats_LME_ICA[region,"ICC"] <- icc$ICC_adjusted
  stats_LME_ICA[region,"N_OCD"]<-nrow(subset(data, Group == 1))
  stats_LME_ICA[region,"N_HC"]<-nrow(subset(data, Group == 0))
}

# Clean up temporary variables
rm(region, region_OBIC_LME, ROI.names_LME_OBIC, s, ICA.test, stdBeta, Cohensd)

# do fdr correction
stats_LME_ICA$p_fdr<-p.adjust(stats_LME_ICA$p, method="fdr")

# Round results to three decimals
stats_LME_ICA <- round(stats_LME_ICA, digits = 3)

# Export summary statistics to Excel
write.xlsx(stats_LME_ICA, file = "Results/RVI_for_SSD_ICA_lmer_results.xlsx")

  # Plot group distributions of RVI using raincloud plots

data$Group_category <- factor(data$Group, labels = c("HC", "OCD"))
  
data %>%
    ggplot(aes(x = Group, y = RVI_mean, fill = Group)) +
    
    # add half-violin from ggdist package
    ggdist::stat_halfeye(
      ## custom bandwidth
      adjust = 0.5,
      ##move geom to the right
      justification = -.1,
      ## remove slab interval
      .width = 0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = .08
    )

data %>%
  ggplot(aes(x = Group, y = RVI_subcortical, fill = Group)) +
  
  # add half-violin from ggdist package
  ggdist::stat_halfeye(
    ## custom bandwidth
    adjust = 0.5,
    ##move geom to the right
    justification = -.1,
    ## remove slab interval
    .width = 0,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .08
  )

#tiff("../Figures/RVI_OCD_vs_HC.tiff", compression = "zip", res = 288, width = 1800, height = 1440)
svglite(filename = "../Figures/RVI_reanalysis_OCD_vs_HC.svg")
data %>%
  ggplot(aes(x = Group_category, y = RVI_mean_reanalysis, fill = Group_category)) +
  ylab("Regional Vulnerability  Index") +
  xlab("Groups") +
  labs(fill = "") +
  theme(legend.position = "none", axis.title = element_text(size = 15), axis.text = element_text(size=15)) +
  # add half-violin from ggdist package
  ggdist::stat_halfeye(
    ## custom bandwidth
    adjust = 0.5,
    ##move geom to the right
    justification = -.06,
    ## remove slab interval
    .width = 0,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .08
  )
invisible(dev.off())

hist(data$RVI_subcortical_reanalysis)

#####
# Relate mean RVI to GWC per ROI from DK atlas
#####
# Plot GWC values using ggseg

detach("package:psych")

# Reduce data frame to variables for group and GWC
#ROI.GWC<-data %>%
#  select(contains("Group"), (contains("h_GWC_")))
#ROI.GWC <- subset(ROI.GWC, select = -c(lh_GWC_unknown, rh_GWC_unknown))

ggplot_ROIs_spaced <- brain_regions(dk)[-4] # Get labels from ggseg and removes corpus callosum from list
#ggplot_ROIs_spaced <- ggplot_ROIs_spaced[-4] # Remove corpus callosum from list

# Estimate the relation between mean RVI as Cohen's d and standardized beta for GWC per ROI to plot in ggseg

ggplot_ROIs <- c("bankssts","caudalanteriorcingulate","caudalmiddlefrontal","cuneus",
                 "entorhinal","frontalpole","fusiform","inferiorparietal","inferiortemporal","insula",
                 "isthmuscingulate","lateraloccipital","lateralorbitofrontal","lingual",
                 "medialorbitofrontal","middletemporal","paracentral","parahippocampal",
                 "parsopercularis","parsorbitalis","parstriangularis","pericalcarine",
                 "postcentral","posteriorcingulate","precentral","precuneus",
                 "rostralanteriorcingulate","rostralmiddlefrontal",
                 "superiorfrontal","superiorparietal","superiortemporal",
                 "supramarginal","temporalpole","transversetemporal")

# Set up dataframe to hold results which are to be plotted by ggseg
results_Cohensd <- tibble(
  region = rep(c(ggplot_ROIs_spaced), 2),
  mean = as.double(NA),
  hemi = c(rep("left", 34), rep("right", 34))
)

# Estimate mean GWC per ROI per group
for(hemi in c("lh","rh")){ # Loop over hemispheres
  i=0 # Sets starting point for counter
  if (hemi == "lh"){ 
    i=i+0 # Set to zero as left is first in the dataframe
  } else if (hemi == "rh"){
    i=i+34
  }
  
  for(region in ggplot_ROIs){ # Loop over regions
    
    i=i+1
    
    data$r <- data[[paste0(hemi,"_GWC_",region)]]
    lmer.test <- lmer(r ~ RVI_mean_reanalysis + Group + Sex + Age + SurfaceHoles + (1 | Site), data = data, REML = FALSE) 
    stdBeta <- effectsize(lmer.test)
    #Cohensd <- lme.dscore(lmer.test, data=data, type = "lme4")
    
    # print(i)
    # print(Cohensd$d[1])
    results_Cohensd[i,"mean"] <- stdBeta$Std_Coefficient[2]
    #results_Cohensd[i,"mean"] <- Cohensd$d[1]
    
    
    #print(i) # For debugging
    #print(paste0(hemi,"_GWC_",region))
  }
}

# Export summary statistics to Excel
write.xlsx(results_Cohensd, file = "Results/RVI_reanalysis_stdb_GWC.xlsx")

# Plot standardized mean difference between groups using ggseg
svglite(filename = "../Figures/RVI_reanalysis_stdB_GWC.svg")
results_Cohensd %>%
  ggplot() +
  geom_brain(atlas = dk,
             position=position_brain(hemi ~ side),
             colour ="black",
             aes(fill = mean)) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  labs(fill="Standardized Beta") +
  ggtitle("Regional GWC is related to Regional Vulnerability Index for OCD") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())
invisible(dev.off())

#results_Cohensd %>%
#  ggseg(mapping=aes(fill=mean), position="stacked", colour="black")

# Clean up temporary variables
rm(i, hemi, group, region, ggplot_ROIs, ggplot_ROIs_spaced, Cohensd, lmer.test)

# Write RVI to Excel
writeClipboard(as.character(RVI_excel))
writeClipboard(as.character(scale(data$Age, center = TRUE, scale = TRUE)))


#####
# Relate GWC component 3 to symptom dimensions in OCD patients
#####

  symptomdim_RVI <- lmer(ICA_Rerun7z_3 ~ Agr_Check + Contam_Clean + Sym_Ordering + Sex_Rel + Hoarding + Sex + Age + SurfaceHoles + (1 | Site), data = subset(data, data$Group == 1), REML = FALSE) 
  summary(symptomdim_RVI)
  
  detach("package:psych")
  
  # Reduce data frame to variables for group and GWC
  #ROI.GWC<-data %>%
  #  select(contains("Group"), (contains("h_GWC_")))
  #ROI.GWC <- subset(ROI.GWC, select = -c(lh_GWC_unknown, rh_GWC_unknown))
  
  ggplot_ROIs_spaced <- brain_regions(dk)[-4] # Get labels from ggseg and removes corpus callosum from list
  #ggplot_ROIs_spaced <- ggplot_ROIs_spaced[-4] # Remove corpus callosum from list
  
  # Estimate the relation between mean RVI as Cohen's d and standardized beta for GWC per ROI to plot in ggseg
  
  ggplot_ROIs <- c("bankssts","caudalanteriorcingulate","caudalmiddlefrontal","cuneus",
                   "entorhinal","frontalpole","fusiform","inferiorparietal","inferiortemporal","insula",
                   "isthmuscingulate","lateraloccipital","lateralorbitofrontal","lingual",
                   "medialorbitofrontal","middletemporal","paracentral","parahippocampal",
                   "parsopercularis","parsorbitalis","parstriangularis","pericalcarine",
                   "postcentral","posteriorcingulate","precentral","precuneus",
                   "rostralanteriorcingulate","rostralmiddlefrontal",
                   "superiorfrontal","superiorparietal","superiortemporal",
                   "supramarginal","temporalpole","transversetemporal")
  
  # Set up dataframe to hold results which are to be plotted by ggseg
  results_Cohensd <- tibble(
    region = rep(c(ggplot_ROIs_spaced), 2),
    mean = as.double(NA),
    hemi = c(rep("left", 34), rep("right", 34))
  )
  
  # Estimate mean GWC per ROI per group
  for(hemi in c("lh","rh")){ # Loop over hemispheres
    i=0 # Sets starting point for counter
    if (hemi == "lh"){ 
      i=i+0 # Set to zero as left is first in the dataframe
    } else if (hemi == "rh"){
      i=i+34
    }
    
    for(region in ggplot_ROIs){ # Loop over regions
      
      i=i+1
      
      data$r <- data[[paste0(hemi,"_GWC_",region)]]
      lmer.test <- lmer(r ~ Agr_Check + Contam_Clean + Sym_Ordering + Sex_Rel + Hoarding + Sex + Age + SurfaceHoles + (1 | Site), data = subset(data, data$Group == 1), REML = FALSE) 
      stdBeta <- effectsize(lmer.test)
      #Cohensd <- lme.dscore(lmer.test, data=data, type = "lme4")
      
      # print(i)
      # print(Cohensd$d[1])
      results_Cohensd[i,"mean"] <- stdBeta$Std_Coefficient[5]
      #results_Cohensd[i,"mean"] <- Cohensd$d[1]
      
      
      #print(i) # For debugging
      #print(paste0(hemi,"_GWC_",region))
    }
  }
  
  # Export summary statistics to Excel
  write.xlsx(results_Cohensd, file = "Results/GWC_ROI_CohensD_results.xlsx")
  
  # Plot standardized mean difference between groups using ggseg
  results_Cohensd %>%
    ggplot() +
    geom_brain(atlas = dk,
               position=position_brain(hemi ~ side),
               colour ="black",
               aes(fill = mean)) +
    scale_fill_gradient2(low="blue", mid = "white", high="red") +
    labs(fill="Standardized Beta") +
    ggtitle("Regional GWC is related to sexual/religious OCD symptoms") +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())


# DELETE LATER
  # Run RVIpkg using the pattern derived from schizophrenia 
  RVI_cortical_SSD <- RVI_func(ID='BIDS', DXcontrol='Group==0', covariates=c('Age','Sex'), resp.range=c(cortical_column1:cortical_column34),EP=ENIGMA_cortical$SSD, data=data)
  RVI_subcortical_SSD <- RVI_func(ID='BIDS', DXcontrol='Group==0', covariates=c('Age','Sex', 'EstimatedTotalIntraCranialVol'), resp.range=c(subcortical_column1:subcortical_column8),EP=ENIGMA_subcortical$SSD, data=data)
  RVI_SSD_mean <- (RVI_cortical_SSD$RVI$RVI + RVI_subcortical_SSD$RVI$RVI)/2
  
  # Checks correlation between RVI for schizophrenia and OCD in the sample
  cor(data$RVI_cortical, RVI_cortical_SSD$RVI$RVI) # Moderately similar, r=0.47
  cor(data$RVI_subcortical, RVI_subcortical_SSD$RVI$RVI) # Nearly identical, r=0.98
  cor(data$RVI_mean, RVI_SSD_mean) # Very similar, r=0.89
  
  summary(lmer(RVI_cortical ~ RVI_cortical_SSD$RVI$RVI + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE))
  summary(lmer(RVI_subcortical ~ RVI_subcortical_SSD$RVI$RVI + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE))
  summary(lmer(RVI_mean ~ RVI_SSD_mean + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE))
  
  summary(lmer(RVI_cortical_SSD$RVI$RVI ~ Group + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE))
  summary(lmer(RVI_subcortical_SSD$RVI$RVI ~ Group + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE))
  summary(lmer(RVI_SSD_mean ~ Group + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE))
  
  boxplot(data$RVI_mean ~ data$Group)
  mean(subset(data$RVI_mean, data$Group == 1))
  mean(subset(data$RVI_mean, data$Group == 0))
  
  data$Group <- factor(data$Group, labels = c("HC", "OCD"))
  
  data %>%
    ggplot(aes(x = Group, y = RVI_SSD_mean, fill = Group)) +
    ylim(-1, 1) +
    # add half-violin from ggdist package
    ggdist::stat_halfeye(
      ## custom bandwidth
      adjust = 0.5,
      ##move geom to the right
      justification = -.1,
      ## remove slab interval
      .width = 0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = .08
    )
  
  data %>%
    ggplot(aes(x = Group, y = RVI_subcortical_SSD$RVI$RVI, fill = Group)) +
    
    # add half-violin from ggdist package
    ggdist::stat_halfeye(
      ## custom bandwidth
      adjust = 0.5,
      ##move geom to the right
      justification = -.1,
      ## remove slab interval
      .width = 0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = .08
    )
  
  data %>%
    ggplot(aes(x = Group, y = RVI_cortical_SSD$RVI$RVI, fill = Group)) +
    
    # add half-violin from ggdist package
    ggdist::stat_halfeye(
      ## custom bandwidth
      adjust = 0.5,
      ##move geom to the right
      justification = -.1,
      ## remove slab interval
      .width = 0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = .08
    )
  
  lmer_RVI_mean_SSD <-lmer(ICA_Rerun7z_3 ~ Group + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE) 
  summary(lmer_RVI_mean_SSD)
  plot_model(lmer_RVI_mean_SSD, title = "Mean cortical and subcortical RVI", type = "std", show.values = TRUE)
  test <-lme.dscore(lmer_RVI_mean_SSD, data=data, type = "lme4")
  
  
  
## TEST
  library(emmeans)
  test.lmer <- lmer(lh_GWC_pericalcarine ~ Group + Sex + Age + SurfaceHoles + (1 | Site), data = data, REML = FALSE)
  summary(test.lmer)
  emmeans(test.lmer, specs = "Group")
  cor(residuals(test.lmer), residuals(test.lm))
  plot(residuals(test.lmer), residuals(test.lm))
  hist(residuals(test.lmer))
  hist(residuals(test.lm))
  
  cor(predict(test.lmer), predict(test.lm))
  plot(predict(test.lmer), predict(test.lm))
  
  
  test.lm <- lm(lh_GWC_pericalcarine ~ Group + Sex + Age + SurfaceHoles + as.factor(Site), data = data)
  summary(test.lm)
  emmeans(test.lm, specs = "Group")
  
  test.lm.noSite <- lm(lh_GWC_pericalcarine  ~ Group + Sex + Age + SurfaceHoles, data = data)
  summary(test.lm.noSite)
  emmeans(test.lm.noSite, specs = "Group")
  
  test <- lmer(lh_GWC_pericalcarine ~ Group + (1 | Site), data = data, REML = FALSE)
  emmeans(test, specs = "Group")
  
  boxplot(data$lh_GWC_pericalcarine ~ data$Group) # Overall median by group
  boxplot(data$lh_GWC_pericalcarine ~ data$Site) # Overall median by group
  boxplot(data$lh_GWC_pericalcarine ~ data$Group:data$Site) # Median by group by site
  
  test.lmer <- lmer(rh_GWC_pericalcarine ~ Group + Sex + Age + SurfaceHoles + (1 | Site), data = data, REML = FALSE)
  summary(test.lmer)
  
  emmeans(test.lmer, specs = "Group", lmer.df = "kenward-roger")
  emmeans(test.lmer, specs = "Group", lmer.df = "satterthwaite")

  estimated_mean <- predict(test.lmer)
  mean(subset(estimated_mean, data$Group == "HC"))
  mean(subset(estimated_mean, data$Group == "OCD"))
  

#
  
  cortical_column1 <- which(colnames(data) == "lh_GWC_bankssts")
  cortical_column34 <- cortical_column1+33
lh_df_bil_GWC <- data[cortical_column1:cortical_column34]  
data[cortical_column1]

# Using intraclass correlation to check absolute agreement between mean RVI for OCD and SSD
RVI_ICC_test<-data.frame(mean_RVI_OCD=matrix(0,848,1),
                         mean_RVI_SSD=matrix(0,848,1))
RVI_ICC_test$mean_RVI_OCD <- data$RVI_mean
RVI_ICC_test$mean_RVI_SSD <- RVI_SSD_mean
library(psych)
ICC(RVI_ICC_test, lmer = FALSE)


test_no_fixed <- lmer(ICA_Rerun7z_3 ~ (1 | Site), data = data, REML = FALSE)
icc(test_no_fixed)

test_fixed <- lmer(ICA_Rerun7z_3 ~ Group + Sex + Age + SurfaceHoles + (1 | Site), data = data, REML = FALSE)
icc(test_fixed)









# Vilde testing
# Aim: test GM and WM per IC3 and GWC
# Estimate Cohen's d for GM and WM per ROI
# Estimate Cohen's d for GM and WM drives GWC IC3?
# Test if contrast is driven by reduction in WM

ggplot_ROIs_spaced <- brain_regions(dk)[-4] # Get labels from ggseg and removes corpus callosum from list


ggplot_ROIs <- c("bankssts","caudalanteriorcingulate","caudalmiddlefrontal","cuneus",
                 "entorhinal","frontalpole","fusiform","inferiorparietal","inferiortemporal","insula",
                 "isthmuscingulate","lateraloccipital","lateralorbitofrontal","lingual",
                 "medialorbitofrontal","middletemporal","paracentral","parahippocampal",
                 "parsopercularis","parsorbitalis","parstriangularis","pericalcarine",
                 "postcentral","posteriorcingulate","precentral","precuneus",
                 "rostralanteriorcingulate","rostralmiddlefrontal",
                 "superiorfrontal","superiorparietal","superiortemporal",
                 "supramarginal","temporalpole","transversetemporal")



# Set up dataframe to hold results which are to be plotted by ggseg
results_Cohensd <- tibble(
  region = rep(c(ggplot_ROIs_spaced), 2),
  mean = as.double(NA),
  hemi = c(rep("left", 34), rep("right", 34))
)

# Estimate mean GWC per ROI per group
for(hemi in c("lh","rh")){ # Loop over hemispheres
  i=0 # Sets starting point for counter
  if (hemi == "lh"){ 
    i=i+0 # Set to zero as left is first in the dataframe
  } else if (hemi == "rh"){
    i=i+34
  }
  
  for(region in ggplot_ROIs){ # Loop over regions
    
    i=i+1
    
    data$r <- data[[paste0(hemi,"_GWC_",region)]]
    lmer.test <- lmer(r ~ ICA_Rerun7z_3 + Sex + Age + SurfaceHoles + (1 | Site), data = data, REML = FALSE) 
    Cohensd <- lme.dscore(lmer.test, data=data, type = "lme4")
    
    # print(i)
    # print(Cohensd$d[1])
    results_Cohensd[i,"mean"] <- Cohensd$d[1]
  }
}

write.xlsx(results_Cohensd_wm, file = "GWC_group_CohensD_results.xlsx")


# WM per IC3
ggplot_ROIs_spaced <- brain_regions(dk)[-4] # Get labels from ggseg and removes corpus callosum from list


ggplot_ROIs <- c("bankssts","caudalanteriorcingulate","caudalmiddlefrontal","cuneus",
                 "entorhinal","frontalpole","fusiform","inferiorparietal","inferiortemporal","insula",
                 "isthmuscingulate","lateraloccipital","lateralorbitofrontal","lingual",
                 "medialorbitofrontal","middletemporal","paracentral","parahippocampal",
                 "parsopercularis","parsorbitalis","parstriangularis","pericalcarine",
                 "postcentral","posteriorcingulate","precentral","precuneus",
                 "rostralanteriorcingulate","rostralmiddlefrontal",
                 "superiorfrontal","superiorparietal","superiortemporal",
                 "supramarginal","temporalpole","transversetemporal")
results_Cohensd <- tibble(
  region = rep(c(ggplot_ROIs_spaced), 2),
  mean = as.double(NA),
  hemi = c(rep("left", 34), rep("right", 34))
)

for(hemi in c("lh","rh")){ # Loop over hemispheres
  i=0 # Sets starting point for counter
  if (hemi == "lh"){ 
    i=i+0 # Set to zero as left is first in the dataframe
  } else if (hemi == "rh"){
    i=i+34
  }
  
  for(region in ggplot_ROIs){ # Loop over regions
    
    i=i+1
    
    data$wm <- data[[paste0(hemi,"_WM_",region)]]
    lmer.test_GWC <- lmer(wm ~ ICA_Rerun7z_3 + Sex + Age + SurfaceHoles + (1 | Site), data = data, REML = FALSE) 
    summary <- summary(lmer.test_GWC)
    Cohensd <- lme.dscore(lmer.test, data=data, type = "lme4")

    # print(i)
    # print(Cohensd$d[1])
    results_Cohensd[i,"mean"] <- Cohensd$d[1]
  }
}

write.xlsx(results_Cohensd, file = "WM_IC3_CohensD_results.xlsx")
coefs <- data.frame(coef(summary(lmer.test)))
coefs$p.z <- 2*(1 -pnorm(abs(coefs$t.value)))
coefs$p.Satt <- coef(summary(lmer.test_GWC)) [, 5]
print(coefs)





# GM per IC3
ggplot_ROIs_spaced <- brain_regions(dk)[-4] # Get labels from ggseg and removes corpus callosum from list


ggplot_ROIs <- c("bankssts","caudalanteriorcingulate","caudalmiddlefrontal","cuneus",
                 "entorhinal","frontalpole","fusiform","inferiorparietal","inferiortemporal","insula",
                 "isthmuscingulate","lateraloccipital","lateralorbitofrontal","lingual",
                 "medialorbitofrontal","middletemporal","paracentral","parahippocampal",
                 "parsopercularis","parsorbitalis","parstriangularis","pericalcarine",
                 "postcentral","posteriorcingulate","precentral","precuneus",
                 "rostralanteriorcingulate","rostralmiddlefrontal",
                 "superiorfrontal","superiorparietal","superiortemporal",
                 "supramarginal","temporalpole","transversetemporal")

results_Cohensd <- tibble(
  region = rep(c(ggplot_ROIs_spaced), 2),
  mean = as.double(NA),
  hemi = c(rep("left", 34), rep("right", 34))
)

for(hemi in c("lh","rh")){ # Loop over hemispheres
  i=0 # Sets starting point for counter
  if (hemi == "lh"){ 
    i=i+0 # Set to zero as left is first in the dataframe
  } else if (hemi == "rh"){
    i=i+34
  }
  
  for(region in ggplot_ROIs){ # Loop over regions
    
    i=i+1
    
    data$gm <- data[[paste0(hemi,"_GM_",region)]]
    lmer.test <- lmer(gm ~ ICA_Rerun7z_3 + Sex + Age + SurfaceHoles + (1 | Site), data = data, REML = FALSE) 
    Cohensd <- lme.dscore(lmer.test, data=data, type = "lme4")
    
    # print(i)
    # print(Cohensd$d[1])
    results_Cohensd[i,"mean"] <- Cohensd$d[1]
  }
}

write.xlsx(results_Cohensd, file = "GM_IC3_CohensD_results.xlsx")



# WM PER group
# Set up dataframe to hold results which are to be plotted by ggseg
results_Cohensd <- tibble(
  region = rep(c(ggplot_ROIs_spaced), 2),
  mean = as.double(NA),
  hemi = c(rep("left", 34), rep("right", 34))
)


for(hemi in c("lh","rh")){ # Loop over hemispheres
  i=0 # Sets starting point for counter
  if (hemi == "lh"){ 
    i=i+0 # Set to zero as left is first in the dataframe
  } else if (hemi == "rh"){
    i=i+34
  }
  
  for(region in ggplot_ROIs){ # Loop over regions
    
    i=i+1
    
    data$wm <- data[[paste0(hemi,"_WM_",region)]]
    lmer.test <- lmer(wm ~ Group + Sex + Age + SurfaceHoles + (1 | Site), data = data, REML = FALSE) 
    Cohensd <- lme.dscore(lmer.test, data=data, type = "lme4")
    
    # print(i)
    # print(Cohensd$d[1])
    results_Cohensd[i,"mean"] <- Cohensd$d[1]
    
  }
}

write.xlsx(results_Cohensd, file = "WM_group_CohensD_results_test.xlsx")

results_Cohensd %>%
  ggplot() +
  geom_brain(atlas = dk,
             position=position_brain(hemi ~ side),
             colour ="black",
             aes(fill = mean)) +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  labs(fill="Cohens's d relative to HC") +
  ggtitle("Regional WM in OCD versus HC") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())


lmer.test <- lmer(ICA_Rerun7z_3 ~ Group + Sex + Age + SurfaceHoles + (1 | Site), data = data, REML = FALSE)




# GM per group

ggplot_ROIs_spaced <- brain_regions(dk)[-4] # Get labels from ggseg and removes corpus callosum from list


ggplot_ROIs <- c("bankssts","caudalanteriorcingulate","caudalmiddlefrontal","cuneus",
                    "entorhinal","frontalpole","fusiform","inferiorparietal","inferiortemporal","insula",
                    "isthmuscingulate","lateraloccipital","lateralorbitofrontal","lingual",
                    "medialorbitofrontal","middletemporal","paracentral","parahippocampal",
                    "parsopercularis","parsorbitalis","parstriangularis","pericalcarine",
                    "postcentral","posteriorcingulate","precentral","precuneus",
                    "rostralanteriorcingulate","rostralmiddlefrontal",
                    "superiorfrontal","superiorparietal","superiortemporal",
                    "supramarginal","temporalpole","transversetemporal")



# Set up dataframe to hold results which are to be plotted by ggseg
results_Cohensd <- tibble(
  region = rep(c(ggplot_ROIs_spaced), 2),
  mean = as.double(NA),
  hemi = c(rep("left", 34), rep("right", 34))
)


for(hemi in c("lh","rh")){ # Loop over hemispheres
  i=0 # Sets starting point for counter
  if (hemi == "lh"){ 
    i=i+0 # Set to zero as left is first in the dataframe
  } else if (hemi == "rh"){
    i=i+34
  }
  
  for(region in ggplot_ROIs){ # Loop over regions
    
    i=i+1
    
    data$gm <- data[[paste0(hemi,"_GM_",region)]]
    lmer.test <- lmer(gm ~ Group + Sex + Age + SurfaceHoles + (1 | Site), data = data, REML = FALSE) 
    Cohensd <- lme.dscore(lmer.test, data=data, type = "lme4")
    
    # print(i)
    # print(Cohensd$d[1])
    results_Cohensd[i,"mean"] <- Cohensd$d[1]
    
  }
}

write.xlsx(results_Cohensd, file = "GM_group_CohensD_results.xlsx")



# CT

ggplot_ROIs_spaced <- brain_regions(dk)[-4] # Get labels from ggseg and removes corpus callosum from list


ggplot_ROIs <- c("bankssts","caudalanteriorcingulate","caudalmiddlefrontal","cuneus",
                 "entorhinal","frontalpole","fusiform","inferiorparietal","inferiortemporal","insula",
                 "isthmuscingulate","lateraloccipital","lateralorbitofrontal","lingual",
                 "medialorbitofrontal","middletemporal","paracentral","parahippocampal",
                 "parsopercularis","parsorbitalis","parstriangularis","pericalcarine",
                 "postcentral","posteriorcingulate","precentral","precuneus",
                 "rostralanteriorcingulate","rostralmiddlefrontal",
                 "superiorfrontal","superiorparietal","superiortemporal",
                 "supramarginal","temporalpole","transversetemporal")



# Set up dataframe to hold results which are to be plotted by ggseg
results_Cohensd <- tibble(
  region = rep(c(ggplot_ROIs_spaced), 2),
  mean = as.double(NA),
  hemi = c(rep("left", 34), rep("right", 34))
)


for(hemi in c("lh","rh")){ # Loop over hemispheres
  i=0 # Sets starting point for counter
  if (hemi == "lh"){ 
    i=i+0 # Set to zero as left is first in the dataframe
  } else if (hemi == "rh"){
    i=i+34
  }
  
  for(region in ggplot_ROIs){ # Loop over regions
    
    i=i+1
    
    data$gm <- data[[paste0(hemi,"_thickness_",region)]]
    lmer.test <- lmer(gm ~ Group + Sex + Age + SurfaceHoles + (1 | Site), data = data, REML = FALSE) 
    Cohensd <- lme.dscore(lmer.test, data=data, type = "lme4")
    
    # print(i)
    # print(Cohensd$d[1])
    results_Cohensd[i,"mean"] <- Cohensd$d[1]
    
  }
}

write.xlsx(results_Cohensd, file = "Thickness_group_CohensD_results.xlsx")




###
#Testing why RVI gives negative values
data_test <- read_sav(file="S:/Project/OBIC_dataset/Backup/Cortical_myelination_FINAL_BACKUP_TEST_07Aug22.sav")

cor(data$lh_thickness_bankssts, data_test$lh_thickness_bankssts)

### Calculating GWC based on ROI signal intensities
data$manual_rh_pericalcarine_GWC <- 100*(data$rh_WM_pericalcarine-data$rh_GM_pericalcarine)/((data$rh_WM_pericalcarine+data$rh_GM_pericalcarine)/2)

