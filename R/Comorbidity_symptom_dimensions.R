# Written by Anders Lillevik Thorsen, May 2022
# This script compares GWC by comorbidity and symptom dimensions in OCD patients in the OBIC dataset

library(tidyverse)
library(lme4)
library(skimr)
library(ggplot2)
library(haven)
library(effectsize)
library(parameters)
library(corrplot)
library(ggseg)
library(psych)
library(lmerTest)

# Setup environment
rm(list=ls()) # Clears variables
options(scipen = 999) # Gives decimals rather than power for large number
setwd("S:/Project/OBIC/R")

# Load data
#data <- read_sav(file="S:/Project/OBIC_dataset/Cortical_myelination_FINAL_08May22.sav")
data <- read_sav(file="S:/Project/OBIC/dataset/Cortical_myelination_FINAL_20July22.sav")

data$Age_sq = data$Age*data$Age

# Plot GWC using ggseg


# Compare OCD patients with and without anxiety comorbidity on ICA_3
data_OCD <- subset(data, Group == 1 & GWC_anxiety < 2) # Select OCD patients only

hist(data_OCD$Group)
hist(data$GWC_anxiety)

anxiety_comorb <- lmer(ICA_Rerun7z_3 ~ GWC_anxiety + Sex + Age + SurfaceHoles + (1 | Site), data = data_OCD, REML = FALSE) 
summary(anxiety_comorb)

anxiety_comorb <- lmer(ICA_Rerun7z_1 ~ GWC_anxiety + Sex + Age + SurfaceHoles + (1 | Site), data = data_OCD, REML = FALSE) 
summary(anxiety_comorb)

# With Pa
test <- subset(data, GWC_anxiety != 999)
hist(test$Group)
hist(test$GWC_anxiety)

anxiety_comorb <- lmer(ICA_Rerun7z_3 ~ as.factor(GWC_anxiety) + Sex + Age + SurfaceHoles + (1 | Site), data = test, REML = FALSE) 
summary(anxiety_comorb)

# Compare OCD patients with and without depressive comorbidity on ICA_3
data_OCD <- subset(data, Group == 1 & GWC_depressive < 2) # Select OCD patients only

depressive_comorb <- lmer(ICA_Rerun7z_3 ~ GWC_depressive + Sex + Age + SurfaceHoles + (1 | Site), data = data_OCD, REML = FALSE) 
summary(depressive_comorb)

depressive_comorb <- lmer(ICA_Rerun7z_4 ~ GWC_depressive + Sex + Age + SurfaceHoles + (1 | Site), data = data_OCD, REML = FALSE) 
summary(depressive_comorb)

# Compare OCD patients with different symptom dimensions on ICA_3
data_OCD <- subset(data, Group == 1) # Select OCD patients only

hist(data_OCD$Agr_Check)
hist(data$Contam_Clean)
hist(data$Sym_Ordering)
hist(data$Sex_Rel)
hist(data$Hoarding)

symptomdim_comorb <- lmer(ICA_Rerun7z_3 ~ Agr_Check + Contam_Clean + Sym_Ordering + Sex_Rel + Hoarding + Sex + Age + SurfaceHoles + (1 | Site), data = data_OCD, REML = FALSE) 
summary(symptomdim_comorb)
plot_model(symptomdim_comorb, title = "Component 3", type = "std", show.values = TRUE)
effectsize(symptomdim_comorb)
lme.dscore(symptomdim_comorb, data=data, type = "lme4")

library(sjPlot)
tab_model(symptomdim_comorb)

symptomdim_comorb <- lmer(ICA_Rerun7z_7 ~ Agr_Check + Contam_Clean + Sym_Ordering + Sex_Rel + Hoarding + Sex + Age + SurfaceHoles + (1 | Site), data = data_OCD, REML = FALSE) 
summary(symptomdim_comorb)

library(performance)
check_collinearity(symptomdim_comorb) # Multicolinearity of symptom dimensions not a problem

# Compare medicated and unmedicated OCD 
data_medOCD_unnmedOCD <- subset(data, Medicated < 2) # Select OCD patients only

hist(data_medOCD_unnmedOCD$Medicated)

med_unmed <- lmer(ICA_Rerun7z_3 ~ Medicated + Sex + Age + SurfaceHoles + (1 | Site), data = data_medOCD_unnmedOCD, REML = FALSE) 
summary(med_unmed)

med_unmed <- lmer(ICA_Rerun7z_4 ~ Medicated + Sex + Age + SurfaceHoles + (1 | Site), data = data_medOCD_unnmedOCD, REML = FALSE) 
summary(med_unmed)

# Compared unmedicated OCD and HC
data_unmedOCD_HC <- subset(data, Medicated != 1) # Select OCD patients only

hist(data_unmedOCD_HC$Medicated)

med_unmed <- lmer(ICA_Rerun7z_3 ~ Medicated + Sex + Age + SurfaceHoles + (1 | Site), data = data_unmedOCD_HC, REML = FALSE) 
summary(med_unmed)

med_unmed <- lmer(ICA_Rerun7z_7 ~ Medicated + Sex + Age + SurfaceHoles + (1 | Site), data = data_unmedOCD_HC, REML = FALSE) 
summary(med_unmed)

# Compared unmedicated OCD and HC
data_medOCD_HC <- subset(data, Medicated != 0) # Select OCD patients only

hist(data_medOCD_HC$Medicated)

med_unmed <- lmer(ICA_Rerun7z_3 ~ Medicated + Sex + Age + SurfaceHoles + (1 | Site), data = data_medOCD_HC, REML = FALSE) 
summary(med_unmed)

med_unmed <- lmer(ICA_Rerun7z_7 ~ Medicated + Sex + Age + SurfaceHoles + (1 | Site), data = data_medOCD_HC, REML = FALSE) 
summary(med_unmed)

# Compared childhod and adult onset OCD
data_onset <- subset(data, Onset_OCD <= 1 ) # Select OCD patients only

hist(data_onset$Onset_OCD)

onset_OCD <- lmer(ICA_Rerun7z_3 ~ Onset_OCD + Sex + Age + SurfaceHoles + (1 | Site), data = data_onset, REML = FALSE) 
summary(onset_OCD)

onset_OCD <- lmer(ICA_Rerun7z_7 ~ Onset_OCD + Sex + Age + SurfaceHoles + (1 | Site), data = data_onset, REML = FALSE) 
summary(onset_OCD)

# Relate RVI to ICA components in OCD
data_OCD <- subset(data, data$Group == "OCD") # Select OCD patients only

# Relate RVI to ICA_Rerun7z_3 with group
lmer_OCD_RVI_mean <- lmer(ICA_Rerun7z_3 ~ RVI_mean + Group + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE) 
summary(lmer_OCD_RVI_mean)
plot_model(lmer_OCD_RVI_mean, title = "Component 3", type = "std", show.values = TRUE)
lme.dscore(lmer_OCD_RVI_mean, data=data, type = "lme4")

lmer_OCD_RVI_mean <- lmer(ICA_Rerun7z_7 ~ RVI_mean + Group + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE) 
summary(lmer_OCD_RVI_mean)

# Within OCD patients only
lmer_OCD_RVI_cort <-lmer(ICA_Rerun7z_7 ~ RVI_cortical + Sex + SurfaceHoles + Age + (1 | Site), data = data_OCD, REML = FALSE) 
summary(lmer_OCD_RVI_cort)

lmer_OCD_RVI_subcort <-lmer(ICA_Rerun7z_7 ~ RVI_subcortical + Sex + SurfaceHoles + Age + (1 | Site), data = data_OCD, REML = FALSE) 
summary(lmer_OCD_RVI_subcort)

lmer_OCD_RVI_mean <- lmer(ICA_Rerun7z_3 ~ RVI_mean + Sex + SurfaceHoles + Age + (1 | Site), data = data_OCD, REML = FALSE) 
summary(lmer_OCD_RVI_mean)

lmer_YBOCS_RVI_mean <- lmer(RVI_mean ~ YBOCS_total + Sex + SurfaceHoles + Age + (1 | Site), data = data_OCD, REML = FALSE) 
summary(lmer_YBOCS_RVI_mean)

symptomdim_RVI <- lmer(RVI_mean ~ Agr_Check + Contam_Clean + Sym_Ordering + Sex_Rel + Hoarding + Sex + Age + SurfaceHoles + (1 | Site), data = data_OCD, REML = FALSE) 
summary(symptomdim_RVI)

med_unmed_RVI <- lmer(RVI_mean ~ Medicated + Sex + Age + SurfaceHoles + (1 | Site), data = data_OCD, REML = FALSE) 
summary(med_unmed_RVI)

anxiety_RVI <- lmer(RVI_mean ~ GWC_anxiety + Sex + Age + SurfaceHoles + (1 | Site), data = data_OCD, REML = FALSE) 
summary(anxiety_RVI)

depressive_RVI <- lmer(RVI_mean ~ GWC_depressive + Sex + Age + SurfaceHoles + (1 | Site), data = data_OCD, REML = FALSE) 
summary(depressive_RVI)


data_OCD_depressive <- subset(data, Group == 1 & GWC_depressive < 2) # Select OCD patients only
depressive_RVI <- lmer(RVI_mean ~ GWC_depressive + Sex + Age + SurfaceHoles + (1 | Site), data = data_OCD_depressive, REML = FALSE) 
summary(depressive_RVI)

onset_OCD <- lmer(RVI_mean ~ Onset_OCD + Sex + Age + SurfaceHoles + (1 | Site), data = data_OCD, REML = FALSE) 
summary(onset_OCD)

data_OCD <- subset(data, Group == 1 & GWC_depressive < 2 & GWC_anxiety < 2) # Select OCD patients only
everything_OCD_RVI <- lmer(RVI_mean ~ Medicated + YBOCS_total + GWC_anxiety + GWC_depressive+ Sex + Age + SurfaceHoles + (1 | Site), data = data_OCD, REML = FALSE) 
summary(everything_OCD_RVI)

summary(lmer(RVI_mean ~ Sex + SurfaceHoles + Age + (1 | Site), data = data_OCD, REML = FALSE)) 


# Testing here
cor(data$RVI_subcortical, data$RVI_mean)
cor(data$RVI_cortical, data$RVI_mean)
cor(data$RVI_subcortical, data$RVI_cortical)

cor(data$RVI_cortical, data$lh_GWC_lateraloccipital)

  #### Status 12 May
# Anxiety related to Comp1
# Depressive related to Comp4, 5, 7
# Multicolinearity of symptom dimensions not a problem
# Sex_rel related to Comp3 (!)
# Contam_Clean related to Comp2, 4, 5, 7
# Hoarding related to Comp7
# Med and unmed OCD differs in Comp1, 4
# Unmed OCD and HC never differs
# MedOCD and HC differs in Comp3 (!) and Comp1
# Onset related to Comp1, 2, 4, 5
# RVI has been coded according to Boedhoe, 2020 (largest ENIGMA-OCD study of cortical thickness and subcortical volume) to coded and used
# RVI estimation run and RVI is significantly higher in OCD than HC. However, mean RVI is very close to 0 for both OCD and HC and Cohen's d is 0.17-0.20
# Mean and subcortical RVI is related to Comp3 (and others) in OCD patients
# RVI is not related to any clinical characteristic
# Needs to confirm RVI estimation and investigate subcortal, cortical and mean RVI to GWC

# Load ENIGMA group differences from RVIpkg, add OCD data from Boedhoe, 2020 Am. J. Psychiatry
ENIGMA_subcortical <- RVIpkg::EP.Subcortical
ENIGMA_cortical <- RVIpkg::EP.GM

ENIGMA_cortical[nrow(ENIGMA_cortical) + 1,1] = "temporalpole" # Add temporal pole as this is missing IRVpkg
#ENIGMA_subcortical<- ENIGMA_subcortical[-c(1), ] # Drop Ventricle as this is missing in OCD

ENIGMA_subcortical$OCD <- c(.11, -0.05, 0.01, 0, 0.09, -0.09, -0.06, -0.03) # Missing Ventricle (row 1)
ENIGMA_cortical$OCD <- c(-0.0275,-0.0175,-0.0845,-0.0375,-0.024,-0.099,-0.14,-0.073,-0.06,-0.073,-0.102,-0.048,-0.093,-0.0955,-0.002,-0.0695,-0.0725,-0.0595,-0.044,0.018,0.0065,-0.064,-0.03,-0.099,-0.0345,-0.093,-0.0535,-0.0545,0.0035,-0.0245,-0.0135,-0.014,-0.065,0.0245)

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

# Test to check if corrected columns are selected for cortical thickness/subcortical volume
cort_cols <- colnames(data[0,401:434]) 
if(cort_cols[1] != "bil_thickness_bankssts" & cort_cols[34] != "bil_thickness_insula"){
  print("ERROR - you have selected the wrong columns")
  } else {
  print("ALL GOOD - you have selected the correct columns")
  }

subcort_cols <- colnames(data[0,437:444])
if(subcort_cols[1] != "bil_vol_lateralventricle" & subcort_cols[8] != "bil_vol_Accumbensarea"){
  print("ERROR - you have selected the wrong columns")
  } else {
    print("ALL GOOD - you have selected the correct columns")
  }

# Run IRVpkg for cortical thickness
library(RVIpkg)
RVI_cortical <- RVI_func(ID='BIDS', DXcontrol='Group==0', covariates=c('Age','Sex'), resp.range=c(401:434),EP=ENIGMA_cortical$OCD, data=data)
data$RVI_cortical <- RVI_cortical$RVI$RVI

# Run IRVpkg for subcortical volume
RVI_subcortical <- RVI_func(ID='BIDS', DXcontrol='Group==0', covariates=c('Age','Sex', 'EstimatedTotalIntraCranialVol'), resp.range=c(437:444),EP=ENIGMA_subcortical$OCD, data=data)
data$RVI_subcortical <- RVI_subcortical$RVI$RVI

# Calculate mean RVI across the cortex and subcortex
data$RVI_mean <- (data$RVI_cortical+data$RVI_subcortical)/2


describeBy(data$RVI_subcortical, group = data$Group) # M HC = 0.01, M OCD = 0.08
describeBy(data$RVI_cortical, group = data$Group) # M HC = 0, M OCD = 0.03
describeBy(data$RVI_mean, group = data$Group) # M HC = 0, M OCD = 0.05

# Test Group difference in RVI
lmer_RVI_cort <-lmer(RVI_cortical ~ Group + Sex + SurfaceHoles + Age +(1 | Site), data = data, REML = FALSE) 
summary(lmer_RVI_cort)
plot_model(lmer_RVI_cort, title = "Cortical RVI", type = "std", show.values = TRUE)
lme.dscore(lmer_RVI_cort, data=data, type = "lme4")

lmer_RVI_subcort <-lmer(RVI_subcortical ~ Group + Sex + SurfaceHoles + Age +(1 | Site), data = data, REML = FALSE) 
summary(lmer_RVI_subcort)
plot_model(lmer_RVI_subcort, title = "Subcortical RVI", type = "std", show.values = TRUE)
lme.dscore(lmer_RVI_subcort, data=data, type = "lme4")

lmer_RVI_mean <-lmer(RVI_mean ~ Group + Sex + SurfaceHoles + Age +(1 | Site), data = data, REML = FALSE) 
summary(lmer_RVI_mean)
plot_model(lmer_RVI_mean, title = "Mean cortical and subcortical RVI", type = "std", show.values = TRUE)
lme.dscore(lmer_RVI_mean, data=data, type = "lme4")

lmer_RVI_mean <-lmer(ICA_Rerun7z_3 ~ RVI_mean + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE) 
summary(lmer_RVI_mean)

lmer_RVI_mean <-lmer(ICA_Rerun7z_3 ~ Group + RVI_mean + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE) 
summary(lmer_RVI_mean)

boxplot(subset(data$RVI_mean, data$Group == 1))
boxplot(subset(data$RVI_mean, data$Group == 0))
boxplot(data$RVI_mean ~ data$Group)



# Calculate mean GWC across cortex
library(fame)


rowMeans(data)

data$lh_mean_GWC <- rowMeans(data[, 118:151])
data$rh_mean_GWC <- rowMeans(data[, 153:186])
data$bil_mean_GWC <- rowMeans(data[, 449:450])

lmer_RVI_cort <-lmer(bil_mean_GWC~ RVI_cortical*Group+ Sex + SurfaceHoles + Age +(1 | Site), data = data, REML = FALSE) 
summary(lmer_RVI_cort)

lmer_mean_GWC <-lmer(bil_mean_GWC ~ Group+ Sex + SurfaceHoles + Age +(1 | Site), data = data, REML = FALSE) 
summary(lmer_mean_GWC)

hist(data$RVI_mean)

library(ggplot2)
RVI_cort_plot <- ggplot(data, aes(x = as.factor(Group), y = bil_mean_GWC)) +
  geom_violin()
RVI_cort_plot
RVI_cort_plot + geom_boxplot(width = 0.3)

RVI_mean_plot <- ggplot(data, aes(x = as.factor(Group), y = RVI_mean)) +
  geom_violin()
RVI_mean_plot
RVI_mean_plot + geom_boxplot(width = 0.3)

library(psych)
describeBy(data$RVI_mean, group = data$Group)
a <- cohen.d(data$ICA_Rerun7z_3, group=as.factor(data$Group))
a$cohen.d

describeBy(data$ICA_Rerun7z_3, group = data$Group)
a <- cohen.d(data$ICA_Rerun7z_3, group=as.factor(data$Group))
a$cohen.d
hist(data$ICA_Rerun7z_3)

a <- cohen.d(data$RVI_mean, group=as.factor(data$Group))
a$cohen.d

# ICV = data$EstimatedTotalIntraCranialVol
which( colnames(data) == "lh_GWC_bankssts")


# Delete below here
data$RVI_BIDS <- RVI_cortical$RVI$BIDS
data$BIDS[603]
data$RVI_BIDS[603]

ROI <- data[0,356:360]
ROI <- colnames(ROI)




# Testing below here
ICA_z_1 <- lmer(ICA_Rerun7z_1 ~ Group + Sex + (1 | Site) + Age + Age_sq, data = OBIC, REML = FALSE) 
summary(ICA_z_1)

ICA_z_2 <- lmer(ICA_Rerun7z_2 ~ Group + Sex + (1 | Site) + Age + Age_sq, data = OBIC, REML = FALSE) 
summary(ICA_z_2)

ICA_z_3 <- lmer(ICA_Rerun7z_3 ~ Group + Sex + SurfaceHoles + Age +  (1 | Site), data = data, REML = FALSE) 
summary(ICA_z_3)
plot_model(ICA_z_3, title = "Component 3", type = "std", show.values = TRUE)

library(EMAtools)
lme.dscore(ICA_z_3, data=data, type = "lme4")

ICA_z_3_Age2 <- lmer(ICA_Rerun7z_3 ~ Group + Sex + SurfaceHoles + (1 | Site) + Age + Age_sq, data = OBIC, REML = FALSE) 
summary(ICA_z_3_Age2)

anova(ICA_z_3, ICA_z_3_Age2)

ICA_5 <- lmer(ICA7_5 ~ Group + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE) 
summary(ICA_5)

ICA_z_4 <- lmer(ICA_Rerun7z_4 ~ Group + Sex + (1 | Site) + Age + Age_sq, data = OBIC, REML = FALSE) 
summary(ICA_z_4)

ICA_z_5 <- lmer(ICA_Rerun7z_5 ~ Group + Sex + (1 | Site) + Age + Age_sq, data = OBIC, REML = FALSE) 
summary(ICA_z_5)

ICA_z_6 <- lmer(ICA_Rerun7z_6 ~ Group + Sex + (1 | Site) + Age + Age_sq, data = OBIC, REML = FALSE) 
summary(ICA_z_6)

ICA_z_7 <- lmer(ICA_Rerun7z_7 ~ Group + Sex + (1 | Site) + Age + Age_sq, data = OBIC, REML = FALSE) 
summary(ICA_z_7)

cor(data$ICA7_5, data$ICA_Rerun7z_3, method = "kendall")

p <- c(0.3237,0.38410,0.00837,0.727035,0.212325,0.976548,0.373037)

fdrp <- p.adjust(p, method="fdr")

data_count_Bergen <- subset(data, Site == 9 & Group == 1)
print(data_count_Bergen$BIDS)

#Test group differences for the seven components"
library(tidyverse)

#data <- data
data$r<-0
ROI.names_LME_OBIC<-data %>%
  select(contains("*ICA"), (contains("_Rerun")), -(contains("Medicated")))
region_OBIC_LME<-names(ROI.names_LME_OBIC)

stats_LME_ICA<-data.frame(t=matrix(0,7,1), p=matrix(0,7,1))
row.names(stats_LME_ICA)<-region_OBIC_LME

for(region in region_OBIC_LME){
  #add the ICA value 
  #data$r <- ROI.names_LME_OBIC[,region]
  data$r<- data[[paste0(ROI.names_LME_OBIC[,region])]
  print(region)

  
  #run the LMER
  ICA.test <- lmer(r ~ Group + Sex + (1 | Site) + Age + Age_sq, data = data, REML = FALSE) 
  s=summary(ICA.test)
  
  #Add values to table
  stats_LME_ICA[region,"t"]<-s$coefficients[2,4]
  stats_LME_ICA[region,"p"]<-s$coefficients[2,5]
  
}


# altered code

for(region in 1:length(region_OBIC_LME)){
  #add the ICA value 
  data$r <- data[region_OBIC_LME[region]]
                print(region)
                
                
                #run the LMER
                ICA.test <- lmer(r ~ Group + Sex + (1 | Site) + Age + Age_sq, data = data, REML = FALSE) 
                s=summary(ICA.test)
                
                #Add values to table
                stats_LME_ICA[region,"t"]<-s$coefficients[2,4]
                stats_LME_ICA[region,"p"]<-s$coefficients[2,5]
                
}


# This works by 12 May
for(region in 1:length(region_OBIC_LME)){
  #add the ICA value 
  #data$r <- data[region_OBIC_LME[region]]
  
  data$r <- data[[paste0("ICA_Rerun7z_", region)]]
  
  #run the LMER
  ICA.test <- lmer(r ~ Group + Sex + Age + SurfaceHoles + (1 | Site), data = data, REML = FALSE) 
  s=summary(ICA.test)
  
  #Add values to table
  stats_LME_ICA[region,"t"]<-s$coefficients[2,4]
  stats_LME_ICA[region,"p"]<-s$coefficients[2,5]

}

# do fdr correction
p_fdr_LME<-p.adjust(stats_LME_ICA$p, method="fdr")
stats_LME_ICA$p_fdr<-p_fdr_LME


# Test if new and old datasets match
Cortical_myelination_with_ICA_ROI_Bruk_denne <- read_sav("Anders/Cortical_myelination_with_ICA_ROI_Bruk_denne.sav")

Cortical_myelination_with_ICA_ROI_Bruk_denne <- Cortical_myelination_with_ICA_ROI_Bruk_denne[order(Cortical_myelination_with_ICA_ROI_Bruk_denne$BIDS),]
cor(data$Group, Cortical_myelination_with_ICA_ROI_Bruk_denne$Group)
cor(data$ICA_Rerun7z_3, Cortical_myelination_with_ICA_ROI_Bruk_denne$ICA_Rerun7z_3)

# Dataet matcher bortsett fra Group (r=0.9902)
# Finn ut hvilke caser som har skiftet Group

for(n in 1:nrow(data)){
  if(data$Group[n] == Cortical_myelination_with_ICA_ROI_Bruk_denne$Group[n]){
    data$same_group[n] = 1
    Cortical_myelination_with_ICA_ROI_Bruk_denne$same_group[n] = 1
  } else {
      data$same_group[n] = 0
      Cortical_myelination_with_ICA_ROI_Bruk_denne$same_group[n] = 0
      }
}

data$same_group
Cortical_myelination_with_ICA_ROI_Bruk_denne$same_group

subjects <- print(subset(data$BIDS, data$same_group == 0))
# Returns "sub-09subject00066_T1w" "sub-09subject00067_T1w" "sub-09subject00068_T1w" "sub-09subject00069_T1w"


Cortical_myelination_with_ICA_ROI_Bruk_denne$Group



data$lh_GWC_bankssts[1]
Cortical_myelination_with_ICA_ROI_Bruk_denne$lh_GWC_bankssts[1]

ICA_z_3 <- lmer(ICA_Rerun7z_3 ~ Group + Sex + (1 | site) + Age, data = Cortical_myelination_with_ICA_ROI_Bruk_denne, REML = FALSE) 
summary(ICA_z_3)


# Count cases
data_OCD <- subset(data, data$Group == 1)
summary(as.factor(data_OCD$Onset_OCD))
hist(data$GWC_anxiety)

group_by(data$Age, group = data$Group) # M HC = 0.01, M OCD = 0.08

grouped_data <- data %>% group_by(Group)
summary(grouped_data)

hist(data$Age)
qqnorm(y = data$Age)
describeBy(data$Age, group = data$Group)
summary(data$Age)
summary(lmer(Age ~ Group + (1 | Site), data = data, REML = FALSE))

describeBy(data$Education, group = data$Group)
summary(data$Education)
hist(data$Education)
qqnorm(y = data$Education); qqline(y = data$Education, col = 10)
summary(lmer(Education ~ Group + Age + Sex + (1 | Site), data = data, REML = FALSE))

data %>% count(Sex, Group)

summary()

summary(glmer(Sex ~ Group + (1 | Site), data = data, family = binomial))

summary(glmer(Group ~ Sex + (1 | Site), data = data, family = binomial))

summary(as.factor(data$Group))
data %>% count(Medicated, Group)
data %>% count(Agr_Check, Group)
data %>% count(Contam_Clean, Group)
data %>% count(Sym_Ordering, Group)
data %>% count(Sex_Rel, Group)
data %>% count(Hoarding, Group)


