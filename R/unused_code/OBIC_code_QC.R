# Written by Anders Lillevik Thorsen, August 2022
# This script will read and merge csv of different types to the main OBIC dataset
# Load relevant packages

# Setup environment
rm(list=ls()) # Clears variables
options(scipen = 999) # Gives decimals rather than power for large number

setwd("S:/Project/OBIC_dataset/Anders")
#setwd("S:/Project/OBIC/R")

# Load data
data <- read_sav(file="S:/Project/OBIC_dataset/Cortical_myelination_FINAL_RVI_07Aug22.sav")

ICA_weights <- read.csv(file="S:/Project/OBIC_dataset/Anders/ICA_weights_ICA.csv", sep = ";", dec = ",")

ICA_weights$ICA_1 <- as.numeric(ICA_weights$ICA_1)
ICA_weights$ICA_2 <- as.numeric(ICA_weights$ICA_2)
ICA_weights$ICA_3 <- as.numeric(ICA_weights$ICA_3)
ICA_weights$ICA_4 <- as.numeric(ICA_weights$ICA_4)
ICA_weights$ICA_5 <- as.numeric(ICA_weights$ICA_5)
ICA_weights$ICA_6 <- as.numeric(ICA_weights$ICA_6)
ICA_weights$ICA_7 <- as.numeric(ICA_weights$ICA_7)

test1 <- summary(lmer(ICA_weights$ICA_1 ~ Group + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE))
test2 <- summary(lmer(ICA_weights$ICA_2 ~ Group + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE))
test3 <- summary(lmer(ICA_weights$ICA_3 ~ Group + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE))
test4 <- summary(lmer(ICA_weights$ICA_4 ~ Group + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE))
test5 <- summary(lmer(ICA_weights$ICA_5 ~ Group + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE))
test6 <- summary(lmer(ICA_weights$ICA_6 ~ Group + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE))
test7 <- summary(lmer(ICA_weights$ICA_7 ~ Group + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE))

pvals<-data.frame(p=matrix(0,7,1))

pvals[1,1] <- test1$coefficients[2,5]
pvals[2,1] <- test2$coefficients[2,5]
pvals[3,1] <- test3$coefficients[2,5]
pvals[4,1] <- test4$coefficients[2,5]
pvals[5,1] <- test5$coefficients[2,5]
pvals[6,1] <- test6$coefficients[2,5]
pvals[7,1] <- test7$coefficients[2,5]

pvals$p_fdr<-p.adjust(pvals$p, method="fdr")

summary(lmer(ICA_weights$ICA_3 ~ Group + RVI_mean + Sex + SurfaceHoles + Age + (1 | Site), data = data, REML = FALSE))

rh_WM_signal <- read.csv(file="S:/Project/OBIC_dataset/Anders/N848_WM.rh.DK.csv", sep = "", dec = ".")

orig_colnames <- colnames(rh_WM_signal)
new_colnames <- paste("rh_WM", colnames(rh_WM_signal), sep = "_")
colnames(rh_WM_signal) <- new_colnames
colnames(rh_WM_signal)[1] <- "ASEG_ID"
