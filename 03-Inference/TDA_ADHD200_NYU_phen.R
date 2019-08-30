#!/usr/bin/env Rscript

# Clear workspace
rm(list = ls())
# Set working directory
# Load 'rstudioapi' package
if(!is.element("rstudioapi",row.names(installed.packages()))){
  install.packages("rstudioapi")
}
library(rstudioapi)
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
# Save graphical parameters
op <- par()

# Read phenotypic files
# Set input path, with phenotypic .csv files from
# https://fcon_1000.projects.nitrc.org/indi/adhd200/
phen_dir <- ""
# Concatenate train and test datasets
phen <- rbind(read.csv(file.path(phen_dir, "NYU_phenotypic.csv")),
              read.csv(file.path(phen_dir, "NYU_TestRelease_phenotypic.csv")))
# Set NAs
phen[phen==""] <- NA
phen[phen=="N/A"] <- NA
phen[phen=="-999"] <- NA
# Read TDA file
# Set directory where sample info is going to be save (same as TDA_ADHD200_NYU_pp2atlas.R)
smp_dir <- ""
tda <- read.csv(file.path(smp_dir, "ADHD200_NYU_ppNIHPD_TDA.csv"))
# Concatenate datasets
datos <- cbind(phen[match(tda$ID,phen$ScanDir.ID),],tda)

# Write data
tda_dir <- file.path(getwd(),"03-Inference")
write.csv(datos, file.path(tda_dir, "ADHD200_NYU_ppNIHPD_phenTDA.csv"), row.names = F)
