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

# Read preprocessed files list
# Set same directory as TDA_ADHD200_NYU_pp.R 'outdir'
pp_dir <- ""
# Set directory where sample info is going to be save (same as TDA_ADHD200_NYU_pp2atlas.R)
smp_dir <- ""
datos <- read.csv(file.path(smp_dir,"ADHD200_NYU_ppNIHPD.csv"), colClasses = "character")

# Atlas labels
atlas <- c("AAL","CC200","P264","CC400")

# Compute area, kurtosis and slope from Rips persistence diagram
# Load 'e1017' package
if(!is.element("e1071",row.names(installed.packages()))){
  install.packages("e1071")
}
library(e1071)
# Add slope function
slope <- function(y){
  x <- c(1:length(y))
  coefficients(lm(y~x))[2]
}
# Read Rips function
source(file.path(getwd(),"02-TDA","TDA_RipsDiagram.R"))

# Read AAL file info
aal <- read.table(file.path(getwd(),"01-Preprocessing","atlas","AAL","ROI_MNI_V4.txt"))
# Generate lobule category
aal$lob <- floor(aal$V3/1000)
aal$lob[which(aal$lob==3)] <- 4 # Add insula to paralimbic areas
aal$lob <- factor(aal$lob)
levels(aal$lob) <- c("Fr","Limb","Occ","Par","Sub","Temp","Cbl")

# Read P264 excell with functional network info
# Read ROI and Network labels
# Load 'readxl' package to read .xlsx file
if(!is.element("readxl",rownames(installed.packages()))){
  install.packages("readxl")  
}
library(readxl)
# ROI coordinates can be found here:
# https://www.jonathanpower.net/2011-neuron-bigbrain.html
xls <- read_xlsx(file.path(getwd(),"01-Preprocessing","atlas","P264","Neuron_consensus_264.xlsx"),skip=1)
net_labs <- factor(xls$X__12)
# Abbreviation
levels(net_labs) <- c("AUD","CBL","CinOp","DMN","DAN","FPN","MEM","SAL","SMN-H","SMN-M","SUB","NO","VAN","VIS")

# Create empty matrix to store area, kurtosis and slope
# ncol = 3*(4 atlases + 7 intra-LOB + 21 inter-LOB + 13 intra-FN + 78 inter-FN) = 369
tda <- matrix(0, nrow = nrow(datos), ncol = 369)
# Column names
tda_lab <- c("area","kurt","slop")
colnames(tda) <- colnames(tda, do.NULL = FALSE)
# Atlas
colnames(tda)[1:12] <- paste0(tda_lab,rep(atlas,each=3))
# Anatomical Lobules - intra
colnames(tda)[13:33] <- paste0(tda_lab,rep(levels(aal$lob),each=3))
# Anatomical Lobules - inter
aal_comb <- combn(1:7,2)
colnames(tda)[34:96] <- paste0(tda_lab,rep(paste0(levels(aal$lob)[aal_comb[1,]],"-",levels(aal$lob)[aal_comb[2,]]),each=3))
# Functional Networks - intra
colnames(tda)[97:135] <- paste0(tda_lab,rep(levels(net_labs)[-12],each=3))
# Functional Networks - inter
p264_comb <- combn(1:13,2)
colnames(tda)[136:369] <- paste0(tda_lab,rep(paste0((levels(net_labs)[-12])[p264_comb[1,]],"-",(levels(net_labs)[-12])[p264_comb[2,]]),each=3))

# Search Rips filtration files and compute TDA metrics
for(ii in 1:nrow(datos)){
    
  # Print progress
  print(paste0("File: ",ii,"/",nrow(datos)))
  
  # Find individual TDA directory
  tda_dir <- file.path(pp_dir, "NYU", datos$ID[ii], "session_1", datos$Rest[ii],
                      "rest.rsfMRIv2_nihpd", "tda")
  if(!dir.exists(tda_dir)) stop("Directory not founded!")
  
  # List files
  tda_ls <- list.files(tda_dir)
  
  # Find AAL atlas
  tda_idx <- grep("AAL",tda_ls)
  if(length(tda_idx)==1){
    
    # Read Rips filtration (add filtration at distance zero)
    rips <- c(simplify2array(read.table(file.path(tda_dir,tda_ls[tda_idx]))),0)
    
    # Store TDA metrics
    # Area
    tda[ii,1] <- sum(rips)
    # Kurtosis
    tda[ii,2] <- kurtosis(rips)
    # Slope
    tda[ii,3] <- slope(rips)
    
    # Re-read connectivity matrix
    cmx <- cor(read.table(file.path(dirname(tda_dir),"ts","pp_woGSR_AAL_ts.txt")))
    if(sum(is.na(cmx)) > 0) cmx[is.na(cmx)] <- 0
    
    # Extract intra-FC sub-networks
    for(jj in 1:7){
      # Extract network positions
      intra_idx <- which(aal$lob==levels(aal$lob)[jj])
      
      # Create sub-connectivity matrix
      smx <- cmx[intra_idx,intra_idx]
      
      # Apply Rips diagram
      rips <- c(rips_hom(smx),0)
      
      # Store TDA metrics
      # Area
      tda[ii,10+jj*3] <- sum(rips)
      # Kurtosis
      tda[ii,11+jj*3] <- kurtosis(rips)
      # Slope
      tda[ii,12+jj*3] <- slope(rips)
    }
    
    # Extract inter-FC sub-networks
    for(jj in 1:ncol(aal_comb)){
      # Extract network positions
      net1_lab <- levels(aal$lob)[aal_comb[1,jj]]
      net2_lab <- levels(aal$lob)[aal_comb[2,jj]]
      net1_idx <- which(aal$lob==net1_lab)
      net2_idx <- which(aal$lob==net2_lab)
      inter_idx <- c(net1_idx,net2_idx)
      
      # Create sub-connectivity matrix
      smx <- cmx[inter_idx,inter_idx]
      
      # Apply Rips diagram
      rips <- c(rips_hom(smx),0)
      
      # Store TDA metrics
      # Area
      tda[ii,31+jj*3] <- sum(rips)
      # Kurtosis
      tda[ii,32+jj*3] <- kurtosis(rips)
      # Slope
      tda[ii,33+jj*3] <- slope(rips)
    }
    

  }
  
  # Find CC200 atlas
  tda_idx <- grep("CC200",tda_ls)
  if(length(tda_idx)==1){
    
    # Read Rips filtration (add filtration at distance zero)
    rips <- c(simplify2array(read.table(file.path(tda_dir,tda_ls[tda_idx]))),0)
    
    # Store TDA metrics
    # Area
    tda[ii,4] <- sum(rips)
    # Kurtosis
    tda[ii,5] <- kurtosis(rips)
    # Slope
    tda[ii,6] <- slope(rips)
    
  }
  
  # Find CC400 atlas
  tda_idx <- grep("CC400",tda_ls)
  if(length(tda_idx)==1){
    
    # Read Rips filtration (add filtration at distance zero)
    rips <- c(simplify2array(read.table(file.path(tda_dir,tda_ls[tda_idx]))),0)
    
    # Store TDA metrics
    # Area
    tda[ii,10] <- sum(rips)
    # Kurtosis
    tda[ii,11] <- kurtosis(rips)
    # Slope
    tda[ii,12] <- slope(rips)
    
  }
  
  # Find P264 atlas
  tda_idx <- grep("P264",tda_ls)
  if(length(tda_idx)==1){
    
    # Read Rips filtration (add filtration at distance zero)
    rips <- c(simplify2array(read.table(file.path(tda_dir,tda_ls[tda_idx]))),0)
    
    # Store TDA metrics
    # Area
    tda[ii,7] <- sum(rips)
    # Kurtosis
    tda[ii,8] <- kurtosis(rips)
    # Slope
    tda[ii,9] <- slope(rips)
    
    # Re-read connectivity matrix
    cmx <- cor(read.table(file.path(dirname(tda_dir),"ts","pp_woGSR_P264_ts.txt")))
    if(sum(is.na(cmx)) > 0) cmx[is.na(cmx)] <- 0
    
    # Extract intra-FC sub-networks
    for(jj in 1:13){
      # Extract network positions
      intra_idx <- which(net_labs==(levels(net_labs)[-12])[jj])
      
      # Create sub-connectivity matrix
      smx <- cmx[intra_idx,intra_idx]
      
      # Apply Rips diagram
      rips <- c(rips_hom(smx),0)
      
      # Store TDA metrics
      # Area
      tda[ii,94+jj*3] <- sum(rips)
      # Kurtosis
      tda[ii,95+jj*3] <- kurtosis(rips)
      # Slope
      tda[ii,96+jj*3] <- slope(rips)
    }
    
    # Extract inter-FC sub-networks
    for(jj in 1:ncol(p264_comb)){
      # Extract network positions
      net1_lab <- (levels(net_labs)[-12])[p264_comb[1,jj]]
      net2_lab <- (levels(net_labs)[-12])[p264_comb[2,jj]]
      net1_idx <- which(net_labs==net1_lab)
      net2_idx <- which(net_labs==net2_lab)
      inter_idx <- c(net1_idx,net2_idx)

      # Create sub-connectivity matrix
      smx <- cmx[inter_idx,inter_idx]
      
      # Apply Rips diagram
      rips <- c(rips_hom(smx),0)
      
      # Store TDA metrics
      # Area
      tda[ii,133+jj*3] <- sum(rips)
      # Kurtosis
      tda[ii,134+jj*3] <- kurtosis(rips)
      # Slope
      tda[ii,135+jj*3] <- slope(rips)
    }
  }
  
}

# Write results
datos <- cbind(datos,tda)
write.csv(datos, file.path(smp_dir,"ADHD200_NYU_ppNIHPD_TDA.csv"), quote = F, row.names = F)
