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

# Find preprocessed resting files ROI timeseries and their corresponding Rips filtrations
# Set same directory as TDA_ADHD200_NYU_pp.R 'outdir'
pp_dir <- ""
ts_ls <- list.files(path = pp_dir, pattern = "_ts.txt", full.names = T, recursive = T)
ts_n <- length(ts_ls)
rips_ls <- list.files(path = pp_dir, pattern = "_Rips.txt", full.names = T, recursive = T)
# Set those without Rips filtration
if(length(rips_ls)==0) rips_no <- 1:ts_n else rips_no <- (1:ts_n)[-match(unlist(strsplit(ts_ls, "_ts.txt")),gsub("/tda/", "/ts/", unlist(strsplit(rips_ls, "_Rips.txt"))))]

# If there are timeseries files without their Rips filtrations, compute them:
if(length(rips_no) > 0){
  
  # Read Rips function
  source(file.path(getwd(),"02-TDA","TDA_RipsDiagram.R"))
  
  # Compute Rips filtration
  for(ii in rips_no){
    
    # Print advance
    print(paste0("File: ",ii,"/",ts_n))
    
    # Compute connectivity matrix (via Pearson's R)
    cmx <- cor(read.table(ts_ls[ii]))
    if(sum(is.na(cmx)) > 0) cmx[is.na(cmx)] <- 0
    
    # Compute Rips filtration
    rips <- rips_hom(cmx)
    
    # Create output directory
    out_dir <- file.path(dirname(dirname(ts_ls[ii])),"tda")
    if(!dir.exists(out_dir)) dir.create(out_dir, recursive = T)
    
    # Write output
    write.table(x = rips, quote = F, row.names = F, col.names = F,
                file = file.path(out_dir,paste0(unlist(strsplit(basename(ts_ls[ii]),"_ts.txt")),"_Rips.txt")))
    
  }
  
}

# Concatenate every Rips filtration
# Extract ID and fMRI scan number from Rips list
id <- basename(dirname(dirname(dirname(dirname(dirname(rips_ls))))))
rest <- basename(dirname(dirname(dirname(rips_ls))))
# Get atlas labels
atlas <- c("AAL","CC200","P264","CC400")
for(ii in 1:length(atlas)){
  
  # Find corresponding indexes
  atlas_idx <- grep(atlas[ii], rips_ls)
  
  # Read and concatenate those files
  rips <- t(sapply(atlas_idx, function(x) scan(rips_ls[x], quiet = T)))
  
  # Create dataframe
  df <- data.frame(ID=id[atlas_idx], rest=rest[atlas_idx])
  df <- cbind(df, rips)
  
  # Write output
  write.csv(df, file.path(getwd(),"02-TDA",paste0("ADHD200_NYU_ppNIHPD_Rips_",atlas[ii],".csv")))
  
}
