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

# Input parameters
# Set input path, with .zip files from
# https://fcon_1000.projects.nitrc.org/indi/adhd200/
in_dir <- ""
# Set output directory
out_dir <- ""

# Read .zip files
gz_ls <- list.files(path = in_dir, pattern = ".tar.gz", full.names = T)
gz_n <- length(gz_ls)
if(gz_n == 0) stop("No .tar.gz files found!!")

# Uncompress files
if(!dir.exists(out_dir)) dir.create(path = out_dir, recursive = T)
for(ii in gz_ls) untar(tarfile = ii, exdir = out_dir)

# List every resting file
res_ls <- list.files(path = out_dir, pattern = "rest.nii.gz", full.names = T, recursive = T)
res_ls <- res_ls[grep(pattern = "rest_1", res_ls)]
# Generate the corresponding T1w filepaths
t1_ls <- file.path(dirname(dirname(res_ls)),"anat_1")
for(ii in 1:length(t1_ls)) if(length(dir(t1_ls[ii])) == 1) t1_ls[ii] <- file.path(t1_ls[ii],dir(t1_ls[ii]))
# Check if they exists
if(!all(file.exists(t1_ls))) stop("T1w filenames do not correspond!!")

# Find not preprocessed files
res_pp <- paste0(unlist(strsplit(res_ls, ".nii.gz")), ".rsfMRIv2_nihpd/stats/res4d_woGSR.nii.gz")
no_pp <- which(!file.exists(res_pp))
# If not preprocessed files left, preprocess those ones
if(length(no_pp) > 0){
  
  # Preprocess MRI files
  # Create log directory
  log_dir <- paste0(out_dir, "log/")
  if(!dir.exists(log_dir)) dir.create(log_dir)
  # Preprocessing script (bash)
  ppSAS <- file.path(getwd(),"01-Preprocessing","ppSAS","ProcessConnRoiv2_INB_nihpd.sh")
  # Launch preprocessing to cluster (it is highly recommended to parallelize this step)
  for(ii in no_pp){
    
    # Set command
    cmd <- paste0(ppSAS," ",res_ls[ii]," ",t1_ls[ii])
    # Send it
    system(cmd)
    
  } 
  
}

# Add ID to check registration pictures
chck_ls <- list.files(path = out_dir, pattern = "checkreg_rest.jpg", full.names = T, recursive = T)
chck_ses <- basename(dirname(dirname(dirname(chck_ls))))
chck_id <- basename(dirname(dirname(dirname(dirname(dirname(chck_ls))))))
# Rename
chck_new <- file.path(dirname(chck_ls), paste0("checkreg_",chck_id,"_",chck_ses,".jpg"))
file.rename(chck_ls, chck_new)
