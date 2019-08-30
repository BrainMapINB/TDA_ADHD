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

# Find preprocessed resting files
# Set same directory as TDA_ADHD200_NYU_pp.R 'outdir'
pp_dir <- ""
pp_ls <- list.files(path = pp_dir, pattern = "res4d_woGSR.nii.gz", full.names = T, recursive = T)

# Read those with bad quality/registration
bad_ls <- scan(file = file.path(getwd(), "01-Preprocessing", "check_reg_ADHD200.txt"),
               what = "character", quiet = T)
# Match IDs
pp_id <- basename(dirname(dirname(dirname(dirname(dirname(pp_ls))))))
pp_r <- basename(dirname(dirname(dirname(pp_ls))))
bad_boys <- match(bad_ls,paste0(pp_id,"_r",substr(pp_r,6,6)), nomatch = F)
# Remove these elements
pp_ls <- pp_ls[-bad_boys]
pp_id <- pp_id[-bad_boys]
pp_r <- pp_r[-bad_boys]

# Set directory where sample info is going to be save
smp_dir <- ""
if(!dir.exists(smp_dir)) dir.create(smp_dir, recursive = T)
# If individual preprocessing info is not stored, create it
if(!file.exists(file.path(smp_dir,"ADHD200_NYU_ppNIHPD.csv"))){
  
  # Load 'fslr' package
  if(!is.element("fslr",row.names(installed.packages()))){
    install.packages("fslr")
  }
  library(fslr)
  
  # Create empty objects in order to record useful info
  pp_n <- length(pp_ls)
  pp_nvol <- vector("integer", pp_n)
  pp_tr <- vector("numeric", pp_n)
  pp_spks <- vector("integer", pp_n)
  pp_avgRMS <- vector("numeric", pp_n)
  # Read files
  for(ii in 1:pp_n){
    
    # Print advance
    print(paste0("File: ",ii,"/",pp_n))
    
    # Read info
    pp_tr[ii] <- scan(file.path(dirname(dirname(pp_ls[ii])), "TR"), quiet = T)
    pp_nvol[ii] <- as.integer(fslval(file = pp_ls[ii], keyword = "dim4", verbose = F))
    pp_spks[ii] <- scan(file.path(dirname(dirname(pp_ls[ii])), "QMov", "NVols_RMSpeaks.txt"), quiet = T)
    pp_avgRMS[ii] <- scan(file.path(dirname(dirname(pp_ls[ii])), "prefiltered_func_data_mcf_rel_mean.rms"),
                          quiet = T)
    
  }
  # Save output
  df <- data.frame(ID = pp_id,
                   Rest = pp_r,
                   TR = pp_tr,
                   NVols = pp_nvol,
                   NSpks = pp_spks,
                   AvgRelRMS = pp_avgRMS)
  write.csv(df, file.path(smp_dir,"ADHD200_NYU_ppNIHPD.csv"),
            quote = F, row.names = F)
} else df <- read.csv(file.path(smp_dir,"ADHD200_NYU_ppNIHPD.csv"))

# Find preprocessed files without ROI's timeseries
ts_ls <- list.files(path = pp_dir, pattern = "pp_woGSR_CC400_ts.txt", full.names = T, recursive = T)
ts_ones <- match(dirname(dirname(ts_ls)),dirname(dirname(pp_ls)))
pp_n <- length(pp_ls)
# If there are preprocessed files without their ROI's timeseries, compute them:
if(length(ts_ones) < pp_n){
  
  # Set indices
  ts_no <- (1:pp_n)[-ts_ones]
  
  # Compute ROI timeseries
  # AFNI's path
  afni_path <- system("which afni", intern = T)
  afni_path <- substr(afni_path, 1, nchar(afni_path)-4)
  # Load 'fslr' package
  if(!is.element("fslr",rownames(installed.packages()))){
    install.packages("fslr")
  }
  library(fslr)
  # Atlas' filenames
  atlas_AAL <- file.path(getwd(),"01-Preprocessing","Atlas","AAL","ROI_MNI_V4.nii")
  atlas_CC200 <- file.path(getwd(),"01-Preprocessing","Atlas","CC200","ADHD200_parcellate_200_2mm.nii.gz")
  atlas_P264 <- file.path(getwd(),"01-Preprocessing","Atlas","P264","power264_2mm.nii.gz")
  atlas_CC400 <- file.path(getwd(),"01-Preprocessing","Atlas","CC400","ADHD200_parcellate_400_2mm.nii.gz")
  atlas <- c(atlas_AAL,atlas_CC200,atlas_P264,atlas_CC400)
  names(atlas) <- c("AAL","CC200","P264","CC400")
  
  for(ii in ts_no){
    
    # Print advance
    print(paste0("File: ",ii,"/",pp_n))
    
    # Create temporary directory
    tmp_dir <- paste0("/tmp/ConnMx_",round(runif(1,100,999)))
    while(dir.exists(tmp_dir)) tmp_dir <- paste0("/tmp/ConnMx_",round(runif(1,100,999)))
    dir.create(tmp_dir)
    
    # Transform Hz to sigma
    hp <- 1/(2*sqrt(2*log(2))*df$TR[ii]*0.01)
    lp <- 1/(2*sqrt(2*log(2))*df$TR[ii]*0.08)
    
    # Band-pass filter
    fslmaths(file = pp_ls[ii],
             outfile = file.path(tmp_dir,"pp_bptf"),
             opts = paste("-bptf",hp,lp))
    
    # Spatial normalization
    fsl_applywarp(infile = file.path(tmp_dir,"pp_bptf"),
                  reffile = file.path(getwd(),"01-Preprocessing","ppSAS","template","t1w_brain_2mm.nii.gz"),
                  warpfile = file.path(dirname(dirname(pp_ls[ii])),"reg","highres2standard_warp.nii.gz"),
                  outfile = file.path(tmp_dir,"pp_norm"),
                  opts = paste0("--premat=",
                                file.path(dirname(dirname(pp_ls[ii])),"reg","meanfunc2highres.mat"),
                                " --interp=trilinear"))
    
    # Extract ROI timeseries
    out_dir <- file.path(dirname(dirname(pp_ls[ii])),"ts")
    if(!dir.exists(out_dir)) dir.create(out_dir, recursive = T)
    
    # Extract temporal series for each region of interest (ROI)
    for(jj in 1:4){
      
      # Print atlas name
      print(names(atlas)[jj])
      
      # Set command line
      cmd <- paste0(
        afni_path,"3dROIstats -mask ",atlas[jj],
        " -quiet ",file.path(tmp_dir,"pp_norm.nii.gz"),
        " > ",file.path(out_dir,paste0("pp_woGSR_",names(atlas)[jj],"_ts.txt"))
      )
      system(cmd)
      
    }
    
    # Remove temporary directory
    unlink(tmp_dir,recursive=T)
  }
  
}
