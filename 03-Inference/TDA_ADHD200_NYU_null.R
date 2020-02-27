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

# Read data
datos <- read.csv(file.path(getwd(),"03-Inference", "ADHD200_NYU_ppNIHPD_phenTDA.csv"))

# Set head motion quality control (based on Satterhwaite et al. 2013)
# Remove those with less than 240s (4min) of clean timeseries
qc_idx <- which(datos$TR*(datos$NVols-datos$NSpks) > 240)
datos <- datos[qc_idx,]
# Remove those with less than 0.55mm of RMS relative displacement
qc_idx <- which(datos$AvgRelRMS < 0.55)
datos <- datos[qc_idx,]
# Remove those with secondary disorder
datos <- datos[which(is.na(datos$Secondary.Dx.)),]
# Remove those with not naive medication status
datos <- datos[-which(datos$Med.Status == 2),]
# Remove those with less than 80 IQ
datos <- datos[-which(datos$Full4.IQ < 80),]
# Remove controls with higher ADHD index
datos <- datos[-which(datos$ADHD.Index > 60 & datos$DX == 0),]
# Remove subjects with uncomplete info
datos <- datos[which(!is.na(datos$Gender)),]

# Set experimental groups
datos$grupo <- datos$DX
datos$grupo[which(datos$grupo>1)] <- 1
datos$grupo <- factor(datos$grupo)
levels(datos$grupo) <- c("TDC","ADHD")
datos$Sex <- factor(datos$Gender)
levels(datos$Sex) <- c("F","M")

####################################################################################################################
# Permutate by network rewiring
####################################################################################################################

# Read timeseries files
ts_dir <- file.path(getwd(),"02-TDA","TS")
if(!dir.exists(ts_dir)){
  # If directory does not exist, then create it
  dir.create(ts_dir)
  
  # Store timeseries files (cured ones)
  # Set same directory as TDA_ADHD200_NYU_pp2atlas.R 'ppdir'
  pp_dir <- 
  for(aa in c("AAL","CC200","P264","CC400")){
    # Create directory for each atlas
    if(!dir.exists(file.path(ts_dir,aa))) dir.create(file.path(ts_dir,aa))

    # Copy files
    for(ii in 1:nrow(datos)){
      # Find file
      id <- as.character(datos$ID[ii])
      if(nchar(id)==5) id <- paste0("00",id)
      ts_file <- file.path(pp_dir,"NYU",id,"session_1",datos$Rest[ii],"rest.rsfMRIv2_nihpd","ts",paste0("pp_woGSR_",aa,"_ts.txt"))
      # Copy file
      file.copy(ts_file, file.path(ts_dir,aa,paste0(id,"_",aa,"_ts.txt")))
    }
  }
}

# Permute
perm_dir <- file.path(getwd(),"02-TDA","PERM")
if(!dir.exists(perm_dir)){
  # If directory does not exist, then create it
  dir.create(perm_dir)
  # In fact, create two more: constrain and Weighted
  permC_dir <- file.path(perm_dir, "Constrain")
  if(!dir.exists(permC_dir)) dir.create(permC_dir)
  permW_dir <- file.path(perm_dir, "Weighted")
  if(!dir.exists(permW_dir)) dir.create(permW_dir)

  # Compute connectivity matrices
  cmx_list <- vector("list", 4)
  names(cmx_list) <- c("AAL","CC200","P264","CC400")
  for(ii in 1:length(cmx_list)){
    # List timeseries files
    ts_files <- list.files(file.path(ts_dir,names(cmx_list)[ii]), full.names = T)
    cmx_list[[ii]] <- sapply(1:length(ts_files), function(x) cor(read.table(ts_files[x])))
    if(sum(is.na(cmx_list[[ii]])) > 0) cmx_list[[ii]][is.na(cmx_list[[ii]])] <- 0
    # Reshape object
    nnodes <- sqrt(dim(cmx_list[[ii]])[1])
    cmx_list[[ii]] <- array(cmx_list[[ii]], dim = c(nnodes,nnodes,length(ts_files)))
  }
    
  # Load 'parallel' library
  library(parallel)
  # Compute TDA features
  # Read Rips function
  source(file.path(getwd(),"02-TDA","TDA_RipsDiagram.R"))
  
  # Generate observed values (permutation 0)
  for(aa in 1:length(cmx_list)){
    print(names(cmx_list)[aa])
    # Sapply rewiring
    cl <- makeCluster(detectCores())
    clusterExport(cl=cl,
                  varlist=c("cmx_list", "ts_files", "aa", "rips_hom", "ripsDiag"),
                  envir=environment())
    tda <- t(parSapply(cl, 1:length(ts_files), function(x){
      # Extract matrix
      cmx <- cmx_list[[aa]][,,x]
      # Compute Rips filtration
      rips_hom(cmx)
    }))
    stopCluster(cl)
    write.table(tda, file.path(permC_dir,paste0("Rips_",names(cmx_list)[aa],"_perm0.csv")),
                sep = ",", quote = F, row.names = F, col.names = F)
    write.table(tda, file.path(permW_dir,paste0("Rips_",names(cmx_list)[aa],"_perm0.csv")),
                sep = ",", quote = F, row.names = F, col.names = F)
  }
  
  # Compute permutations
  set.seed(18900217)
  nperm <- 1000
  # Weighted first
  antes <- Sys.time()
  for(pp in 1:nperm){
    if(pp%%100==0) print(pp)
    for(aa in 1:length(cmx_list)){
      # Sapply rewiring
      cl <- makeCluster(detectCores())
      clusterExport(cl=cl,
                    varlist=c("cmx_list", "ts_files", "aa", "rips_hom", "ripsDiag"),
                    envir=environment())
      tda <- t(parSapply(cl, 1:length(ts_files), function(x){
        # Rewire matrix (Weighted randomization)
        cmx <- cmx_list[[aa]][,,x]
        diag(cmx) <- 0
        cmx[upper.tri(cmx)] <- sample(cmx[upper.tri(cmx)])
        cmx[lower.tri(cmx)] <- 0
        cmx <- cmx + t(cmx)
        diag(cmx) <- 1
        # Compute Rips filtration
        rips_hom(cmx)
      }))
      stopCluster(cl)
      write.table(tda, file.path(permW_dir,paste0("Rips_",names(cmx_list)[aa],"_perm",pp,".csv")),
                  sep = ",", quote = F, row.names = F, col.names = F)
    }# for(aa in 1:length(cmx_list))
  }# for(pp in 1:nperm)
  ahora <- Sys.time()
  print(ahora-antes) # Time difference of 4.539018 days (in one computer Intel Core i7-4790 CPU @ 3.60GHz × 8 with Ubuntu 18.04.3 LTS 64-bit)
  
  # Constrain ones
  antes <- Sys.time()
  # Load 'igraph' package
  if(!is.element("igraph",row.names(installed.packages()))){
    install.packages("igraph")
  }
  library(igraph)
  # Compute permutations
  for(pp in 1:nperm){
    if(pp%%100==0) print(pp)
    for(aa in 1:length(cmx_list)){
      # Sapply rewiring
      cl <- makeCluster(detectCores())
      clusterExport(cl=cl,
                    varlist=c("cmx_list", "ts_files", "aa", "rips_hom", "ripsDiag",
                              "graph.adjacency", "rewire", "keeping_degseq", "vcount",
                              "get.edgelist", "E"),
                    envir=environment())
      tda <- t(parSapply(cl, 1:length(ts_files), function(x){
        # Remove diagonal
        cmx <- cmx_list[[aa]][,,x]
        diag(cmx) <- 0
        cmx[cmx<0] <- 0
        # Convert connectivity matrix to igraph object
        ig <- graph.adjacency(cmx, mode="undirected", weighted=T)
        # Rewire (keeping degree sequence)
        rg <- rewire(ig, with=keeping_degseq(niter=round(vcount(ig)^1.5)))
        # Converto igraph to matrix
        cmx <- matrix(0, vcount(ig), vcount(ig))
        cmx[get.edgelist(rg)] <- E(ig)$weight
        cmx <- cmx + t(cmx)
        diag(cmx) <- 1
        # Compute Rips filtration
        rips_hom(cmx)
      }))
      stopCluster(cl)
      write.table(tda, file.path(permC_dir,paste0("Rips_",names(cmx_list)[aa],"_perm",pp,".csv")),
                  sep = ",", quote = F, row.names = F, col.names = F)
    }# for(aa in 1:length(cmx_list))
  }# for(pp in 1:nperm)
  ahora <- Sys.time()
  print(ahora-antes) # Time difference of 12.5 days (in one computer)
  
  # Calculate permutation averages and CI
  for(aa in 1:length(cmx_list)){
    
    # Read observed Rips filtration
    rips_files <- list.files(path = permW_dir, pattern = paste0(names(cmx_list)[aa],"_perm0.csv"), full.names = T)
    rips_smp <- read.csv(rips_files[1], header = F)
    
    # Mean CI based on bootstrap (B=1000)
    rips_boot <- t(sapply(1:nperm, function(x) colMeans(rips_smp[sample.int(nrow(rips_smp), replace = T),])))
    # Store average and CI
    rips_smp[1,] <- colMeans(rips_smp)
    # Take quantiles from bootstrap replicates.
    rips_smp[2,] <- apply(rips_boot,2,quantile,0.025)
    rips_smp[3,] <- apply(rips_boot,2,quantile,0.975)
    
    # Read permuted Rips filtrations and save averages and CI
    # Weighted
    rips_boot <- t(sapply(1:nperm, function(x){
      colMeans(read.csv(paste0(permW_dir,"/Rips_",names(cmx_list)[aa],"_perm",x,".csv"), header = F))}))
    rips_smp[4,] <- colMeans(rips_boot)
    rips_smp[5,] <- apply(rips_boot,2,quantile,0.025)
    rips_smp[6,] <- apply(rips_boot,2,quantile,0.975)
    
    # Save results
    rips_smp <- rips_smp[1:6,]
    rownames(rips_smp) <- c("OBSavg","OBSlowCI","OBSupCI",
                            "PERMavg","PERMlowCI","PERMupCI")
    write.table(rips_smp, paste0(perm_dir,"/Rips_",names(cmx_list)[aa],"_AvgCI.csv"),
                quote = F, sep = ",", col.names = F)
  }
  
  # Calculate permutated TDA variables
  tda_dir <- file.path(perm_dir,"TDA")
  if(!dir.exists(tda_dir)) dir.create(tda_dir)
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
  # Calculate permutated TDA variables
  for(aa in 1:length(cmx_list)){
    
    # Empty object
    perm_Aw <- perm_Ac <- matrix(0, nrow = nrow(datos), ncol = nperm)
    perm_Kw <- perm_Kc <- matrix(0, nrow = nrow(datos), ncol = nperm)
    perm_Sw <- perm_Sc <- matrix(0, nrow = nrow(datos), ncol = nperm)
    
    # Read permuted Rips filtrations and save TDA variables
    for(pp in 1:nperm){
      # Weighted
      rips_smp <- read.csv(paste0(permW_dir,"/Rips_",names(cmx_list)[aa],"_perm",pp,".csv"), header = F)
      perm_Aw[,pp] <- apply(rips_smp, 1, sum)
      perm_Kw[,pp] <- apply(rips_smp, 1, kurtosis)
      perm_Sw[,pp] <- apply(rips_smp, 1, slope)
      
      # Constrain
      rips_smp <- read.csv(paste0(permC_dir,"/Rips_",names(cmx_list)[aa],"_perm",pp,".csv"), header = F)
      perm_Ac[,pp] <- apply(rips_smp, 1, sum)
      perm_Kc[,pp] <- apply(rips_smp, 1, kurtosis)
      perm_Sc[,pp] <- apply(rips_smp, 1, slope)
    }
    
    # Save results
    write.table(perm_Aw, paste0(tda_dir,"/permW_area_",names(cmx_list)[aa],".csv"),
                quote = F, sep = ",", col.names = F, row.names = F)
    write.table(perm_Kw, paste0(tda_dir,"/permW_kurt_",names(cmx_list)[aa],".csv"),
                quote = F, sep = ",", col.names = F, row.names = F)
    write.table(perm_Sw, paste0(tda_dir,"/permW_slop_",names(cmx_list)[aa],".csv"),
                quote = F, sep = ",", col.names = F, row.names = F)
    write.table(perm_Ac, paste0(tda_dir,"/permC_area_",names(cmx_list)[aa],".csv"),
                quote = F, sep = ",", col.names = F, row.names = F)
    write.table(perm_Kc, paste0(tda_dir,"/permC_kurt_",names(cmx_list)[aa],".csv"),
                quote = F, sep = ",", col.names = F, row.names = F)
    write.table(perm_Sc, paste0(tda_dir,"/permC_slop_",names(cmx_list)[aa],".csv"),
                quote = F, sep = ",", col.names = F, row.names = F)
  }
  
  # Extract permuted p-values
  # Create empty object to store results
  res <- matrix(1, nrow = 5, ncol = 6)
  rownames(res) <- c(names(cmx_list),"TOTAL")
  colnames(res) <- c("areaW","areaC","kurtW","kurtW","slopW","slopC")
  # Compute probabilities versus observed value
  for(aa in 1:length(names(cmx_list))){
    for(bb in 1:ncol(res)){
      # Extract variable labels
      varname <- substr(colnames(res)[bb],1,4)
      nullname <- substr(colnames(res)[bb],5,5)
      # Find files
      tda <- read.csv(paste0(tda_dir,"/perm",nullname,"_",varname,"_",names(cmx_list)[aa],".csv"), header = F)
      idx <- which(names(datos) == paste0(varname,names(cmx_list)[aa]))
      # Count surpassed ones
      res[aa,bb] <- sum(sum((abs(tda)-abs(datos[,idx]))>0)/nrow(tda))/nperm
    }
  }
  # Get variable average
  res[5,] <- colMeans(res[1:length(names(cmx_list)),])
  # Write results
  write.csv(res, file.path(perm_dir,"TDA_PERMpval.csv"), quote = F)

}# if(!dir.exists(perm_dir))
