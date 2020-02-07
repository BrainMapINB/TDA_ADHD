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
  # In fact, create two more: constrain and unconstrain
  permC_dir <- file.path(perm_dir, "Constrain")
  if(!dir.exists(permC_dir)) dir.create(permC_dir)
  permU_dir <- file.path(perm_dir, "Unconstrain")
  if(!dir.exists(permU_dir)) dir.create(permU_dir)

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
    write.table(tda, file.path(permU_dir,paste0("Rips_",names(cmx_list)[aa],"_perm0.csv")),
                sep = ",", quote = F, row.names = F, col.names = F)
  }
  
  # Compute permutations
  set.seed(18900217)
  nperm <- 1000
  # Unconstrain first
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
        # Rewire matrix (unconstrained randomization)
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
      write.table(tda, file.path(permU_dir,paste0("Rips_",names(cmx_list)[aa],"_perm",pp,".csv")),
                  sep = ",", quote = F, row.names = F, col.names = F)
    }# for(aa in 1:length(cmx_list))
  }# for(pp in 1:nperm)
  ahora <- Sys.time()
  print(ahora-antes)
  
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
  print(ahora-antes)
  
  # Calculate permutation averages
  for(aa in 1:length(cmx_list)){
    
    # Read observed Rips filtration
    rips_files <- list.files(path = permU_dir, pattern = paste0(names(cmx_list)[aa],"_perm0.csv"), full.names = T)
    rips_smp <- read.csv(rips_files[1], header = F)
    
    # Mean CI based on bootstrap (B=1000)
    rips_boot <- t(sapply(1:nperm, function(x) colMeans(rips_smp[sample.int(nrow(rips_smp), replace = T),])))
    # Store average and CI
    rips_smp[1,] <- colMeans(rips_smp)
    # Take quantiles from bootstrap replicates.
    rips_smp[2,] <- apply(rips_boot,2,quantile,0.025)
    rips_smp[3,] <- apply(rips_boot,2,quantile,0.975)
    # Save results
    write.table(rips_smp[1:3,], paste0(perm_dir,"/Rips_",names(cmx_list)[aa],"_AvgCI.csv"),
                quote = F, sep = ",", row.names = F, col.names = F)
    
    # Read permuted Rips filtrations and save averages
    # Unconstrain
    rips_boot <- t(sapply(1:nperm, function(x){
      colMeans(read.csv(paste0(permU_dir,"/Rips_",names(cmx_list)[aa],"_perm",x,".csv"), header = F))}))
    write.table(rips_boot, paste0(perm_dir,"/Rips_",names(cmx_list)[aa],"_permU.csv"),
                quote = F, sep = ",", row.names = F, col.names = F)
    # Constrain
    rips_boot <- t(sapply(1:nperm, function(x){
      colMeans(read.csv(paste0(permC_dir,"/Rips_",names(cmx_list)[aa],"_perm",x,".csv"), header = F))}))
    write.table(rips_boot, paste0(perm_dir,"/Rips_",names(cmx_list)[aa],"_permC.csv"),
                quote = F, sep = ",", row.names = F, col.names = F)
  }
  
}# if(!dir.exists(perm_dir))


####################################################################################################################
# Test Betti0 curves stability
####################################################################################################################

# Set a results directory in order to save images (change this path at your convenience)
res_dir <- file.path(getwd(),"04-Results")
if(!dir.exists(res_dir)) dir.create(res_dir)

# Plot Betti0 curves along permutated data
#pdf(file.path(res_dir,"all_permRips.pdf"),7,5)
par(mfrow = c(2,2))
atlas <-c("AAL","CC200","P264","CC400")
for(aa in 1:length(atlas)){
  
  # Read observed data
  rips_smp <- read.csv(paste0(perm_dir,"/Rips_",atlas[aa],"_AvgCI.csv"), header = F)
  plot(unlist(rips_smp[1,]),1:ncol(rips_smp), type = "l", las = 1,
       xlab = "Filtration value", ylab = "B0", main = atlas[aa])
  # Add CI
  lines(unlist(rips_smp[2,]),1:ncol(rips_smp), lty = 2)
  lines(unlist(rips_smp[3,]),1:ncol(rips_smp), lty = 2)
  
  # Unconstrain null model
  rips_smp <- read.csv(paste0(perm_dir,"/Rips_",atlas[aa],"_permU.csv"), header = F)
  # Draw average
  lines(colMeans(rips_smp), 1:ncol(rips_smp), col = "darkolivegreen")
  # Draw CI
  rips_boot <- t(sapply(1:1000, function(x) colMeans(rips_smp[sample.int(nrow(rips_smp), replace = T),])))
  lines(apply(rips_boot,2,quantile,0.025),1:ncol(rips_smp), lty = 2, col = "darkolivegreen")
  lines(apply(rips_boot,2,quantile,0.975),1:ncol(rips_smp), lty = 2, col = "darkolivegreen")
  
  # Constrain null model
  rips_smp <- read.csv(paste0(perm_dir,"/Rips_",atlas[aa],"_permC.csv"), header = F)
  # Draw average
  lines(colMeans(rips_smp), 1:ncol(rips_smp), col = "darkgoldenrod")
  # Draw CI
  rips_boot <- t(sapply(1:1000, function(x) colMeans(rips_smp[sample.int(nrow(rips_smp), replace = T),])))
  lines(apply(rips_boot,2,quantile,0.025),1:ncol(rips_smp), lty = 2, col = "darkgoldenrod")
  lines(apply(rips_boot,2,quantile,0.975),1:ncol(rips_smp), lty = 2, col = "darkgoldenrod")
}
#dev.off()

# Plot Betti0 variables along permutated data
# Compute area, kurtosis and slope from Rips persistence diagram
# Load 'e1017' package
if(!is.element("e1071",row.names(installed.packages()))){
  install.packages("e1071")
}
library(e1071)
library(ggplot2)
library(gridExtra)
# Add slope function
slope <- function(y){
  x <- c(1:length(y))
  coefficients(lm(y~x))[2]
}
# Plot it
#pdf(file.path(res_dir,"all_permTDA.pdf"),7,5)
par(mfrow = c(2,2))
glist <- vector("list", length(atlas))
for(aa in 1:length(atlas)){
  
  # Unconstrain
  rips_smp <- read.csv(paste0(perm_dir,"/Rips_",atlas[aa],"_permU.csv"), header = F)
  # Calculate variables
  Au <- apply(rips_smp, 1, sum)
  Ku <- apply(rips_smp, 1, kurtosis)
  Su <- apply(rips_smp, 1, slope)
  
  # Constrain
  rips_smp <- read.csv(paste0(perm_dir,"/Rips_",atlas[aa],"_permC.csv"), header = F)
  # Calculate variables
  Ac <- apply(rips_smp, 1, sum)
  Kc <- apply(rips_smp, 1, kurtosis)
  Sc <- apply(rips_smp, 1, slope)
  
  # Add to data.frame
  tda_smp <- as.data.frame(cbind(Au, Ac, Ku, Kc, Su, Sc))
  # Scale based on observed data
  idx <- which(names(datos)==paste0("area",atlas[aa]))
  tda_smp$Au <- (tda_smp$Au - mean(datos[,idx]))/sd(datos[,idx])
  tda_smp$Ac <- (tda_smp$Ac - mean(datos[,idx]))/sd(datos[,idx])
  idx <- which(names(datos)==paste0("kurt",atlas[aa]))
  tda_smp$Ku <- (tda_smp$Ku - mean(datos[,idx]))/sd(datos[,idx])
  tda_smp$Kc <- (tda_smp$Ku - mean(datos[,idx]))/sd(datos[,idx])
  idx <- which(names(datos)==paste0("slop",atlas[aa]))
  tda_smp$Su <- (tda_smp$Su - mean(datos[,idx]))/sd(datos[,idx])
  tda_smp$Sc <- (tda_smp$Sc - mean(datos[,idx]))/sd(datos[,idx])
  # Boxplot
  boxplot(tda_smp, border = c("darkolivegreen","darkgoldenrod"), main = atlas[aa], pch = ".")
  abline(h = 0, lty = 2)
  
  # Violin plot
  # Reshape data.frame
  tda_rng <- data.frame(avg = colMeans(tda_smp),
                        min = apply(tda_smp,2,min),
                        max = apply(tda_smp,2,min),
                        tda = c("Au","Ac","Ku","Kc","Su","Sc"),
                        null = rep(c("U","C"),3))
  # Plot
  # Most basic error bar
  glist[[aa]] <- ggplot(tda_rng, aes(x = tda, y = avg)) +
    geom_bar(aes(x=tda, y=avg), stat="identity", fill="gray", alpha=0.3) +
    geom_errorbar(aes(x=tda, ymin=min, ymax=max, color = null), width=0.4, alpha=0.9, size=1.3) +
    geom_hline(yintercept=0, linetype="dashed", color = "black") + 
    ggtitle(atlas[aa]) + xlab("") + ylab("z") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
}
#dev.off()
grid.arrange(glist[[1]], glist[[2]], glist[[3]], glist[[4]])
