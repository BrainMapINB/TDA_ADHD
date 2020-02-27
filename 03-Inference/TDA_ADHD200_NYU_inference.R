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
# Sample inference
####################################################################################################################

# Set experimental groups
datos$grupo <- datos$DX
datos$grupo[which(datos$grupo>1)] <- 1
datos$grupo <- factor(datos$grupo)
levels(datos$grupo) <- c("TDC","ADHD")
datos$Sex <- factor(datos$Gender)
levels(datos$Sex) <- c("F","M")

#------------------------------------------------------------------------------------------------------------------
# Describe dataset
#------------------------------------------------------------------------------------------------------------------

# Load 'psych' package
if(!is.element("psych",row.names(installed.packages()))){
  install.packages("psych")
}
library(psych)

# Descriptive statistics
# N
table(datos$grupo)
# Sex
(t <- table(datos$grupo,datos$Sex))
fisher.test(t) # Odds Ratio and significance
# Age
describeBy(datos$Age, datos$grupo)
t.test(Age~grupo, data=datos, var.equal = T)
DF <- cbind(datos$Age,datos$grupo); colnames(DF) <- c("Age","grupo")
cohen.d(DF,"grupo")$cohen.d
# ADHD index
describeBy(datos$ADHD.Index, datos$grupo)
t.test(ADHD.Index~grupo, data=datos, var.equal = T)
DF <- cbind(datos$ADHD.Index,datos$grupo); colnames(DF) <- c("adhd","grupo")
cohen.d(DF,"grupo")$cohen.d
# IQ
describeBy(datos$Full4.IQ, datos$grupo)
t.test(Full4.IQ~grupo, data=datos, var.equal = T)
DF <- cbind(datos$Full4.IQ,datos$grupo); colnames(DF) <- c("iq","grupo")
cohen.d(DF,"grupo")$cohen.d
# Average RMS Relative Displacement
describeBy(datos$AvgRelRMS, datos$grupo, digits = 4, mat = T)
t.test(AvgRelRMS~grupo, data=datos, var.equal = T)
DF <- cbind(datos$AvgRelRMS,datos$grupo); colnames(DF) <- c("adhd","grupo")
cohen.d(DF,"grupo")$cohen.d

#------------------------------------------------------------------------------------------------------------------
# Kendall Concordance Coefficient (KCC) along parcellations
#------------------------------------------------------------------------------------------------------------------

# Load 'irr' package
if(!is.element("irr",row.names(installed.packages()))){
  install.packages("irr")
}
library(irr)

atlas <- c("AAL","CC200","P264","CC400")
# Betti0-AREA
# All atlases together
(k <- kendall(cbind(datos$areaAAL,datos$areaCC200,datos$areaP264,datos$areaCC400)))
# Pairwise
area_kcc <- array(0,dim=c(4,4))
# AAL-CC200
(k <- kendall(cbind(datos$areaAAL,datos$areaCC200)))
area_kcc[1,2] <- k$value
# AAL-P264
(k <- kendall(cbind(datos$areaAAL,datos$areaP264)))
area_kcc[1,3] <- k$value
# AAL-CC400
(k <- kendall(cbind(datos$areaAAL,datos$areaCC400)))
area_kcc[1,4] <- k$value
# CC200-P264
(k <- kendall(cbind(datos$areaCC200,datos$areaP264)))
area_kcc[2,3] <- k$value
# CC200-CC400
(k <- kendall(cbind(datos$areaCC200,datos$areaCC400)))
area_kcc[2,4] <- k$value
# P264-CC400
(k <- kendall(cbind(datos$areaP264,datos$areaCC400)))
area_kcc[3,4] <- k$value

# Betti0-KURTOSIS
# All atlases together
(k <- kendall(cbind(datos$kurtAAL,datos$kurtCC200,datos$kurtP264,datos$kurtCC400)))
# Pairwise
kurt_kcc <- array(0,dim=c(4,4))
# AAL-CC200
(k <- kendall(cbind(datos$kurtAAL,datos$kurtCC200)))
kurt_kcc[1,2] <- k$value
# AAL-P264
(k <- kendall(cbind(datos$kurtAAL,datos$kurtP264)))
kurt_kcc[1,3] <- k$value
# AAL-CC400
(k <- kendall(cbind(datos$kurtAAL,datos$kurtCC400)))
kurt_kcc[1,4] <- k$value
# CC200-P264
(k <- kendall(cbind(datos$kurtCC200,datos$kurtP264)))
kurt_kcc[2,3] <- k$value
# CC200-CC400
(k <- kendall(cbind(datos$kurtCC200,datos$kurtCC400)))
kurt_kcc[2,4] <- k$value
# P264-CC400
(k <- kendall(cbind(datos$kurtP264,datos$kurtCC400)))
kurt_kcc[3,4] <- k$value

# Betti0-Slope
# All atlases together
k <- (kendall(cbind(datos$slopAAL,datos$slopCC200,datos$slopP264,datos$slopCC400)))
# Pairwise
slop_kcc <- array(0,dim=c(4,4))
# AAL-CC200
(k <- kendall(cbind(datos$slopAAL,datos$slopCC200)))
slop_kcc[1,2] <- k$value
# AAL-P264
(k <- kendall(cbind(datos$slopAAL,datos$slopP264)))
slop_kcc[1,3] <- k$value
# AAL-CC400
(k <- kendall(cbind(datos$slopAAL,datos$slopCC400)))
slop_kcc[1,4] <- k$value
# CC200-P264
(k <- kendall(cbind(datos$slopCC200,datos$slopP264)))
slop_kcc[2,3] <- k$value
# CC200-CC400
(k <- kendall(cbind(datos$slopCC200,datos$slopCC400)))
slop_kcc[2,4] <- k$value
# P264-CC400
(k <- kendall(cbind(datos$slopP264,datos$slopCC400)))
slop_kcc[3,4] <- k$value

#------------------------------------------------------------------------------------------------------------------
# Plot results
# Load 'corrplot' package
if(!is.element("corrplot",row.names(installed.packages()))){
  install.packages("corrplot")
}
library(corrplot)

# Set a results directory in order to save images (change this path at your convenience)
res_dir <- file.path(getwd(),"04-Results")
if(!dir.exists(res_dir)) dir.create(res_dir)
#pdf(file.path(res_dir,"KCC.pdf"),5,5)
#svg(file.path(res_dir,"KCC.svg"),5,5)
par(mfrow=c(1,3))
# Area
colnames(area_kcc) <- rownames(area_kcc) <- atlas
colorcito <- c(heat.colors(50),rev(heat.colors(50)))
corrplot(area_kcc,method = "square", type = "upper", diag=F, mar=c(2,0,2,0),
         col = colorcito, tl.col = "black", addCoef.col = "black",
         tl.srt = 45, is.corr = F, cl.lim = c(0,1), cl.ratio = .3)
# Kurtosis
colnames(kurt_kcc) <- rownames(kurt_kcc) <- atlas
corrplot(kurt_kcc,method = "square", type = "upper", diag=F, mar=c(2,0,2,0),
         col = colorcito, tl.col = "black", addCoef.col = "black",
         tl.srt = 45, is.corr = F, cl.lim = c(0,1), cl.ratio = .3)
# Slope
colnames(slop_kcc) <- rownames(slop_kcc) <- atlas
corrplot(slop_kcc,method = "square", type = "upper", diag=F, mar=c(2,0,2,0),
         col = colorcito, tl.col = "black", addCoef.col = "black",
         tl.srt = 45, is.corr = F, cl.lim = c(0,1), cl.ratio = .3)
# Save plot
#dev.off()
par(op)


####################################################################################################################
# Group inferences
####################################################################################################################
#------------------------------------------------------------------------------------------------------------------
# Groups comparison - logistic regression
#------------------------------------------------------------------------------------------------------------------

# Create empty object to save every result
res <- array(0, dim=c(20,6))
rownames(res) <- c("AAL-logOR","AAL-logORer","AAL-OR","AAL-z","AAL-p",
                   "CC200-logOR","CC200-logORer","CC200-OR","CC200-z","CC200-p",
                   "P264-logOR","P264-logORer","P264-OR","P264-z","P264-p",
                   "CC400-logOR","CC400-logORer","CC400-OR","CC400-z","CC400-p")
colnames(res) <- c("area","kurt","slop","sex","age","avgFD")

# AAL
fit <- summary(glm(grupo ~ scale(areaAAL) + scale(kurtAAL) + scale(slopAAL) + Sex + scale(Age) + scale(AvgRelRMS),
                   data = datos, family=binomial("logit")))
res[1,] <- fit$coefficients[2:7,1]
res[2,] <- fit$coefficients[2:7,2]
res[3,] <- exp(fit$coefficients[2:7,1])
res[4,] <- fit$coefficients[2:7,3]
res[5,] <- fit$coefficients[2:7,4]
# CC200
fit <- summary(glm(grupo ~ scale(areaCC200) + scale(kurtCC200) + scale(slopCC200) + Sex + scale(Age) + scale(AvgRelRMS),
                   data = datos, family=binomial("logit")))
res[6,] <- fit$coefficients[2:7,1]
res[7,] <- fit$coefficients[2:7,2]
res[8,] <- exp(fit$coefficients[2:7,1])
res[9,] <- fit$coefficients[2:7,3]
res[10,] <- fit$coefficients[2:7,4]
# P264
fit <- summary(glm(grupo ~ scale(areaP264) + scale(kurtP264) + scale(slopP264) + Sex + scale(Age) + scale(AvgRelRMS),
                   data = datos, family=binomial("logit")))
res[11,] <- fit$coefficients[2:7,1]
res[12,] <- fit$coefficients[2:7,2]
res[13,] <- exp(fit$coefficients[2:7,1])
res[14,] <- fit$coefficients[2:7,3]
res[15,] <- fit$coefficients[2:7,4]

# CC400
fit <- summary(glm(grupo ~ scale(areaCC400) + scale(kurtCC400) + scale(slopCC400) + Sex + scale(Age) + scale(AvgRelRMS),
                   data = datos, family=binomial("logit")))
res[16,] <- fit$coefficients[2:7,1]
res[17,] <- fit$coefficients[2:7,2]
res[18,] <- exp(fit$coefficients[2:7,1])
res[19,] <- fit$coefficients[2:7,3]
res[20,] <- fit$coefficients[2:7,4]

# Results
signif(res,3)

# Print results
atlas <- c("AAL","CC200","P264","CC400")
perm_dir <- file.path(getwd(),"02-TDA","PERM")
#pdf(file.path(res_dir,"grp_all_rips.pdf"),7,5)
#svg(file.path(res_dir,"grp_all_rips.svg"),7,5)
par(mfrow = c(2,2))
for(aa in 1:length(atlas)){
  
  # Add Null model first
  rips <- read.csv(paste0(perm_dir,"/Rips_",atlas[aa],"_AvgCI.csv"), header = F, row.names = 1)
  rips$Y <- 0
  # (Empty) plot
  y <- 1:ncol(rips); x <- y/length(y)
  plot(x,y, type = "n", xlim = c(0,0.6), las = 1, frame.plot = F, axes = F,
       ylab = "Betti 0", xlab = "Filtration Value", main = atlas[aa])
  axis(side = 1, lwd = 2)
  axis(side = 2, lwd = 2, las = 1)
  # Average + CI
  polygon(c(unlist(rips[5,]),rev(unlist(rips[6,]))),c(y,rev(y)),col = "gray25", border = FALSE)
  lines(unlist(rips[4,]),y, col = "gray25", lwd = 1)
  
  # Observed data
  filename <- file.path(getwd(),"02-TDA",paste0("ADHD200_NYU_ppNIHPD_Rips_",atlas[aa],".csv"))
  # Read Rips filtration for each participant
  rips <- read.csv(filename)
  # Match IDs and rest scan
  rips <- rips[match(datos$ID,rips$ID),]
  rips <- rips[,-(1:3)]
  # Get group averages
  rips$Y <- rep(0,nrow(rips)) # Add distance zero
  ripsAVG <- describeBy(rips, datos$grupo, mat = T)
  # Draw c.i. light lines
  colorcito <- c("blue", "red")
  # Draw group means
  lines(ripsAVG$mean[which(ripsAVG$group1=="ADHD")], y, col="red", lwd=1.5)
  lines(ripsAVG$mean[which(ripsAVG$group1=="TDC")], y, col="blue", lwd=1.5)
  # ADHD
  mADHD <- ripsAVG$mean[which(ripsAVG$group1=="ADHD")]
  ciADHD <- ripsAVG$se[which(ripsAVG$group1=="ADHD")]*1.96
  lines(mADHD, y, col="red", lwd=1.5)
  #lines(m-ci, y, col="red", lty=2)
  #lines(m+ci, y, col="red", lty=2)
  polygon(c(mADHD-ciADHD,rev(mADHD+ciADHD)),c(y,rev(y)),col = scales::alpha("red", 0.35), border = FALSE)
  # TDC
  mTDC <- ripsAVG$mean[which(ripsAVG$group1=="TDC")]
  ciTDC <- ripsAVG$se[which(ripsAVG$group1=="TDC")]*1.96
  lines(mTDC, y, col="blue", lwd=1.5)
  #lines(m-ci, y, col="blue", lty=2)
  #lines(m+ci, y, col="blue", lty=2)
  polygon(c(mTDC-ciTDC,rev(mTDC+ciTDC)),c(y,rev(y)),col = scales::alpha("blue", 0.35), border = FALSE)
  # Add a legend
  #legend("bottomleft", legend = c("ADHD","TDC"), col = c("red", "blue"), lty=1, cex=0.8)
  # Add dotted rectangle
  #rect(xleft = 0.3, xright = 0.4, ybottom = 9, ytop = 41, lty = 2)
}
#dev.off()
par(op)

# Add forest plots
# Load 'ggplot2' package
if(!is.element("ggplot2",row.names(installed.packages()))){
  install.packages("ggplot2")
}
library(ggplot2)
gglist <- vector("list",length(atlas))
for(gg in 1:length(gglist)){
  # Generate odds ratio data.frame
  ORdf <- data.frame(OR = res[3+5*(gg-1),1:3],
                     lowCI = exp(res[1+5*(gg-1),1:3]-1.96*res[2+5*(gg-1),1:3]),
                     upCI = exp(res[1+5*(gg-1),1:3]+1.96*res[2+5*(gg-1),1:3]),
                     var = c("A","K","S"))
  # Genearte plot
  gglist[[gg]] <- ggplot(ORdf, aes(x = OR, y = var)) +
    geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
    geom_errorbarh(aes(xmax = upCI, xmin = lowCI), size = .5, height = 
                     .2, color = "gray50") +
    geom_point(size = 2, color = "red", shape = 19) +
    scale_y_discrete(limits = rev(unique(sort(ORdf$var))), position = "right") +
    scale_x_continuous(position = "top") +
    theme_bw(base_rect_size = 0.1) + ylab(atlas[gg])
}
# Arrange plots
# Load 'gridExtra' package
if(!is.element("gridExtra",rownames(installed.packages()))){
  install.packages("gridExtra")
}
library(gridExtra)
g5 <- grid.arrange(grobs=gglist, nrow = 2, ncol=2)
# Save plot
outfile <- file.path(res_dir,"grp_all_forest.svg")
#ggsave(outfile, plot = g5, device = "svg", width = 3, height = 3)

# Compute edge-wise proportion test
if(file.exists(file.path(getwd(),"03-Inference","CMX_AAL_3D.rds"))){
  cmx <- readRDS(file.path(getwd(),"03-Inference","CMX_AAL_3D.rds"))
} else{
  # Generate connectivity matrices
  cmx <- array(0, dim=c(116,116,nrow(datos)))
  # Set same directory as TDA_ADHD200_NYU_pp2atlas.R 'ppdir'
  pp_dir <- 
    for(ii in 1:nrow(datos)){
      
      # Read time series files
      id <- as.character(datos$ID[ii])
      if(nchar(id)==5) id <- paste0("00",id)
      ts_file <- file.path(pp_dir,"NYU",id,"session_1",datos$Rest[ii],"rest.rsfMRIv2_nihpd","ts","pp_woGSR_AAL_ts.txt")
      # Compute connectivity matrix
      cmx[,,ii] <- cor(read.table(ts_file))
      
    }
}

# Apply threshold
cmx_dim <- dim(cmx)
cmx_th <- array(as.integer(cmx>=0.65), cmx_dim)
# Create empty objects
res <- array(0, cmx_dim[1:2])
resP <- array(1, cmx_dim[1:2])
# Apply proportion test
for(ii in 1:cmx_dim[1]){
  for(jj in 1:cmx_dim[2]){
    if(length(unique(cmx_th[ii,jj,]))>1){
      # Apply proportion test
      pt <- prop.test(table(datos$grupo,cmx_th[ii,jj,]))
      # Save results
      res[ii,jj] <- diff(pt$estimate)
      resP[ii,jj] <- pt$p.value
    }
  }
}
# Plot results
aal <- read.table(file.path(getwd(),"01-Preprocessing","atlas","AAL","ROI_MNI_V4.txt"))
colnames(res) <- aal$V3; rownames(res) <- aal$V3
# Use jet pallette
# Load 'matlab' package
if(!is.element("matlab",row.names(installed.packages()))){
  install.packages("matlab")
}
library(matlab)
par(op)
corrplot(res, col=rev(jet.colors(64)),is.corr = F, method = "square", mar=c(3,1,1,1)+0.1,
         cl.cex = 0.7, cl.lim = c(-0.3,0.3), tl.cex = 0.5, tl.col = "black")
# Only significant results
corrplot(res*(resP<0.01), col=rev(jet.colors(64)),is.corr = F, method = "square", mar=c(3,1,1,1)+0.1,
         cl.cex = 0.7, cl.lim = c(-0.3,0.3), tl.cex = 0.5, tl.col = "black")
# Generate connectivity matrix for BrainViewer
write.table((-1)*res*(resP<0.01), file.path(res_dir,"prop65_sig01_AAL.edge"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Generate chord diagram
# Load packages
if(sum(is.na(match(c("ggraph","igraph","tidyverse","RColorBrewer"),row.names(installed.packages()))))>0){
  install.packages(c("ggraph","igraph","tidyverse","RColorBrewer"))
}
# Libraries
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)

# create a data frame giving the hierarchical structure of your individuals
# Generate lobule category
aal$lob <- floor(aal$V3/1000)
aal$lob[which(aal$lob==3)] <- 4 # Add insula to paralimbic areas
aal$lob <- factor(aal$lob)
levels(aal$lob) <- c("FRONTAL","LIMBIC","OCCIPITAL","PARIETAL","SUBCORTICAL","TEMPORAL","CEREBELLUM")
# Create graph from adjacency matrix
ig <- graph.adjacency(res*(resP<0.01), mode="undirected", weighted=T)
# A coord diagram
ggraph(ig, layout = 'linear', circular = TRUE) + 
  theme_void() +
  geom_edge_arc(aes(colour = E(ig)$weight), edge_width = 2) +
  scale_edge_colour_gradientn(colours = rev(jet.colors(64)), space = "Lab", limits=c(-0.25,0.25),
                              na.value = "grey50", guide = "edge_colourbar") +
  geom_node_point(aes(x = x*1.07, y=y*1.07, colour=aal$lob), shape = 19, size=4) +
  scale_color_manual(values=c("blue3","brown3","chartreuse4","chocolate2",
                              "darkcyan","darkred","darkgoldenrod2"))
#ggsave(file.path(res_dir,"chord_prop65_sig01_AAL.svg"), device = "svg", width = 7.5, height = 6)

# Check edge-wise diagrams at different threshold distances
# Apply thresholds
for(hh in c(.5,.25,0)){
  
  # Apply threshold
  cmx_th <- array(as.integer(cmx>=hh), cmx_dim)
  # Create empty objects
  res <- array(0, cmx_dim[1:2])
  resP <- array(1, cmx_dim[1:2])
  # Apply proportion test
  for(ii in 1:cmx_dim[1]){
    for(jj in 1:cmx_dim[2]){
      if(length(unique(cmx_th[ii,jj,]))>1){
        # Apply proportion test
        pt <- prop.test(table(datos$grupo,cmx_th[ii,jj,]))
        # Save results
        res[ii,jj] <- diff(pt$estimate)
        resP[ii,jj] <- pt$p.value
      }
    }
  }
  
  # Generate connectivity matrix for BrainViewer
  #write.table((-1)*res*(resP<0.01), file.path(res_dir,paste0("prop",hh*100,"_sig01_AAL.edge")),
  #            quote = F, sep = "\t", row.names = F, col.names = F)
  
  # Create graph from adjacency matrix
  ig <- graph.adjacency(res*(resP<0.01), mode="undirected", weighted=T)
  # A coord diagram
  ggraph(ig, layout = 'linear', circular = TRUE) + 
    theme_void() +
    geom_edge_arc(aes(colour = E(ig)$weight), edge_width = 2) +
    scale_edge_colour_gradientn(colours = rev(jet.colors(64)), space = "Lab", limits=c(-0.35,0.35),
                                na.value = "grey50", guide = "edge_colourbar") +
    geom_node_point(aes(x = x*1.07, y=y*1.07, colour=aal$lob), shape = 19, size=5) +
    scale_color_manual(values=c("blue3","brown3","chartreuse4","chocolate2",
                                "darkcyan","darkred","darkgoldenrod2"))
  ggsave(file.path(res_dir,paste0("chord_prop",hh*100,"_sig01_AAL.pdf")), device = "pdf", width = 7.5, height = 6)
  #ggsave(file.path(res_dir,paste0("chord_prop",hh*100,"_sig001_AAL.svg")), device = "svg", width = 7.5, height = 6)
}

#------------------------------------------------------------------------------------------------------------------
# Check lobular level - logistic regression
#------------------------------------------------------------------------------------------------------------------

# Standardize variables
datos[,42:398] <- scale(datos[,42:398])

# Abbreviation
lob_lvl <- c("Fr","Limb","Occ","Par","Sub","Temp","Cbl")
lob_n <- length(lob_lvl)
tda <- c("area","kurt","slop")

# Logistic regression - Inter and inter module
# Create pairwise combinations
inter_comb <- combn(1:lob_n,2)
lob_lvl_comb <- sapply(1:ncol(inter_comb), function(x) paste0(lob_lvl[inter_comb[1,x]],".",lob_lvl[inter_comb[2,x]]))
lob_lvl_comb <- c(lob_lvl,lob_lvl_comb)
lob_comb_n <- length(lob_lvl_comb)
# Create empty matrix in order to store results
res <- array(rep(rep(c(1,0,1,1),each=lob_comb_n),3), dim=c(lob_comb_n,4*length(tda)))
# Row and column names
rownames(res) <- lob_lvl_comb
colnames(res) <- c(paste0(c("OR-","z-","p-","pNBS-"),rep(c("area","kurt","slop"),each=4)))
# Logistic model
for(ii in 1:lob_comb_n){
  logformula <- paste0("grupo ~ area",lob_lvl_comb[ii]," + kurt",lob_lvl_comb[ii],
                       " + slop",lob_lvl_comb[ii],
                       " + Sex + Age + AvgRelRMS")
  fit <- summary(glm(as.formula(logformula),data = datos, family=binomial("logit")))
  for(jj in 1:3){
    coefpos <- which(rownames(fit$coefficients)==paste0(tda[jj],lob_lvl_comb[ii]))
    if(length(coefpos)==1) res[ii,c(jj*4-3,jj*4-2,jj*4-1)] <- fit$coefficients[coefpos,c(1,3,4)]
  }
}
# Compute odds ratio
res[,4*(0:2)+1] <- exp(res[,4*(0:2)+1])

# Multiple comparison correction
# Network-based statistics (NBS)
# Define significant threshod (bi-sided)
# (check multiple comparison only in area because is the only with more than one significant edge)
zth <- abs(qnorm(0.025))
sig_links <- which(abs(res[,2])>zth)
mOBS <- matrix(0, nrow = lob_n, ncol = lob_n)
if(length(sig_links)>0){
  for(kk in sig_links){
    if(kk <= lob_n){
      mOBS[kk,kk] <- abs(res[kk,2])
    } else{
      mOBS[inter_comb[1,kk-lob_n],inter_comb[2,kk-lob_n]] <- abs(res[kk,2])
      mOBS[inter_comb[2,kk-lob_n],inter_comb[1,kk-lob_n]] <- abs(res[kk,2])
    }
  }
}

# Set 'igraph' object
gOBS <- graph.adjacency(mOBS, mode = "undirected", weighted = T)
# Permutate (nperm = 10000)
nperm <- 10000
# Empty objecto to store null distributions
null_d <- matrix(0, nrow = nperm, ncol = 2)
set.seed(18900217)
for(pp in 1:nperm){
  
  # Every 100 print message
  if((pp %% 100) == 0) cat(paste0(pp," "))
  
  # Create empty vector in order to store results
  resPERM <- vector("numeric", lob_comb_n)
  # Resample group labels
  datos$grpPERM <- sample(datos$grupo)
  # Logistic model
  for(ii in 1:lob_comb_n){
    logformula <- paste0("grpPERM ~ area",lob_lvl_comb[ii]," + kurt",lob_lvl_comb[ii],
                         " + slop",lob_lvl_comb[ii],
                         " + Sex + Age + AvgRelRMS")
    fit <- summary(glm(as.formula(logformula),data = datos, family=binomial("logit")))
    coefpos <- which(rownames(fit$coefficients)==paste0(tda[1],lob_lvl_comb[ii]))
    if(length(coefpos)==1) resPERM[ii] <- fit$coefficients[coefpos,3]
  }
  # Eval graph components
  sig_links <- which(abs(resPERM)>zth)
  if(length(sig_links)>0){
    mPERM <- matrix(0, nrow = lob_n, ncol = lob_n)
    for(kk in sig_links){
      if(kk <= lob_n){
        mPERM[kk,kk] <- abs(resPERM[kk])-zth
      } else{
        mPERM[inter_comb[1,kk-lob_n],inter_comb[2,kk-lob_n]] <- abs(resPERM[kk])-zth
        mPERM[inter_comb[2,kk-lob_n],inter_comb[1,kk-lob_n]] <- abs(resPERM[kk])-zth
      }
    }
    # Generate igraph object
    gPERM <- graph.adjacency(mPERM, mode = "undirected", weighted = T)
    gPERM_comp <- components(gPERM)
    maxPERM <- c(0,0)
    # Find components bigger than one
    if(sum(gPERM_comp$csize > 1)>0){
      mPERM[lower.tri(mPERM)] <- 0
      for(hh in 1:gPERM_comp$no){
        if(gPERM_comp$csize[hh]>1){
          sub_pos <- which(gPERM_comp$membership == hh)
          sub_mPERM <- mPERM[sub_pos,sub_pos]
          # Updated number of links maximum
          if(sum(c(sub_mPERM>0))>maxPERM[1]) maxPERM[1] <- sum(c(sub_mPERM>0))
          # Updated weighted number of links maximum
          if(sum(c(sub_mPERM))>maxPERM[2]) maxPERM[2] <- sum(c(sub_mPERM))
        } #if(gPERM_comp$csize[hh]>1)
      } #for(hh in 1:gPERM_comp$no)
      # Store maximum values
      null_d[pp,] <- maxPERM
    } #if(sum(gPERM_comp$csize > 1)>0)
  } #if(length(sig_links)>0){
}#for(pp in 1:nperm)

# Compute FWE p-value for observed data
gOBS_comp <- components(gOBS)
# Find components bigger than one
if(sum(gOBS_comp$csize > 1)>0){
  mOBS[lower.tri(mOBS)] <- 0
  maxOBS <- matrix(0, nrow = gOBS_comp$no, ncol = 4)
  for(hh in 1:gOBS_comp$no){
    if(gOBS_comp$csize[hh]>1){
      sub_pos <- which(gOBS_comp$membership == hh)
      sub_mOBS <- mOBS[sub_pos,sub_pos]
      # Store number of links (and FWE p-value)
      maxOBS[hh,1] <- sum(c(sub_mOBS>0))
      maxOBS[hh,3] <- sum(null_d[,1] >= maxOBS[hh,1])/nperm
      # Store weighted number of links
      maxOBS[hh,2] <- sum(c(sub_mOBS))-zth*maxOBS[hh,1]
      maxOBS[hh,4] <- sum(null_d[,2] >= maxOBS[hh,2])/nperm
      # Plot p-value stability
      pval_est <- cumsum(null_d[,2] >= maxOBS[hh,2])/(1:nperm)
      plot(pval_est, type="l", ylim = c(0,0.06), las = 1,
           xlab = "Permutation index", ylab = "p-value",
           main = paste0("Permuted p-value ",signif(maxOBS[hh,4],3)))
      abline(h=0.05, col="red", lty=2)
      # Marginal error
      pval_me <- 2*sqrt(pval_est*(1-pval_est)/1:nperm)
      lines(pval_est+pval_me, col = "chartreuse4")
      lines(pval_est-pval_me, col = "chartreuse4")
    } #if(gOBS_comp$csize[hh]>1)
    # Find FWE significant links in the sub-matrix
    comp_idx <- which(sub_mOBS>0, arr.ind = T)
    # Realocate indices to the entire matrix
    comp_idx[,1] <- sub_pos[comp_idx[,1]]
    comp_idx[,2] <- sub_pos[comp_idx[,2]]
    # Store FWE p-values
    for(gg in 1:nrow(comp_idx)){
      if(comp_idx[gg,1]==comp_idx[gg,2]){
        res[comp_idx[gg,1],4] <- maxOBS[hh,4]
      } else{
        res[intersect(
          which(inter_comb[1,]==comp_idx[gg,1]),
          which(inter_comb[2,]==comp_idx[gg,2])
        )+lob_n,4] <- maxOBS[hh,4]
      }
    } #for(gg in 1:length(comp_idx))
  } #for(hh in 1:gOBS_comp$no)
} #if(sum(gOBS_comp$csize > 1)>0)

# Find significant rows
sig_idx <- unique(which(res[,4*(1:3)]<=0.05,arr.ind = T)[,1])
signif(res[sig_idx,],4)
#resL <- res
#null_dL <- null_d

# Plot significant intra and inter-module effects (NBS corrected)
sig_mat <- matrix(0,nrow = lob_n, ncol = lob_n)
sig_diag <- which(res[1:lob_n,4]<=0.05)
diag(sig_mat)[sig_diag] <- log(res[sig_diag,1])
sig_tri <- which(res[-(1:lob_n),4]<=0.05)
for(ii in sig_tri) sig_mat[inter_comb[1,ii],inter_comb[2,ii]] <- log(res[ii+lob_n,1])
for(ii in sig_tri) sig_mat[inter_comb[2,ii],inter_comb[1,ii]] <- log(res[ii+lob_n,1])
rownames(sig_mat) <- c("FRT", "LIMB", "OCC", "PAR", "SUB", "TEM", "CBL")
colnames(sig_mat) <- c("FRT", "LIMB", "OCC", "PAR", "SUB", "TEM", "CBL")

# Corrplot
library(corrplot)
#pdf(file.path(res_dir,"cplot_AAL_lob.pdf"),4,4)
#svg(file.path(res_dir,"cplot_AAL_lob.svg"),4,4)
library(matlab)
corrplot(sig_mat, col=rev(jet.colors(100)),is.corr = F,type = "lower", tl.col = "black",
         method = "square", cl.cex = 0.7, cl.lim = c(-1,0), mar=c(3,1,1,1)+0.1)
#p.mat = 1-sig_fdr, insig = "pch", pch = "*", pch.col = "white")
title(xlab = "log(odds ratio) (pNBS < 0.05)")
#dev.off()

# Chord diagram of significant (uncorrected) intra and inter-module effects
# Load 'circlize' package
if(!is.element("circlize",rownames(installed.packages()))){
  install.packages("circlize")
}
library(circlize)

# Create network factor
net_fac <- as.factor(1:lob_n)
levels(net_fac) <- c("FRONTAL","LIMBIC","OCCIPITAL","PARIETAL","SUBCORT","TEMPORAL","CEREBEL")
# Initialize the plot.
pdf(file.path(res_dir,"chord_AAL_lob.pdf"),4,4)
#svg(file.path(res_dir,"chord_AAL_lob.svg"),4,4)
par(mar = c(1, 1, 1, 1) )
circos.clear()
circos.initialize(factors = net_fac, xlim=c(-1,1))
# Build the regions of track #1
circos.trackPlotRegion(factors = net_fac, ylim = c(-1,1),
                       bg.col = "aliceblue",
                       bg.border = "black")
# Add labels
fac_pos <- ux(34, "mm")*(1:lob_n)
circos.text(fac_pos,rep(0,lob_n),
            labels = levels(net_fac),
            facing = "bending.inside", cex=1.2)
# Add a links between a point and another
# Load 'scales' package
if(!is.element("scales",row.names(installed.packages()))){
  install.packages("scales")
}
library(scales)
# Find significant (NBS corrected) weigths
sig_tri <- which(res[-(1:lob_n),4]<=0.05)

# Draw links
for(ii in sig_tri){
  # Set central position
  p1 <- c(0.20,0.5); p2 <- c(-0.5,-0.2)
  # Set color
  t <- sig_mat[inter_comb[1,ii],inter_comb[2,ii]]
  tp <- round(100*(t+1)/2)
  tp_col <- rev(jet.colors(100))[tp]
  # Draw link
  circos.link(net_fac[inter_comb[1,ii]], p1,
              net_fac[inter_comb[2,ii]], p2,
              col = scales::alpha(tp_col,.9), border = scales::alpha(tp_col,.4),
              h.ratio=0.8)
}

# Draw intra-network above NBS corrected links
for(ii in sig_diag){
  
  # Set corner positions
  p1 <- c(-0.9,-0.6); p2 <- c(0.6,0.9)
  # Set color
  t <- sig_mat[ii,ii]
  tp <- round(100*(t+1)/2)
  tp_col <- rev(jet.colors(100))[tp]
  # Draw link
  circos.link(net_fac[ii], p1,
              net_fac[ii], p2,
              col = tp_col,
              #border = "magenta",
              lwd = 2,
              h.ratio=0.8)
  
}
dev.off()
# Clear cicle parameters
circos.clear()
par(op)

#------------------------------------------------------------------------------------------------------------------
# Check functional network level - logistic regression
#------------------------------------------------------------------------------------------------------------------

# Abbreviation
net_lvl <- c("AUD","CBL","CinOp","DMN","DAN","FPN","MEM","SAL","SMN.H","SMN.M","SUB","VAN","VIS")
net_n <- length(net_lvl)
tda <- c("area","kurt","slop")

# Logistic regression - Inter and inter module
# Create pairwise combinations
inter_comb <- combn(1:net_n,2)
net_lvl_comb <- sapply(1:ncol(inter_comb), function(x) paste0(net_lvl[inter_comb[1,x]],".",net_lvl[inter_comb[2,x]]))
net_lvl_comb <- c(net_lvl,net_lvl_comb)
net_comb_n <- length(net_lvl_comb)
# Create empty matrix in order to store results
res <- array(rep(rep(c(1,0,1,1),each=net_comb_n),3), dim=c(net_comb_n,4*length(tda)))
# Row and column names
rownames(res) <- net_lvl_comb
colnames(res) <- c(paste0(c("OR-","z-","p-","pNBS-"), rep(c("area","kurt","slop"),each=4)))
# Logistic model
for(ii in 1:net_comb_n){
  logformula <- paste0("grupo ~ area",net_lvl_comb[ii]," + kurt",net_lvl_comb[ii],
                       " + slop",net_lvl_comb[ii],
                       " + Sex + Age + AvgRelRMS")
  fit <- summary(glm(as.formula(logformula),data = datos, family=binomial("logit")))
  for(jj in 1:3){
    coefpos <- which(rownames(fit$coefficients)==paste0(tda[jj],net_lvl_comb[ii]))
    if(length(coefpos)==1) res[ii,c(jj*4-3,jj*4-2,jj*4-1)] <- fit$coefficients[coefpos,c(1,3,4)]
  }
}
# Compute odds ratio
res[,4*(0:2)+1] <- exp(res[,4*(0:2)+1])

# Multiple comparison correction
# Network-based statistics (NBS)
# Define significant threshod (bi-sided)
zth <- abs(qnorm(0.025))
sig_links <- which(abs(res[,2])>zth)
mOBS <- matrix(0, nrow = net_n, ncol = net_n)
if(length(sig_links)>0){
  for(kk in sig_links){
    if(kk <= net_n){
      mOBS[kk,kk] <- abs(res[kk,2])
    } else{
      mOBS[inter_comb[1,kk-net_n],inter_comb[2,kk-net_n]] <- abs(res[kk,2])
      mOBS[inter_comb[2,kk-net_n],inter_comb[1,kk-net_n]] <- abs(res[kk,2])
    }
  }
}
# Generate igraph object
gOBS <- graph.adjacency(mOBS, mode = "undirected", weighted = T)
# Permutate (nperm = 10000)
nperm <- 10000
# Empty objecto to store null distributions
null_d <- matrix(0, nrow = nperm, ncol = 2)
for(pp in 1:nperm){
  
  # Every 100 print message
  if((pp %% 100) == 0) cat(paste0(pp," "))
  
  # Create empty vector in order to store results
  resPERM <- vector("numeric", net_comb_n)
  # Resample group labels
  datos$grpPERM <- sample(datos$grupo)
  # Logistic model
  for(ii in 1:net_comb_n){
    logformula <- paste0("grpPERM ~ area",net_lvl_comb[ii]," + kurt",net_lvl_comb[ii],
                         " + slop",net_lvl_comb[ii],
                         " + Sex + Age + AvgRelRMS")
    fit <- summary(glm(as.formula(logformula),data = datos, family=binomial("logit")))
    coefpos <- which(rownames(fit$coefficients)==paste0(tda[1],net_lvl_comb[ii]))
    if(length(coefpos)==1) resPERM[ii] <- fit$coefficients[coefpos,3]
  }
  # Eval graph components
  sig_links <- which(abs(resPERM)>zth)
  if(length(sig_links)>0){
    mPERM <- matrix(0, nrow = net_n, ncol = net_n)
    for(kk in sig_links){
      if(kk <= net_n){
        mPERM[kk,kk] <- abs(resPERM[kk])-zth
      } else{
        mPERM[inter_comb[1,kk-net_n],inter_comb[2,kk-net_n]] <- abs(resPERM[kk])-zth
        mPERM[inter_comb[2,kk-net_n],inter_comb[1,kk-net_n]] <- abs(resPERM[kk])-zth
      }
    }
    # Generate igraph object
    gPERM <- graph.adjacency(mPERM, mode = "undirected", weighted = T)
    gPERM_comp <- components(gPERM)
    maxPERM <- c(0,0)
    # Find components bigger than one
    if(sum(gPERM_comp$csize > 1)>0){
      mPERM[lower.tri(mPERM)] <- 0
      for(hh in 1:gPERM_comp$no){
        if(gPERM_comp$csize[hh]>1){
          sub_pos <- which(gPERM_comp$membership == hh)
          sub_mPERM <- mPERM[sub_pos,sub_pos]
          # Updated number of links maximum
          if(sum(c(sub_mPERM>0))>maxPERM[1]) maxPERM[1] <- sum(c(sub_mPERM>0))
          # Updated weighted number of links maximum
          if(sum(c(sub_mPERM))>maxPERM[2]) maxPERM[2] <- sum(c(sub_mPERM))
        } #if(gPERM_comp$csize[hh]>1)
      } #for(hh in 1:gPERM_comp$no)
      # Store maximum values
      null_d[pp,] <- maxPERM
    } #if(sum(gPERM_comp$csize > 1)>0)
  } #if(length(sig_links)>0){
}#for(pp in 1:nperm)

# Compute FWE p-value for observed data
gOBS_comp <- components(gOBS)
# Find components bigger than one
if(sum(gOBS_comp$csize > 1)>0){
  mOBS[lower.tri(mOBS)] <- 0
  maxOBS <- matrix(0, nrow = gOBS_comp$no, ncol = 4)
  for(hh in 1:gOBS_comp$no){
    if(gOBS_comp$csize[hh]>1){
      sub_pos <- which(gOBS_comp$membership == hh)
      sub_mOBS <- mOBS[sub_pos,sub_pos]
      # Store number of links (and FWE p-value)
      maxOBS[hh,1] <- sum(c(sub_mOBS>0))
      maxOBS[hh,3] <- sum(null_d[,1] >= maxOBS[hh,1])/nperm
      # Store weighted number of links
      maxOBS[hh,2] <- sum(c(sub_mOBS))-zth*maxOBS[hh,1]
      maxOBS[hh,4] <- sum(null_d[,2] >= maxOBS[hh,2])/nperm
      # Plot p-value stability
      pval_est <- cumsum(null_d[,2] >= maxOBS[hh,2])/(1:nperm)
      plot(pval_est, type="l", ylim = c(0,0.06), las = 1,
           xlab = "Permutation index", ylab = "p-value",
           main = paste0("Permuted p-value ",signif(maxOBS[hh,4],3)))
      abline(h=0.05, col="red", lty=2)
      # Marginal error
      pval_me <- 2*sqrt(pval_est*(1-pval_est)/1:nperm)
      lines(pval_est+pval_me, col = "chartreuse4")
      lines(pval_est-pval_me, col = "chartreuse4")
    } #if(gOBS_comp$csize[hh]>1)
    # Find FWE significant links in the sub-matrix
    comp_idx <- which(sub_mOBS>0, arr.ind = T)
    # Realocate indices to the entire matrix
    comp_idx[,1] <- sub_pos[comp_idx[,1]]
    comp_idx[,2] <- sub_pos[comp_idx[,2]]
    # Store FWE p-values
    for(gg in 1:nrow(comp_idx)){
      if(comp_idx[gg,1]==comp_idx[gg,2]){
        res[comp_idx[gg,1],4] <- maxOBS[hh,4]
      } else{
        res[intersect(
          which(inter_comb[1,]==comp_idx[gg,1]),
          which(inter_comb[2,]==comp_idx[gg,2])
        )+net_n,4] <- maxOBS[hh,4]
      }
    } #for(gg in 1:length(comp_idx))
  } #for(hh in 1:gOBS_comp$no)
} #if(sum(gOBS_comp$csize > 1)>0)

# Find significant rows (NBS corrected)
sig_idx <- unique(which(res[,4*(1:3)]<=0.05,arr.ind = T)[,1])
signif(res[sig_idx,],4)
#resP <- res
#null_dN <- null_d

# Display results
# Plot significant (NBS corrected) intra and inter-module effects
sig_mat <- matrix(0,nrow = net_n, ncol = net_n)
sig_diag <- which(res[1:net_n,4]<=0.05)
diag(sig_mat)[sig_diag] <- log(res[sig_diag,1])
sig_tri <- which(res[-(1:net_n),4]<=0.05)
for(ii in sig_tri) sig_mat[inter_comb[1,ii],inter_comb[2,ii]] <- log(res[ii+net_n,1])
for(ii in sig_tri) sig_mat[inter_comb[2,ii],inter_comb[1,ii]] <- log(res[ii+net_n,1])
rownames(sig_mat) <- net_lvl; colnames(sig_mat) <- net_lvl
library(corrplot)
#pdf(file.path(res_dir,"cplot_P264_FN.pdf"),4,4)
#svg(file.path(res_dir,"cplot_P264_FN.svg"),4,4)
#par(mfrow=c(1,2))
library(matlab)
corrplot(sig_mat, col=rev(jet.colors(100)),is.corr = F,type = "lower", tl.col = "black",
         method = "square", cl.cex = 0.7, cl.lim = c(-1,0), mar=c(3,1,1,1)+0.1)
title(xlab = "log(OR) (p < 0.05; NBS-FWE)")
#dev.off()

# Chord diagram of significant (NBS corrected) intra and inter-module effects
# Load 'circlize' package
# Create network factor
net_fac <- as.factor(1:net_n)
levels(net_fac) <- net_lvl
# Initialize the plot.
#pdf(file.path(res_dir,"chord_P264_FN.pdf"),4,4)
#svg(file.path(res_dir,"chord_P264_FN.svg"),4,4)
par(mar = c(1, 1, 1, 1) )
circos.clear()
circos.initialize(factors = net_fac, xlim=c(-1,1))
# Build the regions of track #1
circos.trackPlotRegion(factors = net_fac, ylim = c(-1,1),
                       bg.col = "aliceblue",
                       bg.border = "black")
# Add labels
circos.text(2.24*(1:net_n),rep(0,net_n),
            labels = net_lvl,
            facing = "bending.inside", cex=1.2)
# Add a links between a point and another
library(scales)
# Find significant (uncorrected) weigths
#sig_z <- res[net_n+sig_tri,1]
# Draw links
for(ii in sig_tri){
  # Set random position
  p1 <- runif(1,-1,1)+c(-0.1,0.1); p2 <- runif(1,-1,1)+c(-0.1,0.1)
  #p1 <- c(-0.1,0.1); p2 <- c(-0.1,0.1)
  # Set color
  t <- sig_mat[inter_comb[1,ii],inter_comb[2,ii]]
  tp <- round(100*(t+1)/2)
  tp_col <- rev(jet.colors(100))[tp]
  # Draw link
  circos.link(net_fac[inter_comb[1,ii]], p1,
              net_fac[inter_comb[2,ii]], p2,
              col = scales::alpha(tp_col,.8),
              border = scales::alpha(tp_col,.5),
              h.ratio=0.8)
}
# Draw intra-network uncorrected links
for(ii in sig_diag){
  
  # Set corner positions
  p1 <- c(-0.9,-0.6); p2 <- c(0.6,0.9)
  # Set color
  t <- sig_mat[ii,ii]
  tp <- round(100*(t+1)/2)
  tp_col <- rev(jet.colors(100))[tp]
  # Draw link
  circos.link(net_fac[ii], p1,
              net_fac[ii], p2,
              col = scales::alpha(tp_col,.8),
              border = scales::alpha(tp_col,.5),
              h.ratio=0.8)
}
#dev.off()
# Clear cicle parameters
circos.clear()
par(op)
