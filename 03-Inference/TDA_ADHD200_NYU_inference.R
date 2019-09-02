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
# ADHD index
describeBy(datos$ADHD.Index, datos$grupo)
t.test(ADHD.Index~grupo, data=datos, var.equal = T)
# IQ
describeBy(datos$Full4.IQ, datos$grupo)
t.test(Full4.IQ~grupo, data=datos, var.equal = T)
# Average RMS Relative Displacement
describeBy(datos$AvgRelRMS, datos$grupo, digits = 4, mat = T)
t.test(AvgRelRMS~grupo, data=datos, var.equal = T)

#------------------------------------------------------------------------------------------------------------------
# Kendall Concordance Coefficient (KCC) along parcellations
#------------------------------------------------------------------------------------------------------------------

# Load 'irr' package
if(!is.element("irr",row.names(installed.packages()))){
  install.packages("irr")
}
library(irr)

# Betti0-AREA
# All atlases together
icc(cbind(scale(datos$areaAAL),scale(datos$areaCC200),scale(datos$areaP264),scale(datos$areaCC400)))
kendall(cbind(datos$areaAAL,datos$areaCC200,datos$areaP264,datos$areaCC400))

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
icc(cbind(scale(datos$kurtAAL),scale(datos$kurtCC200),scale(datos$kurtP264),scale(datos$kurtCC400)))
kendall(cbind(datos$kurtAAL,datos$kurtCC200,datos$kurtP264,datos$kurtCC400))

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
icc(cbind(scale(datos$slopAAL),scale(datos$slopCC200),scale(datos$slopP264),scale(datos$slopCC400)))
kendall(cbind(datos$slopAAL,datos$slopCC200,datos$slopP264,datos$slopCC400))

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

# Plot results
# Load 'corrplot' package
if(!is.element("corrplot",row.names(installed.packages()))){
  install.packages("corrplot")
}
library(corrplot)

# Set a results directory in order to save images
res_dir <- ""
if(!dir.exists(res_dir)) dir.create(res_dir)
#pdf(file.path(res_dir,"KCC.pdf"),5,5)
#svg(file.path(res_dir,"KCC.svg"),5,5)
par(mfrow=c(1,3))
# Area
colnames(area_kcc) <- c("AAL","CC200","P264","CC400")
rownames(area_kcc) <- c("AAL","CC200","P264","CC400")
colorcito <- c(heat.colors(50),rev(heat.colors(50)))
corrplot(area_kcc,method = "square", type = "upper", diag=F, mar=c(2,0,2,0),
         col = colorcito, tl.col = "black", addCoef.col = "black",
         tl.srt = 45, is.corr = F, cl.lim = c(0,1), cl.ratio = .3)
# Kurtosis
colnames(kurt_kcc) <- c("AAL","CC200","P264","CC400")
rownames(kurt_kcc) <- c("AAL","CC200","P264","CC400")
corrplot(kurt_kcc,method = "square", type = "upper", diag=F, mar=c(2,0,2,0),
         col = colorcito, tl.col = "black", addCoef.col = "black",
         tl.srt = 45, is.corr = F, cl.lim = c(0,1), cl.ratio = .3)
# Slope
colnames(slop_kcc) <- c("AAL","CC200","P264","CC400")
rownames(slop_kcc) <- c("AAL","CC200","P264","CC400")
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
# Groups comparison - single simplex filtration values
#------------------------------------------------------------------------------------------------------------------

# Read Rips filtrations for each participant
fval <- matrix(0, nrow = 2, ncol = 4)
atlas <- c("AAL","CC200","P264","CC400")
for(ii in 1:length(atlas)){
  
  # Read file
  filename <- file.path(getwd(),"02-TDA",paste0("ADHD200_NYU_ppNIHPD_Rips_",atlas[ii],".csv"))
  rips <- read.csv(filename)
  
  # Match subjects
  rips <- rips[match(datos$ID,rips$ID),]
  
  # Describe by group
  show(db <- describeBy(rips$X1, datos$grupo, mat = T, digits = 4))
  fval[,ii] <- db$mean
  show(t.test(rips$X1 ~ datos$grupo, var.equal = T))
  
  # Compute Cohen's d
  show(d <- (db$mean[1] - db$mean[2])/sd(rips$X1))
  
}
colnames(fval) <- atlas
rownames(fval) <- c("TDC", "ADHD")
show(fval)

#------------------------------------------------------------------------------------------------------------------
# Groups comparison - logistic regression
#------------------------------------------------------------------------------------------------------------------

# Create empty object to save every result
res <- array(0, dim=c(12,6))
rownames(res) <- c("AAL-OR","AAL-z","AAL-p",
                   "CC200-OR","CC200-z","CC200-p",
                   "P264-OR","P264-z","P264-p",
                   "CC400-OR","CC400-z","CC400-p")
colnames(res) <- c("area","kurt","slop","sex","age","avgFD")

# AAL
fit <- summary(glm(grupo ~ scale(areaAAL) + scale(kurtAAL) + scale(slopAAL) + Sex + scale(Age) + scale(AvgRelRMS),
                   data = datos, family=binomial("logit")))
res[1,] <- exp(fit$coefficients[2:7,1])
res[2,] <- fit$coefficients[2:7,3]
res[3,] <- fit$coefficients[2:7,4]
# CC200
fit <- summary(glm(grupo ~ scale(areaCC200) + scale(kurtCC200) + scale(slopCC200) + Sex + scale(Age) + scale(AvgRelRMS),
                   data = datos, family=binomial("logit")))
res[4,] <- exp(fit$coefficients[2:7,1])
res[5,] <- fit$coefficients[2:7,3]
res[6,] <- fit$coefficients[2:7,4]
# P264
fit <- summary(glm(grupo ~ scale(areaP264) + scale(kurtP264) + scale(slopP264) + Sex + scale(Age) + scale(AvgRelRMS),
                   data = datos, family=binomial("logit")))
res[7,] <- exp(fit$coefficients[2:7,1])
res[8,] <- fit$coefficients[2:7,3]
res[9,] <- fit$coefficients[2:7,4]
# CC400
fit <- summary(glm(grupo ~ scale(areaCC400) + scale(kurtCC400) + scale(slopCC400) + Sex + scale(Age) + scale(AvgRelRMS),
                   data = datos, family=binomial("logit")))
res[10,] <- exp(fit$coefficients[2:7,1])
res[11,] <- fit$coefficients[2:7,3]
res[12,] <- fit$coefficients[2:7,4]

# Results
signif(res,3)

# Print results (AAL atlas)
#pdf(file.path(res_dir,"grp_all_rips.pdf"),7,5)
#svg(file.path(res_dir,"grp_all_rips.svg"),7,5)
par(mfrow = c(2,2))

for(ii in 1:length(atlas)){
  
  # Filepath
  filename <- file.path(getwd(),"02-TDA",paste0("ADHD200_NYU_ppNIHPD_Rips_",atlas[ii],".csv"))
  
  # Read Rips filtration for each participant
  rips <- read.csv(filename)
  # Match IDs and rest scan
  rips <- rips[match(datos$ID,rips$ID),]
  rips <- rips[,-(1:3)]
  # Get group averages
  rips$Y <- rep(0,nrow(rips)) # Add distance zero
  ripsAVG <- describeBy(rips, datos$grupo, mat = T)
  # (Empty) plot
  y <- 1:ncol(rips); x <- y/length(y)
  plot(x,y, type = "n", xlim = c(0,0.6), las = 1, frame.plot = F,
       ylab = "Betti 0", xlab = "Filtration Value", main = atlas[ii])
  axis(side = 1, lwd = 2)
  axis(side = 2, lwd = 2, las = 1)
  # Draw c.i. light lines
  colorcito <- c("blue", "red")
  # Draw group means and dashed confidence intervals
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
  # Add filtration value for betti0 equal to 1
  # TDC
  lines(x = c(fval[1,ii],fval[1,ii]), y = c(1, ncol(rips)/5), col = "blue", lwd = 2)
  # ADHD
  lines(x = c(fval[2,ii],fval[2,ii]), y = c(1, ncol(rips)/5), col = "red", lwd = 2)
  
}
  
#dev.off()
par(op)

# Add gruped-boxplots
# Load 'ggplot2' package
if(!is.element("ggplot2",row.names(installed.packages()))){
  install.packages("ggplot2")
}
library(ggplot2)
theme_set(theme_classic())
g1 <- ggplot(datos, aes(x=grupo, y=scale(areaAAL), fill=grupo, colour=grupo))+
  geom_boxplot(alpha = .5, width = .2, notch=TRUE)+
  scale_colour_manual(values=c("blue", "red"))+
  scale_fill_manual(values=c("blue", "red"))+
  ylab("AREA (z)")
g2 <- ggplot(datos, aes(x=grupo, y=scale(areaCC200), fill=grupo, colour=grupo))+
  geom_boxplot(alpha = .5, width = .2, notch=TRUE)+
  scale_colour_manual(values=c("blue", "red"))+
  scale_fill_manual(values=c("blue", "red"))+
  ylab("AREA (z)")
g3 <- ggplot(datos, aes(x=grupo, y=scale(areaP264), fill=grupo, colour=grupo))+
  geom_boxplot(alpha = .5, width = .2, notch=TRUE)+
  scale_colour_manual(values=c("blue", "red"))+
  scale_fill_manual(values=c("blue", "red"))+
  ylab("AREA (z)")
g4 <- ggplot(datos, aes(x=grupo, y=scale(areaCC400), fill=grupo, colour=grupo))+
  geom_boxplot(alpha = .5, width = .2, notch=TRUE)+
  scale_colour_manual(values=c("blue", "red"))+
  scale_fill_manual(values=c("blue", "red"))+
  ylab("AREA (z)")
# Combine plots
# Load 'gridExtra' package
if(!is.element("gridExtra",rownames(installed.packages()))){
  install.packages("gridExtra")
}
library(gridExtra)
g5 <- grid.arrange(g1, g2, g3, g4, nrow = 2, ncol=2)
# Save plot
outfile <- file.path(res_dir,"grp_all_box.svg")
#ggsave(outfile, plot = g5, device = "svg")

# Add Cohen's d
areaAVG <- describeBy(datos$areaAAL, datos$grupo)
mu <- mean(datos$areaAAL[which(datos$grupo=="TDC")]) - mean(datos$areaAAL[which(datos$grupo=="ADHD")])
pools <- sqrt(
  (((sum(datos$grupo=="TDC")-1) * var(datos$areaAAL[which(datos$grupo=="TDC")])) +
    (sum(datos$grupo=="ADHD")-1)*var(datos$areaAAL[which(datos$grupo=="ADHD")]))/
    (length(datos$grupo)-2))
d <- mu/pools

# Compute edge-wise proportion test
# Generate connectivity matrices
cmx <- array(0, dim=c(116,116,nrow(datos)))
# Set same directory as TDA_ADHD200_NYU_pp.R 'outdir'
pp_dir <- ""
for(ii in 1:nrow(datos)){
  
  # Read time series files
  id <- as.character(datos$ID[ii])
  if(nchar(id)==5) id <- paste0("00",id)
  ts_file <- file.path(pp_dir,"NYU",id,"session_1",datos$Rest[ii],"rest.rsfMRIv2_nihpd","ts","pp_woGSR_AAL_ts.txt")
  # Compute connectivity matrix
  cmx[,,ii] <- cor(read.table(ts_file))
  
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
par(op)
library(matlab)
corrplot(res, col=rev(jet.colors(64)),is.corr = F, method = "square", mar=c(3,1,1,1)+0.1,
         cl.cex = 0.7, cl.lim = c(-0.3,0.3), tl.cex = 0.5, tl.col = "black")
# Only significant results
corrplot(res*(resP<0.01), col=rev(jet.colors(64)),is.corr = F, method = "square", mar=c(3,1,1,1)+0.1,
         cl.cex = 0.7, cl.lim = c(-0.3,0.3), tl.cex = 0.5, tl.col = "black")
# Generate connectivity matrix for BrainViewer
write.table((-1)*res*(resP<0.01), file.path(res_dir,"prop65_sig01_AAL.edge"),
            quote = F, sep = "\t", row.names = F, col.names = F)

# Generate chord diagram
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
  ig <- graph.adjacency(res*(resP<0.001), mode="undirected", weighted=T)
  # A coord diagram
  ggraph(ig, layout = 'linear', circular = TRUE) + 
    theme_void() +
    geom_edge_arc(aes(colour = E(ig)$weight), edge_width = 2) +
    scale_edge_colour_gradientn(colours = rev(jet.colors(64)), space = "Lab", limits=c(-0.35,0.35),
                                na.value = "grey50", guide = "edge_colourbar") +
    geom_node_point(aes(x = x*1.07, y=y*1.07, colour=aal$lob), shape = 19, size=5) +
    scale_color_manual(values=c("blue3","brown3","chartreuse4","chocolate2",
                                "darkcyan","darkred","darkgoldenrod2"))
  ggsave(file.path(res_dir,paste0("chord_prop",hh*100,"_sig001_AAL.svg")), device = "svg", width = 7.5, height = 6)
  
}

#------------------------------------------------------------------------------------------------------------------
# Check lobular level - logistic regression
#------------------------------------------------------------------------------------------------------------------

# Standardzise variables
datos[,42:398] <- scale(datos[,42:398])

# Abbreviation
lob_lvl <- c("Fr","Limb","Occ","Par","Sub","Temp","Cbl")
lob_n <- length(lob_lvl)
tda <- c("area","kurt","slop")

# Logistic regression - Intra-modular
# Create empty object to save every result
res <- array(rep(rep(c(0,1,1),each=lob_n),3), dim=c(lob_n,9))
rownames(res) <- lob_lvl
colnames(res) <- c(paste0(c("z-","p-","pFDR-"),rep(c("area","kurt","slop"),each=3)))
# Logistic model
for(ii in 1:lob_n){
  logformula <- paste0("grupo ~ area",lob_lvl[ii]," + kurt",lob_lvl[ii]," + slop",lob_lvl[ii],
                       " + Sex + Age + AvgRelRMS")
  fit <- summary(glm(as.formula(logformula),data = datos, family=binomial("logit")))
  for(jj in 1:3){
    coefpos <- which(rownames(fit$coefficients)==paste0(tda[jj],lob_lvl[ii]))
    if(length(coefpos)==1) res[ii,c(jj*3-2,jj*3-1)] <- fit$coefficients[coefpos,3:4]
  }
}
# Multiple comparison correction
res[,c(3,6,9)] <- cbind(p.adjust(res[,2],"fdr"),
                        p.adjust(res[,5],"fdr"),
                        p.adjust(res[,8],"fdr"))

# Logistic regression - Inter and inter module
# Create pairwise combinations
inter_comb <- combn(1:lob_n,2)
lob_lvl_comb <- sapply(1:ncol(inter_comb), function(x) paste0(lob_lvl[inter_comb[1,x]],".",lob_lvl[inter_comb[2,x]]))
lob_lvl_comb <- c(lob_lvl,lob_lvl_comb)
lob_comb_n <- length(lob_lvl_comb)
# Area
res <- array(rep(rep(c(1,0,1,1),each=lob_comb_n),3), dim=c(lob_comb_n,12))
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
res[,c(1,5,9)] <- exp(res[,c(1,5,9)])
# Multiple comparison correction
res[,c(4,8,12)] <- cbind(p.adjust(res[,3],"fdr"),
                         p.adjust(res[,7],"fdr"),
                         p.adjust(res[,11],"fdr"))
rownames(res) <- lob_lvl_comb
colnames(res) <- c("area-OR","area-z","area-p","area-pFDR",
                   "kurt-OR","kurt-z","kurt-p","kurt-pFDR",
                   "slop-OR","slop-z","slop-p","slop-pFDR")
#write.csv(res,file.path(res_dir,"logit_module.csv"))
# Find significant rows
sig_idx <- unique(which(res[,c(4,8,12)]<=0.05,arr.ind = T)[,1])
signif(res[sig_idx,],4)
resL <- res

# Plot significant intra and inter-module effects
sig_mat <- matrix(0,nrow = lob_n, ncol = lob_n)
sig_diag <- which(res[1:lob_n,3]<=0.05)
diag(sig_mat)[sig_diag] <- log(res[sig_diag,1])
sig_tri <- which(res[-(1:lob_n),3]<=0.05)
for(ii in sig_tri) sig_mat[inter_comb[1,ii],inter_comb[2,ii]] <- log(res[ii+lob_n,1])
for(ii in sig_tri) sig_mat[inter_comb[2,ii],inter_comb[1,ii]] <- log(res[ii+lob_n,1])
rownames(sig_mat) <- c("FRT", "LIMB", "OCC", "PAR", "SUB", "TEM", "CBL")
colnames(sig_mat) <- c("FRT", "LIMB", "OCC", "PAR", "SUB", "TEM", "CBL")
# Get FDR corrected
sig_fdr <- matrix(1,nrow = lob_n, ncol = lob_n)
sig_diag <- which(res[1:lob_n,4]<=0.05)
diag(sig_fdr)[sig_diag] <- res[sig_diag,4]
sig_tri <- which(res[-(1:lob_n),4]<=0.05)
for(ii in sig_tri) sig_fdr[inter_comb[1,ii],inter_comb[2,ii]] <- res[ii+lob_n,4]
for(ii in sig_tri) sig_fdr[inter_comb[2,ii],inter_comb[1,ii]] <- res[ii+lob_n,4]
rownames(sig_fdr) <- lob_lvl; colnames(sig_fdr) <- lob_lvl

par(op)
# Corrplot
# Load 'corrplot' package
if(!is.element("corrplot",rownames(installed.packages()))){
  install.packages("corrplot")  
}
library(corrplot)
# Load 'matlab' package
if(!is.element("matlab",rownames(installed.packages()))){
  install.packages("matlab")  
}
library(matlab)
#pdf(file.path(res_dir,"cplot_AAL_lob.pdf"),4,4)
#svg(file.path(res_dir,"cplot_AAL_lob.svg"),4,4)
corrplot(sig_mat, col=rev(jet.colors(100)),is.corr = F,type = "lower", tl.col = "black",
         method = "square", cl.cex = 0.7, cl.lim = c(-1,0), mar=c(3,1,1,1)+0.1,
         p.mat = 1-sig_fdr, insig = "pch", pch = "*", pch.col = "white")
title(xlab = "log(odds ratio) (p < 0.05)")
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
if(!is.element("scales",rownames(installed.packages()))){
  install.packages("scales")
}
library(scales)
# Find significant (uncorrected) weigths
sig_tri <- which(res[-(1:lob_n),3]<=0.05)
sig_fdr <- which(res[-(1:lob_n),4]<=0.05)
# Draw links
for(ii in sig_tri){
  # If uncorrected, do not draw borders
  if(!is.element(ii,sig_fdr)){
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
}
# Draw above FDR corrected links
for(ii in sig_fdr){

  # Set random position
  p1 <- runif(1,-0.9,0.9)+c(-0.15,0.15)#; p2 <- runif(1,-0.9,0.9)+c(-0.15,0.15)
  p2 <- c(-0.15,0.15)
  #p1 <- c(-0.25,0.25); p2 <- c(-0.25,0.25)
  # Set color
  t <- sig_mat[inter_comb[1,ii],inter_comb[2,ii]]
  tp <- round(100*(t+1)/2)
  tp_col <- rev(jet.colors(100))[tp]
  # Draw link
  circos.link(net_fac[inter_comb[1,ii]], p1,
              net_fac[inter_comb[2,ii]], p2,
              col = tp_col,
              #border = "magenta",
              lwd = 2,
              h.ratio=0.8)
  
}
# Draw intra-network above FDR corrected links
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
#dev.off()
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

# Logistic regression - Intra-modular
# Create empty object to save every result
res <- array(rep(rep(c(1,0,1,1),each=net_n),3), dim=c(net_n,12))
rownames(res) <- net_lvl
colnames(res) <- c(paste0(c("OR-","z-","p-","pFDR-"),rep(c("area","kurt","slop"),each=4)))
# Logistic model
for(ii in 1:net_n){
  logformula <- paste0("grupo ~ area",net_lvl[ii]," + kurt",net_lvl[ii]," + slop",net_lvl[ii],
                       " + Sex + Age + AvgRelRMS")
  fit <- summary(glm(as.formula(logformula),data = datos, family=binomial("logit")))
  for(jj in 1:3){
    coefpos <- which(rownames(fit$coefficients)==paste0(tda[jj],net_lvl[ii]))
    if(length(coefpos)==1) res[ii,(jj*4-3):(jj*4-1)] <- fit$coefficients[coefpos,c(1,3:4)]
  }
}
# Compute odds ratio
res[,c(1,5,9)] <- exp(res[,c(1,5,9)])
# Multiple comparison correction
res[,c(4,8,12)] <- cbind(p.adjust(res[,3],"fdr"),
                        p.adjust(res[,7],"fdr"),
                        p.adjust(res[,11],"fdr"))

# Logistic regression - Inter and inter module
# Create pairwise combinations
inter_comb <- combn(1:net_n,2)
net_lvl_comb <- sapply(1:ncol(inter_comb), function(x) paste0(net_lvl[inter_comb[1,x]],".",net_lvl[inter_comb[2,x]]))
net_lvl_comb <- c(net_lvl,net_lvl_comb)
net_comb_n <- length(net_lvl_comb)
# Area
res <- array(rep(rep(c(1,0,1,1),each=net_comb_n),3), dim=c(net_comb_n,12))
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
res[,c(1,5,9)] <- exp(res[,c(1,5,9)])
# Multiple comparison correction
res[,c(4,8,12)] <- cbind(p.adjust(res[,3],"fdr"),
                         p.adjust(res[,7],"fdr"),
                         p.adjust(res[,11],"fdr"))
rownames(res) <- net_lvl_comb
colnames(res) <- c("area-OR","area-z","area-p","area-pFDR",
                   "kurt-OR","kurt-z","kurt-p","kurt-pFDR",
                   "slop-OR","slop-z","slop-p","slop-pFDR")
#write.csv(res, file.path(res_dir,"logit_module.csv"))
# Find significant rows
#sig_idx <- unique(which(res[,c(4,8,12)]<=0.05,arr.ind = T)[,1]) #FDR
sig_idx <- unique(which(res[,c(3,8,12)]<=0.05,arr.ind = T)[,1]) #uncorr
signif(res[sig_idx,],4)
resP <- res

# Display results
# Plot significant (uncorrected) intra and inter-module effects
sig_mat <- matrix(0,nrow = net_n, ncol = net_n)
sig_diag <- which(res[1:net_n,3]<=0.05)
diag(sig_mat)[sig_diag] <- log(res[sig_diag,1])
sig_tri <- which(res[-(1:net_n),3]<=0.05)
for(ii in sig_tri) sig_mat[inter_comb[1,ii],inter_comb[2,ii]] <- log(res[ii+net_n,1])
for(ii in sig_tri) sig_mat[inter_comb[2,ii],inter_comb[1,ii]] <- log(res[ii+net_n,1])
rownames(sig_mat) <- net_lvl; colnames(sig_mat) <- net_lvl
library(corrplot)
#colorcito <- colorRampPalette(c("cyan","blue","white","white","white","white"))
#pdf(file.path(res_dir,"cplot_P264_FN.pdf"),4,4)
#svg(file.path(res_dir,"cplot_P264_FN.svg"),4,4)
#par(mfrow=c(1,2))
library(matlab)
corrplot(sig_mat, col=rev(jet.colors(100)),is.corr = F,type = "lower", tl.col = "black",
         method = "square", cl.cex = 0.7, cl.lim = c(-1,0), mar=c(3,1,1,1)+0.1)
title(xlab = "log(OR) (p < 0.05; uncorr.)")
#dev.off()

# Chord diagram of significant (uncorrected) intra and inter-module effects
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
# Draw intra-network above significant links
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
