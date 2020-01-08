#!/usr/bin/Rscript

# This function generates a covariate column for each one in a vector.

filename <- commandArgs(TRUE)
a1=filename[1];
a2=filename[2];
colpeaks <- function(a1,a2) {
  picos <- read.table(a1)
  
  if(sum(picos)>0){
    idx_picos <- which(picos == 1)
    mat_picos <- matrix(0,nrow(picos),length(idx_picos))
    for(ii in 1:length(idx_picos)) mat_picos[idx_picos[ii],ii] <- 1
    write.table(mat_picos,a2,row.names=F,col.names=F,sep='\t')
  }

}
colpeaks(a1,a2)
