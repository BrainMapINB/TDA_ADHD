# Load TDA package
if(!is.element("TDA",row.names(installed.packages()))){
  # Likely, this libraries will be needed: libgmp3-dev, 
  install.packages("TDA")
}
library(TDA)

# Define function
rips_hom <- function(matriz){
  # Compute distance matrix
  D=1-matriz
  # Take number of nodes
  l=dim(D)[1]
  # Maximum value of the rips filtration
  mscale <- 1
  # Maximum Betti number
  mdim <- 1
  DiagTri <- ripsDiag(D, mdim, mscale, dist = "arbitrary", printProgress = F)
  # Save Rips filtration
  betti0 <- which(DiagTri$diagram[,1]==0)
  muerte <- DiagTri$diagram[betti0[-1],3]
  
  # Return results
  return(muerte)
}