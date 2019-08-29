#!/usr/bin/Rscript

# Esta función calcula las derivadas de las columnas, y el cuadrado de las columnas y las derivadas
# uso: Rscript dp2.R input output

filename <- commandArgs(TRUE)
a1=filename[1];
a2=filename[2];
dp2 <- function(a1,a2) {
  #filename='prefiltered_func_data_mcf.par'
  datos <- read.table(a1)
  
  datos <- cbind(datos,apply(datos,2,function(x) c(0,diff(x))))
  datos <- cbind(datos,apply(datos,2,function(x) x^2))

  # Exporta resultados en otro archivo
  write.table(datos,a2,row.names=F,col.names=F,sep='\t')
}

dp2(a1,a2)