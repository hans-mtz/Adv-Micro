# Clear workspace
rm(list = ls())

# Check to see if packages are installed
list.of.packages <- c("rcdd", "R.matlab")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load libraries
library(rcdd)
library(R.matlab)

# Load matrix
data <- readMat("A_Vrepresentation.mat")

# Set to matrix
A <- data[["A"]]

# Transpose 
A  <- t(A)

# Size of matrixBB
H <- nrow(A)

# Include indicator columns for A1
A   <-  cbind( matrix(rep(0,H)), matrix(rep(0,H)),A)

# Tell R that A is V representation
attr(A, 'representation') = "V" 

# Return H representation
B   <- scdd(A)

# Remove attributes
attr(B$output, "representation") <- NULL

# Save vector indicating equalities
Beq <- B$output[,c(1)]
Beq <- matrix(Beq)

# Delete first two columns
A  <-  A[,-c(1,2)]
B  <-  B$output[,-c(1,2)]
B  <- -B

# Write to a .mat
writeMat("B_Hrepresentation.mat",A=t(A),B=B,Beq=Beq)
