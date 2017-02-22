
######################### Load the Packages ########################
library(optimx);


######################### Define the Simplified Negative Log Likelihood Function ########################
NLL <- function(theta){
  NLL.value <- log(theta[1]) + 1/n*sum(log(theta[2]/p*NLL.Lambda+1)) + 
               (1/n/theta[1])*sum((NLL.yTilde.sq)/(theta[2]/p*NLL.Lambda+1));
  return(as.numeric(NLL.value));
};
cat("Before computing NLL, do the following preparations:\n");
cat("Specify and save NLL.X, NLL.y, NLL.U, NLL.Lambda, NLL.yTilde and NLL.yTilde.sq");


######################### Define the Score and Hessian Functions ########################
Gradient.NLL <- function(theta) {
  S1 <- 1/theta[1] - 1/(n*(theta[1])^2) * sum((NLL.yTilde.sq*p)/(theta[2]*NLL.Lambda+p));
  S2 <- 1/n * sum(NLL.Lambda/(theta[2]*NLL.Lambda+p)) - 1/(n*theta[1]) * sum((NLL.yTilde.sq*NLL.Lambda*p)/((theta[2]*NLL.Lambda+p)^2));
  return(matrix(c(S1,S2), nrow=2, ncol=1));
};

Hessian.NLL <- function(theta) {
  H11 <- -1/(theta[1]^2) + 2/(n*theta[1]^3) * sum((NLL.yTilde.sq*p)/(theta[2]*NLL.Lambda+p));
  H12 <- 1/(n*theta[1]^2) * sum((NLL.yTilde.sq*NLL.Lambda*p)/((theta[2]*NLL.Lambda+p)^2));
  H22 <- -1/n * sum((NLL.Lambda^2)/((theta[2]*NLL.Lambda + p)^2)) + 2/(n*theta[1]) * sum((NLL.yTilde.sq*(NLL.Lambda^2)*p)/((theta[2]*NLL.Lambda+p)^3));
  return(matrix(c(H11,H12,H12,H22), nrow=2, ncol=2));
};