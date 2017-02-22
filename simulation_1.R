############ Step 0: Prelimiaries ############
## 0-1 Load package
source("Load_NLL.R");

## 0-2 Specify the parameters
## Basic parameters
n <- 300;  #____;
ncol.Fixed <- 1500;  #____;
ncol.Perm <- ncol.Fixed - 1;
p <- ncol.Fixed + ncol.Perm;
tauSq.Fixed <- 10;  #____;
tauSq.Perm <- 0; #____;
sigmaSq <- 2; #____;
etaSq <- (tauSq.Fixed+tauSq.Perm) / sigmaSq;
alpha <- 100; # (1-alpha)% of the coefficients set at zero;
cat("True sigmaSq is ", sigmaSq, ";\n", sep="");
cat("True etaSq is ", etaSq, ";\n", sep="");

## Simulation parameters
seed.base <- 1000; #____;
num.dataset <- 100;  #____;
d.start <- 1;
num.permutation <- 500;  #____;
d.end <- d.start + num.dataset - 1;

## Determine the predictor coefficients
Beta <- rep(sqrt(tauSq.Fixed / ncol.Fixed), ncol.Fixed);
Gamma <- rep(sqrt(tauSq.Perm / ncol.Perm), ncol.Perm);

## 0-3 Specify the output dataset
Output.True <- data.frame(Theta1=NA, Theta2=NA, value=NA, fevals=NA, gevals=NA,  niter=NA, convcode=NA, kkt1=NA, kkt2=NA, xtimes=NA);
Output.Perm <- data.frame(SEED=NA, Theta1=NA, Theta2=NA, convcode=NA);
Output.PVal <- data.frame(SEED=NA, Empirical.p=NA);
q <- 0;

for (d in d.start:d.end) {
############ Step 1: Dataset Preparation ############
tempus <- Sys.time();
seed.num <- seed.base + d;
set.seed(seed.num);

## 1-1 Randomly generate predictor variables with indepdent normal
X <- matrix(rnorm(n=n*ncol.Fixed), nrow=n);
Z <- matrix(rnorm(n=n*ncol.Perm), nrow=n);  #____;
## 1-2 Randomely generate random error from iid normal distribution
Error <- rnorm(n, mean=0, sd=sigmaSq^0.5);
## 1-3 Deterministically generate the y response
y <- X %*% Beta + Z %*% Gamma + Error;
## 1-4 Deterministically write y, X and Z into a dataset
Data1 <- data.frame(y,X,Z);
colnames(Data1) <- c("y", paste("X", 1:ncol.Fixed, sep=""), paste("Z", 1:ncol.Perm, sep=""));

set.seed(NULL);

############ Step 2: MLE Estimation of Variance and Signal to Noise Ratio ############
## 2-1 Prespecify theta near the true parameters
theta0 <- c(sigmaSq+runif(1,-1,1), etaSq+runif(1,-1,1));

## 2-2 Process the dataset in for function NLL()
NLL.y <- Data1[ ,1];
NLL.XZ <- as.matrix(Data1[ ,-1]);
NLL.DC <- eigen(NLL.XZ %*% t(NLL.XZ));
NLL.U <- NLL.DC$vectors;
NLL.Lambda <- NLL.DC$values;
NLL.yTilde <- t(NLL.U) %*% NLL.y;
NLL.yTilde.sq <- NLL.yTilde^2;
## 2-3 Output the estimated sigmaSq and etaSq
True.NLL <- optimx(par=theta0, fn=NLL, gr=Gradient.NLL, hess=Hessian.NLL, method="BFGS");
Output.True[d, ] <- True.NLL[1, ];
rownames(Output.True)[d] <- seed.num;

############ Step 3: Permutation Based Test ############
for (pm in 1:num.permutation) {
  q <- q+1;
  ## 3-1 Prepare the permutation matrix and vectors
  NLL.y <- y;
  Z.Perm <- Z[sample(nrow(Z)), ];
  NLL.XZ <- cbind(X, Z.Perm);
  NLL.DC <- eigen(NLL.XZ %*% t(NLL.XZ));
  NLL.U <- NLL.DC$vectors;
  NLL.Lambda <- NLL.DC$values;
  NLL.yTilde <- t(NLL.U) %*% NLL.y;
  NLL.yTilde.sq <- NLL.yTilde^2;
  ## 3-2 Compute the sigmaSq and etaSq then save those values
  Perm.NLL <- optimx(par=theta0, fn=NLL, gr=Gradient.NLL, hess=Hessian.NLL, method="BFGS");
  Output.Perm[q,1] <- seed.num;
  Output.Perm[q,2:4] <- Perm.NLL[1, c(1,2,7)];
};
## 3-3 Write in files;
FileName <- paste("output_perm/", seed.num, ".dat", sep="");
write.table(Output.Perm, file=FileName, row.names=F);
## 3-4 Assess the p value of the test
emp <- mean(c(Output.Perm$Theta1,Output.True$Theta1[d]) <= Output.True$Theta1[d]);
Output.PVal[d,1] <- seed.num;
Output.PVal[d,2] <- emp;
cat("Output p.values of round", d, "in", num.dataset, "\n", sep=" ");
print(Output.PVal[d,2]);
print("------------------------");
cat("\n", Sys.time()-tempus, "\n");
};

write.table(x=Output.True, file=paste("output_true", seed.num, ".dat", sep=""), row.names=TRUE);
write.table(x=Output.PVal, file=paste("output_pval", seed.num, ".dat", sep=""), row.names=FALSE);



