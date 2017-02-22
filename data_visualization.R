library(ggplot2);
Seed.Base <- 1000;  #____
Num.Dataset <- 100;  #____
Ggtitle <- "Simulation 1"; #____
SigmaSq.True <- 2;  #____
EtaS.True <- 5;  #____
## Load output data
Result.True <- read.table("output_true1100.dat", header=TRUE);
Result.Pval <- read.table("output_pval1100.dat", header=TRUE);
## Plot estimates of sigmaSq  
hp1 <- ggplot(data=Result.True);
hp1 + geom_histogram(aes(x=Theta1, fill=factor(convcode)), bins=30, alpha=0.3) +
      geom_vline(xintercept=SigmaSq.True, lty=1, colour="grey") + xlab("sigmaSq") + ggtitle(Ggtitle);
ggsave(paste("plot/histogram_sigmaSq_", Seed.Base+Num.Dataset, ".png", sep=""), plot=last_plot(), width=5, height=3);
## Plot estimate of etaSq  
hp1 + geom_histogram(aes(x=Theta2, fill=factor(convcode)), bins=30, alpha=0.3) +
      geom_vline(xintercept=EtaS.True, lty=1, colour="grey") + xlab("etaSq") + ggtitle(Ggtitle);
ggsave(paste("plot/histogram_etaSq_", Seed.Base+Num.Dataset, ".png", sep=""), plot=last_plot(), width=5, height=3);
## Plot empirical p-values
Result.Pval$convcode <- Result.True$convcode;
Reject.Rate <- mean(Result.Pval$Empirical.p <= 0.05);
hp2 <- ggplot(data=Result.Pval);
hp2 + geom_histogram(aes(x=Empirical.p, fill=factor(convcode)), bins=30, alpha=0.3) +
      geom_vline(xintercept=0.05, lty=1, colour="grey") + xlab("p") + ggtitle(paste(Ggtitle, " RR = ", Reject.Rate, sep=""));
ggsave(paste("plot/histogram_p-val_", Seed.Base+Num.Dataset, ".png", sep=""), plot=last_plot(), width=5, height=3);



