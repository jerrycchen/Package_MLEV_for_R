master <- read.table("1100.dat", header=TRUE);
for (seed in 1001:1100) {
  
  OutputName <- paste("output_perm/", seed, ".dat", sep="");
  
  Data <- subset(master, SEED==seed);
  
  write.table(x=Data, file=OutputName, row.names = FALSE);
  
};