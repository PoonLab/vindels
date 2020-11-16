args <- commandArgs(trailingOnly=T)

if (length(args) < 2){

	print("USAGE: Rscript bayes-factor.R [input dir 1] [input dir 2] ...")
	quit()
}
for (i in 1:length(args)){
	if (!endsWith(args[i],"/")){
		args[i] <- paste0(args[i],"/")
	}
}
commonfiles <- c()
for (i in 1:length(args)){
  logfiles <- c()
  logfiles <- Sys.glob(paste0(args[i],"*.log"))
  logfiles <- unname(sapply(logfiles, basename))
  #print(logfiles)
  if (i == 1){
    commonfiles <- logfiles
  }else{
    toRemove <- sapply(commonfiles, function(x){if (!(x %in% logfiles)) x })
    toRemove <- -c(which(toRemove %in% commonfiles))
    if (length(toRemove) > 0){
      commonfiles <- commonfiles[toRemove]
    }
  }
}
print(commonfiles)
bayes <- data.frame(filename=commonfiles, stringsAsFactors = F)
runs <- c()
for (j in 1:length(args)){
  logfiles <- paste0(args[j],commonfiles)
  means <- c()
  
  print(paste0("DIRECTORY ",j," STARTING..."))

  # ITERATOR: DIRECTORY NUMBER 1
  for (k in 1:length(logfiles)){
    filename <- basename(logfiles[k])
    print(filename)
  
    # DIRECTORY NUMBER 1 PROCESSING
    data <- read.table(logfiles[k], sep="\t", skip=4, head=T)
    loglen <- nrow(data) -1
    #print(loglen)
    
    # calculates the start and end interval of MCMC steps AFTER the burn in (assuming last 90%)
    interval <- c(loglen*0.1+1,loglen+1)
    #print(interval)

    llh <- mean(data$likelihood[interval[1]:interval[2]])
    #print(llh)
    means[k] <- llh
  }
  bayes[,j+1] <- means
  runs[j] <- args[j]
}

bayes$best <- apply(bayes[,1:j+1], 1, function(x) which(x == max(x)))
colnames(bayes) <- c("filename", runs, "best")
write.csv(bayes, "~/PycharmProjects/hiv-withinhost/6_hm/bf-comparison2.csv")
for (x in 1:j){
  print(paste0("DIR",x," Count: ", sum(bayes$best == x)))
}

#print(paste0("DIR2 Count: ", sum(bayes$better == 2)))
