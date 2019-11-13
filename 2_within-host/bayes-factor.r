args <- commandArgs(trailingOnly=T)

if (length(args)!= 2){

	print("USAGE: Rscript bayes-factor.R [input dir 1] [input dir 2] ...")
	quit()
}
for (i in 1:length(args)){
	if (!endsWith(args[i],"/")){
		args[i] <- paste0(args[i],"/")
	}
}
commonfiles <- c()
for (dtry in args){
  logfiles <- c()
  logfiles <- Sys.glob(paste0(dtry,"*.log"))
  
}

for (dtry in args){
  logfiles <- Sys.glob(paste0(dtry,"*.log"))
  smeans <- c()
  rmeans <- c()
  files <- c()
  
  # ITERATOR: DIRECTORY NUMBER 1 
  for (fullpath in logfiles){
    filename <- basename(fullpath)
    
    # DIRECTORY NUMBER 2
    strpath <- paste0(args[2], filename)
    #print(filename)
    if (file.exists(strpath)){
      #print(filename)
      #print(strpath)
      
      # DIRECTORY NUMBER 1 PROCESSING 
      relaxed <- read.table(fullpath, sep="\t", skip=3, head=T)
      loglen <- nrow(relaxed) -1
      #print(loglen)
      
      # calculates the start and end interval of MCMC steps AFTER the burn in (assuming last 90%)
      interval <- c(loglen*0.1+1,loglen+1)
      #print(interval)
      
      rhm <- mean(relaxed$likelihood[interval[1]:interval[2]])
      
      # DIRECTORY NUMBER 2 PROCESSING
      strict <- read.table(strpath, sep="\t", skip=3, head=T)
      loglen <- nrow(strict) -1
      #print(loglen)
      if (loglen != 10000){
        next
      }
      # calculates the start and end interval of MCMC steps AFTER the burn in (assuming last 90%)
      interval <- c(loglen*0.1+1,loglen+1)
      #print(interval)
      
      shm <- mean(strict$likelihood[interval[1]:interval[2]])
      #print(strict$likelihood) 
      
      #print(rhm)
      #print(shm)
      
      smeans <- c(smeans, shm)
      rmeans <- c(rmeans, rhm)
      files <- c(files, filename)
    }
    
  }
  
}

bayes <- data.frame(filename=files, DIR1=rmeans, DIR2=smeans, diff=rmeans-smeans)
bayes$better <- mapply(function(x,y){if (x > y){1} else 2}, bayes$DIR1, bayes$DIR2)
print(bayes)
print(paste0("DIR1 Count: ", sum(bayes$better == 1)))
print(paste0("DIR2 Count: ", sum(bayes$better == 2)))
