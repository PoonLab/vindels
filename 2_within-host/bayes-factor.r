logfiles <- Sys.glob("~/PycharmProjects/hiv-withinhost/6_4_2_BEASTout-lognorm/*.log")
smeans <- c()
rmeans <- c()
files <- c()
for (fullpath in logfiles){
  filename <- basename(fullpath)
  strpath <- paste0("~/11BEAST-constant/output/", filename)
  #print(filename)
  if (file.exists(strpath)){
    print(filename)
    print(strpath)
    strict <- read.table(strpath, sep="\t", skip=3, head=T)
    loglen <- nrow(strict) -1
    print(loglen)
    if (loglen != 10000){
      next
    }
    # calculates the start and end interval of MCMC steps AFTER the burn in (assuming last 90%)
    interval <- c(loglen*0.1+1,loglen+1)
    print(interval)
    
    shm <- 1/mean(1/strict$likelihood[interval[1]:interval[2]])
    #print(strict$likelihood) 
    
    relaxed <- read.table(fullpath, sep="\t", skip=3, head=T)
    loglen <- nrow(relaxed) -1
    print(loglen)
    
    # calculates the start and end interval of MCMC steps AFTER the burn in (assuming last 90%)
    interval <- c(loglen*0.1+1,loglen+1)
    print(interval)
    
    rhm <- 1/mean(1/relaxed$likelihood[interval[1]:interval[2]])
   
    print(rhm)
    print(shm)
     
    smeans <- c(smeans, shm)
    rmeans <- c(rmeans, rhm)
    files <- c(files, filename)
  }

}

bayes <- data.frame(filename=files, skygrid=rmeans, constant=smeans)

