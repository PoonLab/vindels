infolder <- Sys.glob("~/23_newclock_psample/output/*.log")

all.median <- c()
all.stdev <- c()
i <- 0
names <- c()
for (f in infolder){
  i <- i + 1
  logname <- basename(f)
  
  # uses log file name to find and read BEAST log file
  logfile <- read.csv(f, sep="\t", skip=4)
  
  print(logname)
  
  names[i] <- logname
  #counts the number of MCMC steps
  loglen <- nrow(logfile) -1
  print(loglen)
  
  # calculates the start and end interval of MCMC steps AFTER the burn in (assuming last 90%)
  interval <- c(as.integer(loglen*0.1+1),as.integer(loglen+1))
  
  print(interval)
  
  values <- logfile$ucld.mean[interval[1]:interval[2]]
  print(summary(values))
  print(quantile(values, c(0.025,0.975)))
  # calculates the rescale factor using the median of the UCLD.MEAN column (can check that this matches UCLD.MEDIAN on tracer)
  um.median <- median(values)
  um.stdev <- sd(values)

  all.median <- c(all.median, um.median)
  all.stdev <- c(all.stdev, um.stdev)

}

