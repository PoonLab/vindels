
infolder <- Sys.glob("~/PycharmProjects/hiv-withinhost/6_6_e5clock/*.log")

all.median <- c()
all.stdev <- c()

for (f in infolder){
  
  logname <- basename(f)
  
  # uses log file name to find and read BEAST log file
  logfile <- read.csv(f, sep="\t", skip=3)
  
  print(logname)
  
  #counts the number of MCMC steps
  loglen <- nrow(logfile) -1
  print(loglen)
  
  # calculates the start and end interval of MCMC steps AFTER the burn in (assuming last 90%)
  interval <- c(loglen*0.1+1,loglen+1)
  print(interval)
  
  print(summary(log(logfile$ucld.mean[interval[1]:interval[2]])))
  # calculates the rescale factor using the median of the UCLD.MEAN column (can check that this matches UCLD.MEDIAN on tracer)
  um.median <- median(logfile$ucld.mean[interval[1]:interval[2]])
  um.stdev <- sd(logfile$ucld.mean[interval[1]:interval[2]])
  
  all.median <- c(all.median, um.median)
  all.stdev <- c(all.stdev, um.stdev)
  
}
