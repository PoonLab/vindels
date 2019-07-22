
logfiles <- Sys.glob("~/PycharmProjects/hiv-withinhost/6BEASTout/*.log")
logfiles <- logfiles[-5]
smeans <- c()
rmeans <- c()
files <- c()
for (fullpath in logfiles){
  filename <- basename(fullpath)
  strpath <- paste0("~/PycharmProjects/hiv-withinhost/6_5_BEASTstrict/", filename)
  if (file.exists(strpath)){
    strict <- read.table(strpath, sep="\t", header=T)
    shm <- 1/mean(1/strict$likelihood)
    
    relaxed <- read.table(fullpath, sep="\t", header=T)
    rhm <- 1/mean(1/relaxed$likelihood)
    
    smeans <- c(smeans, shm)
    rmeans <- c(rmeans, rhm)
    files <- c(files, filename)
  }

}

bayes <- data.frame(filename=files, relaxed=rmeans, strict=smeans)

