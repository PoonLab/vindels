# Tree modifications 
# used for modifying different attributes of a tree and returning the result

require(ape)

#args <- commandArgs(trailingOnly = T)

#for (i in 1:length(args)){
#  if (!endsWith(args[i],"/")){
#    args[i] <- args[i] + "/"
#  }
#}

# input directory of sampled BEAST trees
# relies on the presence of a "prelim" folder being present
args <- c()
args[1] <- "~/PycharmProjects/hiv-withinhost/7_5_MCC/prelim/"
args[2] <- "~/PycharmProjects/hiv-withinhost/7_5_MCC/rescaled/"


infolder <- Sys.glob(paste0(args[1],"*.tree"))
outpath <- args[2]

# for (f in folders){
#   trees <- Sys.glob(paste0(f, "/*.tree.sample"))
#   foldername <- paste0(basename(f),"/") 
#   
#   dir.create(paste0(path,"rescaled_multi/",foldername), showWarnings = FALSE)
#   
for (treefile in infolder){
  
  print(treefile)
  
  intree <- read.nexus(treefile)
  filename <- basename(treefile)

  # edits the tree file name to get the BEAST log file name
  logname <- paste0(strsplit(filename, "\\.")[[1]][1], ".log")

  # uses log file name to find and read BEAST log file
  logfile <- read.csv(paste0("~/PycharmProjects/hiv-withinhost/6BEASTout-comb/",logname), sep="\t", skip=3)

  print(logname)

  #counts the number of MCMC steps
  loglen <- nrow(logfile) -1
  print(loglen)

  # calculates the start and end interval of MCMC steps AFTER the burn in (assuming last 90%)
  interval <- c(loglen*0.1+1,loglen+1)
  print(interval)

  # calculates the rescale factor using the median of the UCLD.MEAN column (can check that this matches UCLD.MEDIAN on tracer)
  rescale.factor <- median(logfile$ucld.mean[interval[1]:interval[2]])
  print(rescale.factor)

  # rescales all the edge lengths of the tree
  intree$edge.length <- (intree$edge.length * rescale.factor)

  # writes the rescaled tree to a new folder called "rescaled"
  write.tree(intree,paste0(outpath, filename))
}
#}

