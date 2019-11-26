# Tree modifications 
# used for modifying different attributes of a tree and returning the result

require(ape)

args <- commandArgs(trailingOnly = T)

# input directory of sampled BEAST trees
# relies on the presence of a "prelim" folder being present
if (length(args) != 1){
  print("USAGE: Rscript 7_5_mcc_tree_mod.r [working directory]")
  quit()
}
for (i in 1:length(args)){
  if (!endsWith(args[i],"/")){
    args[i] <- paste0(args[i],"/") 
  }
}

treefolder <- args[1]

if (!dir.exists(paste0(treefolder,"prelim/"))){
  quit(status="USAGE: Rscript 7_5_mcc_tree_mod.r [working directory]")
}

csv <- read.csv("~/PycharmProjects/hiv-withinhost/6_hm/bayes-comparison2.csv", row.names=1,stringsAsFactors = F)

logfiles <- csv[,"filename"]
treefiles <- unname(sapply(logfiles,function(x){paste0(strsplit(x,"\\.")[[1]][1],".tree")}))
#treefiles <- unname(sapply(treefiles,function(x){sub("-original","",x)}))

best <- csv[,"best"]
choices <-  c("24BEAST-constant-final/","30BEAST-skygrid-narrow/","33BEAST-skygrid-10/","37BEAST-skygrid-30/")[best]
basepath <- "~/PycharmProjects/hiv-withinhost/6_hm/"

# for writing BEAST log files to a final folder 
for (i in 1:length(logfiles) ){
  file.copy(paste0(basepath,choices[i],logfiles[i]), paste0(basepath, "final/"), overwrite=T)
}


dir.create(paste0(treefolder,"final/"), showWarnings = FALSE)
infolder <- Sys.glob(paste0(treefolder,"prelim/*.tree"))

# for (f in folders){
#   trees <- Sys.glob(paste0(f, "/*.tree.sample"))
#   foldername <- paste0(basename(f),"/") 
#   
#   dir.create(paste0(treefolder,"final_multi/",foldername), showWarnings = FALSE)
#   
r.vec <- c()
for (i in 1:length(treefiles)){
  intree <- paste0("~/PycharmProjects/hiv-withinhost/7_5_MCC/prelim/",choices[i],treefiles[i])
  tre <- read.nexus(intree)
  
  print("INPUT")
  print(intree)
  
  inlog <- paste0("~/PycharmProjects/hiv-withinhost/6_hm/",choices[i],logfiles[i])
  # uses log file name to find and read BEAST log file
  logfile <- read.csv(inlog, sep="\t", skip=4)
  
  print("LOGFILE")
  print(inlog)
  
  #counts the number of MCMC steps
  loglen <- nrow(logfile) -1
  print(loglen)
  
  if (loglen != 10000){
    print(status="INCOMPLETE FILE")
    next()
  }
  
  # calculates the start and end interval of MCMC steps AFTER the burn in (assuming last 90%)
  interval <- c(loglen*0.1+1,loglen+1)
  print(interval)

  # calculates the rescale factor using the median of the UCLD.MEAN column (can check that this matches UCLD.MEDIAN on tracer)
  means <- logfile$ucld.mean[interval[1]:interval[2]]

  rescale.factor <- median(means)
  print(rescale.factor)
  
  r.vec[i] <- rescale.factor

  # rescales all the edge lengths of the tree
  tre$edge.length <- (tre$edge.length * rescale.factor)

  # writes the rescaled tree to a new folder called "final"
  write.tree(tre,paste0(treefolder,"final/",treefiles[i]))
}
#}

