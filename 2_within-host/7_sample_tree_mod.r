# Tree modifications 
# used for modifying different attributes of a tree and returning the result

require(ape)

args <- commandArgs(trailingOnly = T)

if (length(args) != 3){
  print("USAGE: Rscript 7_sample_tree_mod.r [log dir] [input dir] [output dir]")
  quit()
}

for (i in 1:length(args)){
  if (!endsWith(args[i],"/")){
    args[i] <- paste0(args[i],"/") 
  }
}

# input directory of sampled BEAST trees
# relies on the presence of a "prelim" folder being present
logpath <- args[1]
tpath <- args[2]
outfolder <- args[3]

dir.create(outfolder, showWarnings = F)
trefolder <- Sys.glob(paste0(tpath, "*"))

patnames <- unique(sapply(basename(trefolder), function(x) strsplit(x, "_")[[1]][1]))

for (i in 1:length(patnames)){
  
  filename <- paste0(patnames[i],".log")
  print(filename)
   
  # uses log file name to find and read BEAST log file
  logfile <- read.csv(paste0(logpath,filename), sep="\t", skip=4)
  pat_trees <- Sys.glob(paste0(tpath,patnames[i],"*"))
  
  for (t in 1:length(pat_trees)){
    treename <- basename(pat_trees[t])
    state <- as.numeric(paste0(strsplit(strsplit(treename, "\\.")[[1]][1], "_")[[1]][3],"00000"))
  
    intree <- read.tree(pat_trees[t])
    rescale.factor <- logfile[which(logfile$state == state),"ucld.mean"]

    #print(rescale.factor)
  
    # rescales all the edge lengths of the tree
    intree$edge.length <- (intree$edge.length * rescale.factor)
    
    #writes the rescaled tree to a new folder called "rescaled"
    write.tree(intree,paste0(outfolder, gsub(".sample","",treename)))
  }
}

