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
logfolder <- Sys.glob(paste0(logpath,"*.log"))
heights <- c()
for (t in 1:length(trefolder)){
  
  print(trefolder[t])
  filename <- basename(logfolder[l])
  pat <- strsplit(filename, "\\.")[[1]][1]

  pat_trees <- Sys.glob(paste0(tpath,pat,"*"))
  
  for (tree in pat_trees){
    # treename1 <- paste0(strsplit(basename(hfolder[h]), "_recon")[[1]][1], ".tree.sample")
    # 
    # treename1 <- paste0(tfolder, treename1)
    # 
    # intree <- read.tree(treename1)
    # 
    # heights[h] <- max(node.depth.edgelength(intree))
    # uses log file name to find and read BEAST log file
    logfile <- read.csv(logfolder[l], sep="\t", skip=4)

    treename <- basename(tree)
    state <- as.numeric(paste0(strsplit(strsplit(treename, "\\.")[[1]][1], "_")[[1]][3],"00000"))
  
    intree <- read.tree(tree)
    rescale.factor <- logfile[which(logfile$state == state),"ucld.mean"]

    #print(rescale.factor)
  
    # rescales all the edge lengths of the tree
    intree$edge.length <- (intree$edge.length * rescale.factor)
    
    #writes the rescaled tree to a new folder called "rescaled"
    write.tree(intree,paste0(outfolder, gsub(".sample","",treename)))
  }
}

