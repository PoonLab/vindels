require(ape)
args <- commandArgs(trailingOnly = T)

path <- args[1]
if (!endsWith(path,"/")){
  path <- paste0(path,"/") 
}

trefolder <- Sys.glob(paste0(path,"RAxML_bestTree*"))
dir.create(paste0(path,"rooted_trees/"), showWarnings = FALSE)
dir.create(paste0(path,"guide_trees/"), showWarnings = FALSE)


for (file in trefolder){
  tre <- read.tree(file)
  filename <- strsplit(basename(file), "bestTree.") [[1]][2]
  print(filename)

  tip.dates <- sapply(tre$tip.label, function(x){strsplit(x, "_")[[1]][2]})

  rtd <- rtt(tre, as.numeric(tip.dates))

  write.tree(rtd, file=paste0(path, "rooted_trees/", filename))
}

rtdfolder <- Sys.glob(paste0(path,"rooted_trees/*.tree"))

rsqr <- c()
names <- c()
ranges <- c()
count <- 0
subtype <- c()
n <- 0
vn <- 0
for (file in rtdfolder){
  n <- n + 1
  filename <- basename(file)
 
  name <- strsplit(filename, "\\.")[[1]][1]
  names <- c(names,name)
  rtd <- read.tree(file)
  lens <- node.depth.edgelength(rtd)[1:Ntip(rtd)]
  tip.dates <- as.numeric(unname(sapply(rtd$tip.label, function(x) strsplit(x, "_")[[1]][2])))
  rng <- range(tip.dates)
  ranges <- c(ranges, rng)
  
  if (rng < 300){
    print(filename)
  }
  
  # # create a linear model and save it 
  # linear <- lm(lens ~ tip.dates)
  # x <- summary(linear)
  # rsqr <- c(rsqr, x$r.squared)
  # count <- count + length(rtd$tip.label)
  # 
  # subtype <- c(subtype, sapply(rtd$tip.label, function(x)strsplit(x,"\\.")[[1]][1]))
  # # create a figure and save it 
  # png(file=paste("~/vindels/Figures/root-to-tip/recomb-filter/",name,"-rtt.png",sep=""),width=800,height=600, res=120)
  # plot(jitter(lens) ~ jitter(tip.dates), main=name, xlab="Collection Date (Days since a start point)", ylab="Root to tip branch length (Expected subs/site)")
  # abline(linear)
  # dev.off()
  
  
  write.tree(rtd, file=paste0(path,"guide_trees/", filename))
}

