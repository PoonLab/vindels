require(ape)

trefolder <- Sys.glob("~/PycharmProjects/hiv-withinhost/4_5_Raxml/signal/RAxML_bestTree*")

for (file in trefolder){
  tre <- read.tree(file)
  filename <- strsplit(basename(file), "bestTree.") [[1]][2]
  print(filename)

  tip.dates <- sapply(tre$tip.label, function(x){strsplit(x, "_")[[1]][2]})

  rtd <- rtt(tre, as.numeric(tip.dates))

  write.tree(rtd, file=paste0("~/PycharmProjects/hiv-withinhost/4_5_Raxml/signal/rooted_trees/", filename))
}

rtdfolder <- Sys.glob("~/PycharmProjects/hiv-withinhost/4_5_Raxml/signal/rooted_trees/*.tree")

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
  ranges <- c(ranges, max(tip.dates))
  # create a linear model and save it 
  linear <- lm(lens ~ tip.dates)
  x <- summary(linear)
  rsqr <- c(rsqr, x$r.squared)
  count <- count + length(rtd$tip.label)
  if (n > 25){
    print(filename) 
    vn <- vn + length(rtd$tip.label)
  }
  
  subtype <- c(subtype, sapply(rtd$tip.label, function(x)strsplit(x,"\\.")[[1]][1]))
  # create a figure and save it 
  png(file=paste("~/vindels/Figures/root-to-tip/",name,"-rtt.png",sep=""),width=800,height=600, res=120)
  plot(jitter(lens) ~ jitter(tip.dates), main=name, xlab="Collection Date (Days since a start point)", ylab="Root to tip branch length (Expected subs/site)")
  abline(linear)
  dev.off()
  
  
  #write.tree(rtd, file=paste0("~/PycharmProjects/hiv-withinhost/4_5_Raxml/guide_trees/", filename))
}
print(table(subtype))
names(ranges) <- names
names(rsqr) <- names

print(vn)
print(count)
