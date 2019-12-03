require(ape)
args <- commandArgs(trailingOnly = T)
args[1] <- "~/PycharmProjects/hiv-withinhost/4_5_Raxml/100BS/"
#args[1] <- "~/PycharmProjects/hiv-withinhost/4_5_Raxml/100BS-all/"
path <- args[1]
if (!endsWith(path,"/")){
  path <- paste0(path,"/") 
}

trefolder <- Sys.glob(paste0(path,"RAxML_bestTree*"))
dir.create(paste0(path,"rooted_trees/"), showWarnings = FALSE)
dir.create(paste0(path,"guide_trees/"), showWarnings = FALSE)

print("ROOTING TREES ...")
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
rates <- c()
count <- 0
subtype <- c()
n <- 0
vn <- 0
treeroot <- c()
daterange <- c()
print("CREATING GUIDE TREES ...")
for (file in rtdfolder){
  n <- n + 1
  filename <- basename(file)
 
  name <- strsplit(filename, "\\.")[[1]][1]
  names[n] <- name
  rtd <- read.tree(file)
  lens <- node.depth.edgelength(rtd)[1:Ntip(rtd)]
  tip.dates <- as.numeric(unname(sapply(rtd$tip.label, function(x) strsplit(x, "_")[[1]][2])))

  print(filename)
  daterange[n] <- diff(range(tip.dates))
  
  # create a linear model and save it
  linear <- lm(lens ~ tip.dates)
  rates[n] <- coef(linear)[[2]]
  # extract rsqrd value
  x <- summary(linear)
  rsqr <- c(rsqr, x$r.squared)
  
  # counts
  count <- count + length(rtd$tip.label)
  subtype <- c(subtype, sapply(rtd$tip.label, function(x)strsplit(x,"\\.")[[1]][1]))
  
  # create a figure and save it
  png(file=paste("~/vindels/Figures/root-to-tip/all-hm/",name,"-rtt.png",sep=""),width=800,height=600, res=120)
  plot(lens ~ tip.dates, main=name, xlab="Collection Date (Days since a start point)", ylab="Root to tip branch length (Expected subs/site)")
  abline(linear)
  dev.off()
  
  # calculate the tree root height for use in BEAST Skygrid
  xint <- -coef(linear)[[1]]/coef(linear)[[2]]
  treeroot[n] <- ceiling((max(tip.dates) - xint)* 1.25)
  
  write.tree(rtd, file=paste0(path,"guide_trees/", filename))
}
#write.csv(data.frame(file=names, root_height=treeroot), "~/PycharmProjects/hiv-withinhost/4_5_Raxml/100BS/root-heights.csv",quote=F,row.names = F)

#write.csv(data.frame(file=names,date_range=daterange),"~/PycharmProjects/hiv-withinhost/date-ranges2.csv",quote=F,row.names = F)