# This script will read in NEWICK phylogenetic trees from RAXML 
# and use the RTT function in APE to root them. 
# Rooted trees will be published to the "rooted_trees" directory
# These rooted trees will then be read again in order to generate 
# guide trees that have all their branch lengths scaled to 1.

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

# print("ROOTING TREES ...")
# for (file in trefolder){
#   tre <- read.tree(file)
#   filename <- strsplit(basename(file), "bestTree.") [[1]][2]
#   print(filename)
# 
#   tip.dates <- sapply(tre$tip.label, function(x){strsplit(x, "_")[[1]][2]})
# 
#   rtd <- rtt(tre, as.numeric(tip.dates))
# 
#   #write.tree(rtd, file=paste0(path, "rooted_trees/", filename,2))
# }

rtdfolder <- Sys.glob(paste0(path,"rooted_trees/*.tree"))

rsqr <- c()
names <- c()
upper <- c()
lower <- c()
rates <- c()
count <- 0
subtype <- c()
n <- 0
vn <- 0
treeroot <- c()
daterange <- c()
print("CREATING GUIDE TREES ...")
for (i in 1:length(rtdfolder)){
  n <- n + 1
  filename <- basename(rtdfolder[i])
 
  name <- strsplit(filename, "\\.")[[1]][1]
  names[n] <- name
  rtd <- read.tree(rtdfolder[i])
  lens <- node.depth.edgelength(rtd)[1:Ntip(rtd)]
  tip.dates <- as.numeric(unname(sapply(rtd$tip.label, function(x) strsplit(x, "_")[[1]][2])))

  #print(filename)
  daterange[n] <- diff(range(tip.dates))
  
  # create a linear model and save it
  linear <- lm(lens ~ tip.dates)
  rates[n] <- coef(linear)[[2]]
  # extract rsqrd value
  x <- summary(linear)
  rsqr <- c(rsqr, x$r.squared)
  lower[i] <- confint(linear)[2]
  upper[i] <- confint(linear)[4]
  
  # counts
  count <- count + length(rtd$tip.label)
  subtype <- c(subtype, sapply(rtd$tip.label, function(x)strsplit(x,"\\.")[[1]][1]))
  
  # create a figure and save it
  #png(file=paste("~/vindels/Figures/root-to-tip/all-hm/",name,"-rtt.png",sep=""),width=800,height=600, res=120)
  #plot(lens ~ tip.dates, main=name, xlab="Collection Date (Days since a start point)", ylab="Root to tip branch length (Expected subs/site)")
  #abline(linear)
  #dev.off()
  
  # calculate the tree root height for use in BEAST Skygrid
  # when the PROJECTED start of infection took place (plus a 25% buffer)
  xint <- -coef(linear)[[1]]/coef(linear)[[2]]
  treeroot[n] <- ceiling((max(tip.dates) - xint)* 1.25)
  
  #write.tree(rtd, file=paste0(path,"guide_trees/", filename,2))
}
sapply(1:length(names),function(x){
  print(names[x])
  print(rsqr[x])
  print(lower[x])
  print(upper[x])
  
})

#write.csv(data.frame(file=names, root_height=treeroot), "~/PycharmProjects/hiv-withinhost/4_5_Raxml/100BS/root-heights.csv",quote=F,row.names = F)

