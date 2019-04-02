# Tree modifications 
# used for modifying different attributes of a tree and returning the result
 

require(ape)
tree.folder <- Sys.glob("~/PycharmProjects/hiv-withinhost/7SampleTrees/prelim/*.tree.sample")
for (treefile in tree.folder){
  path <- strsplit(treefile,"prelim")[[1]][1]
  intree <- read.tree(treefile)
  filename <- basename(treefile)
  
  logname <- paste0(strsplit(filename, "_")[[1]][1], ".log")
  
  logfile <- read.csv(paste0("~/PycharmProjects/hiv-withinhost/6BEASTout/",logname), sep="\t", skip=3)
  
  loglen <- nrow(logfile) -1
  interval <- c(loglen*0.1+1,loglen+1)
  rescale.factor <- median(logfile$ucld.mean[interval[1]:interval[2]])
  intree$edge.length <- (intree$edge.length * rescale.factor)
  print(intree)
  write.tree(intree,paste0(path,"rescaled/",filename))
}
#modify all tree edge lengths using the value from the log file 


#output 
write.tree(intree,"30631-a.tree")


# ----------------------
#bifurcating only
intree <- multi2di(intree)

#adding in numeric node labels  (1 to Nnode)
nNode <- intree$Nnode
intree$node.label <- as.character(seq(from=1,to=nNode))

#find all cherries 
ntips <- Ntip(intree)
numtips <- tabulate(intree$edge[which(intree$edge[,2] <= ntips),1])
is.cherry <- which(numtips == 2)

print(is.cherry)
c.nodes <- data.frame()

print(intree$edge)
print(intree$edge.length)
for (node in is.cherry){
	#print(node)
  idx.cherry <- which(intree$edge[,1]==node)
	tips <- intree$edge[idx.cherry,2]
	lengths <- intree$edge.length[idx.cherry]
	#print(tips)
	node.name <- paste0("(",intree$tip.label[tips[1]] , ":",lengths[1],",", intree$tip.label[tips[2]],":",lengths[2], ")")
	print(node.name)  
	#c.nodes <- rbind(c.nodes, data.frame(name=node, ,]
  
}
  



#write.tree(intree2,'~/historian/data/F1_relabel.tree')
