
csv <- read.csv("~/PycharmProjects/hiv-withinhost/6_hm/bayes-comparison5.csv", row.names=1,stringsAsFactors = F)

logfiles <- csv[,"filename"]
names <- unname(sapply(logfiles,function(x){strsplit(x,"\\.")[[1]][1]}))
treefiles <- paste0(names,".time.trees")
opsfiles <- paste0(names,".ops")
#treefiles <- unname(sapply(treefiles,function(x){sub("-original","",x)}))

final <- csv[,"final"]
choices <-  c("24BEAST-constant-final/","30BEAST-skygrid-narrow/",
              "33BEAST-skygrid-10/","37BEAST-skygrid-30/", 
              "41BEAST-rwalk/","43BEAST-skygrid-10/","44BEAST-skygrid-30/")[final]
basepath <- "~/PycharmProjects/hiv-withinhost/6_hm/"

dir.create(paste0(basepath,"final/trees/"), showWarnings = F)
dir.create(paste0(basepath, "final/operators/"), showWarnings = F)

# for writing BEAST log files to a final folder 
for (i in 1:length(logfiles) ){
  file.copy(paste0(basepath,choices[i],logfiles[i]), paste0(basepath, "final/"), overwrite=T)
  file.copy(paste0(basepath,choices[i],"trees/",treefiles[i]), paste0(basepath, "final/trees/"), overwrite=T)
  file.copy(paste0(basepath,choices[i],"operators/", opsfiles[i]), paste0(basepath, "final/operators/"), overwrite=T)
}
