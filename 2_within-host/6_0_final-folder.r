args <- commandArgs(trailingOnly=T)

if (length(args) <= 1){
    print("USAGE: Rscript 6_0_final_folder [csv file] [dir 1] [dir2]")
    quit()
}

csv <- read.csv(args[1], row.names=1,stringsAsFactors = F)

logfiles <- csv[,"filename"]
names <- unname(sapply(logfiles,function(x){strsplit(x,"\\.")[[1]][1]}))
treefiles <- paste0(names,".time.trees")
opsfiles <- paste0(names,".ops")
#treefiles <- unname(sapply(treefiles,function(x){sub("-original","",x)}))

dir.create("~/work/6_output/final/", showWarnings=F)

final <- csv[,"final"]
choices <- args[2:length(args)][final]

basepath <- "~/work/6_output/"

dir.create(paste0(basepath,"final/trees/"), showWarnings = F)
dir.create(paste0(basepath, "final/operators/"), showWarnings = F)

# for writing BEAST log files to a final folder 
for (i in 1:length(logfiles) ){
  file.copy(paste0(basepath,choices[i],logfiles[i]), paste0(basepath, "final/"), overwrite=T)
  file.copy(paste0(basepath,choices[i],"trees/",treefiles[i]), paste0(basepath, "final/trees/"), overwrite=T)
  file.copy(paste0(basepath,choices[i],"operators/", opsfiles[i]), paste0(basepath, "final/operators/"), overwrite=T)
}
