require(ape)
args <- commandArgs(trailingOnly=T)

# if (length(args) != 1){
#   print("USAGE: Rscript treetime-parsing.r [input directory]")
#   quit()
# }
# for (i in 1:length(args)){
#   if (!endsWith(args[i],"/")){
#     args[i] <- paste0(args[i],"/") 
#   }
# }

treefolder <- Sys.glob(paste0("~/PycharmProjects/hiv-withinhost/4_5_Raxml/100BS/RAxML_bestTree*"))
dir.create("~/PycharmProjects/hiv-withinhost/4_5_Raxml/100BS/treetime_test2/dates/")
for (infile in treefolder){
  filename <- strsplit(basename(infile),"\\.")[[1]][2]
  print(filename)
  tre <- read.tree(infile)
  
  # dates <- unname(sapply(tre$tip.label, function(x){
  #   fields <- strsplit(x,"_")[[1]]
  #   
  #   days <- as.double(fields[2])/365.25
  #   year <- strsplit(fields[1],"\\.")[[1]][3]
  #   if (year == "-"){
  #     year <- 0
  #   }else{
  #     year <- as.numeric(year)
  #   }
  #   year + days
  # }))
  dates <- unname(sapply(tre$tip.label, function(x)strsplit(x,"_")[[1]][2]))
  
  output <- data.frame("name"=tre$tip.label, "date"=dates)
  #print(head(output))
  write.table(output, paste0("~/PycharmProjects/hiv-withinhost/4_5_Raxml/100BS/treetime_test2/dates/",filename,".tsv"),quote=F,sep="\t", row.names=F)
  
}
