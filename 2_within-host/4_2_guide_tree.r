# RESCALE THE GUIDE TREE FOR BEAST ANALYSIS
# USAGE
# Arg 1 : Input tree file to be reformatted
# Arg 2 : Output tree file to write 

require(ape)  # if ape is not installed, run `install.packages("ape")`
args <- commandArgs(trailingOnly = T)

treefolder <- Sys.glob("~/PycharmProjects/hiv-withinhost/4_5_Raxml/RAxML_bestTree.*")

for (infile in treefolder){
    filename <- basename(infile)
    tr <- read.tree(infile)
    tr$edge.length <- rep(1, times=length(tr$edge.length))
    write.tree(tr, paste0("~/PycharmProjects/hiv-withinhost/4_5_Raxml/guide_trees/",filename))
}

