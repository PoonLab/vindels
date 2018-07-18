require(ape)
require(ggtree)
setwd('~/git/vindels')

tr <- read.tree('7Trees/01_AE_MSA2_.fasta.tree')
dt <- read.csv('7Trees/01_AE_dates.txt', header=F)
dt$V2 <- as.character(dt$V2)

source('date-ranges.R')
dt$midpt <- apply(sapply(dt$V2, get.range), 2, mean)

df <- fortify(tr)
index <- match(df$label, dt$V1)
df$days <- dt$midpt[index]

tr2 <- rtt(tr, df$days[df$isTip])
y <- node.depth.edgelength(tr2)[1:Ntip(tr2)]
x <- df$days[df$isTip]
plot(x,y, xlab='Days since 1970', ylab='Distance from root', cex.lab=1.2)
