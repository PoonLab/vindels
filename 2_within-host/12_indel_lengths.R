# LENGTH ANALYSIS
# --------------------------------------------------------
# script to process indel lengths 


categorize <- function(seqList){
  seqs <- strsplit(seqList, ",")[[1]]
  if (identical(seqs, character(0))){
    return(c())
  }
  lengths <- nchar(seqs)
  category <- c(unname(sapply(lengths, function(x){
    if (x <3){
      "1-2"
    }else if(x == 3){
      "3"
    }else if(x == 4 | x == 5){
      "4-5"
    }else if(x == 6){
      "6"
    }else if(x == 7 | x == 8){
      "7-8"
    }else if(x == 9){
      "9"
    }else if(x == 10 | x == 11){
      "10-11"
    }else if(x == 12){
      "12"
    }else{
      ">12"
    }
  })))
  category
}

path <- "~/PycharmProjects/hiv-withinhost/"
#path <- "~/Lio/"
ins <- read.csv(paste0(path,"12_lengths/all/ins-all.csv"), row.names=1, stringsAsFactors = F)
del <- read.csv(paste0(path,"12_lengths/all/del-all.csv"), row.names=1, stringsAsFactors = F)
 

ins$len <- nchar(ins$indel)
del$len <- nchar(del$indel)

ins$bin <- sapply(ins$indel,categorize)
del$bin <- sapply(del$indel,categorize)

#del <- del[nchar(del$indel) < 200, ]

# verify that  no commas are found within the ins and del data frames 
"," %in% ins$indel
"," %in% del$indel

# order the ins and del dataframes 
labels <- c(">12","12","10-11","9","7-8", "6", "4-5","3", "1-2")
labels <- labels[length(labels):1]  # TO REVERSE THE ORDER 

ins$bin <- factor(ins$bin,levels=labels)
del$bin <- factor(del$bin,levels=labels)

# table manipulation for data display
itab <- as.matrix(table(ins$bin, ins$vloop)) / 20
dtab <- as.matrix(table(del$bin, del$vloop) )/ 20

require(vcd)

idf <- as.data.frame(itab)
ddf <- as.data.frame(dtab)

colnames(idf) <- c("bin", "vloop", "count")
colnames(ddf) <- colnames(idf)

iround <- idf
dround <- ddf
iround$count <- round(iround$count)
dround$count <- round(dround$count)


# --- Mosaic Plot -----
data <- dround
df <- data.frame(bin=factor(rep(data$bin, data$count),levels=labels), vloop = rep(data$vloop, data$count))

# reorder the data frame 
df$bin <- factor(df$bin, levels=c("1-2","3","4-5","6","7-8","9","10-11","12",">12"))
df <- df[order(df$bin),]

require(vcd)
mosaic(~ bin + vloop,
       data = df,
       shade=T, main=NULL,
       spacing=spacing_equal(sp = unit(0.7, "lines")),
       residuals_type="Pearson", direction="v",
       margins=c(2,2,6,2),
       labeling_args = list(tl_labels = c(F,T), 
                            tl_varnames=c(F,T),
                            gp_labels=gpar(fontsize=24),
                            gp_varnames=gpar(fontsize=28),
                            set_varnames = c(vloop="Variable Loop", 
                                             bin="Indel Length (nt)"),
                            offset_labels=c(0,0,0,0),rot_labels=c(0,0,0,0), just_labels=c("center","center","center","center")),
       legend=legend_resbased(fontsize = 20, fontfamily = "",
                              x = unit(0.5, "lines"), y = unit(2,"lines"),
                              height = unit(0.8, "npc"),
                              width = unit(1, "lines"), range=c(-10,10)),
       set_labels=list(Variable.Loop=c("V1","V2","V3","V4","V5")))

# ---- Significance ---- 
# add in the significance level column
idf$sign <- rep(2,45)
# Higher
idf$sign[c(9,19,35,38,40)] <- 3  
# Lower
idf$sign[c(45)] <- 1    

ddf$sign <- rep(2,45)
# Higher
ddf$sign[c(9,10,15,19,35,38,40)] <- 3
# Lower
ddf$sign[c(1,6,17,31,44,45)] <- 1

# Proportion of frameshift indels 
x <- nchar(ins$indel)
sum(x[x%%3 != 0]) / sum(x)

x <- nchar(del$indel)
sum(x[x%%3 != 0]) / sum(x)



labels <- labels[length(labels):1]

# STACK BAR PLOT 
# -------------------------------------------
require(RColorBrewer)
pal <- c("gray28", "cyan","brown","blue4",  'tomato', 'dodgerblue',  'red',  "skyblue", 'darkred' )
pal <- pal[length(pal):1]

data <- ddf

ymx <- 510


#png(filename="~/vindels/Figures/within-host/finalized/del-length-v2", width=1200, height=700)
par(mar=c(6,7,2,8), xpd=F)
ax <- 1.7
lab <- 2.1
plot(NA, xlim=c(0,5), 
     ylim=c(0,ymx), 
     xaxt="n",
     xaxs="i",
     yaxs="i",
     xlab="",
     ylab="",
     cex.lab=lab, cex.axis=ax, las=1)
axis(1,at=seq(0.5,4.5), 
     labels=c("V1", "V2", "V3", "V4", "V5"),
     cex.axis=ax)
title(ylab="Average Deletion Counts", cex.lab=lab, line=4.2)
title(xlab="Variable Loop", cex.lab=lab, line=3.5)
abline(v=seq(0.5,4.5)-0.3, lty=1, col="gray68")
abline(v=seq(0.5,4.5)+0.3, lty=1, col="gray68")
#abline(h=seq(0,150,50), col="gray68")


for (i in seq(0.5,4.5)){
  pos <- 0
  d <- data[data$vloop==i+0.5,]
  #data <- data[nrow(data):1,]
  
  for (j in 1:9){
    s <- d[j,"sign"]
    n <- d[j,"count"]
    
    rect(i-0.15*s, pos, i+0.15*s, pos+n, col=pal[j])
    pos <- pos + n
  }
}
par(xpd=NA)
text(-0.48,ymx,"b)",cex=2)
pal <- pal[length(pal):1]
# legend(5.08,ymx*0.7,
#        legend=labels,cex=1.5, 
#        pch=22,pt.cex=5,pt.bg=pal,
#        y.intersp=1.3,
#        x.intersp=1.0,
#        title="Lengths")


require(ggplot2)
iplot <- ggplot() + 
  geom_bar(aes(x=vloop, y=count, fill=bin), data=idf, stat='identity') + 
  scale_fill_manual(values=rep(pal,4)) + 
  labs(x="Variable Loop", 
       y="Frequency", title="Insertion Lengths") +
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = 1.3, r = 1, b = 0.7, l = 1.5, unit = "cm"),
        axis.line = element_line(colour = "black"), 
        axis.title=element_text(size=18,margin=margin(t = 0, r = 3, b = 0, l = 12)),
        axis.text = element_text(size=16, colour="black"),
        plot.title = element_text(size=22, hjust = 0.5),
        legend.text=element_text(size=16), 
        legend.background=element_rect(colour="black"),
        legend.title=element_text(size=18))
iplot


dplot <- ggplot() + 
  geom_bar(aes(x=vloop, y=count, fill=bin), data=ddf, stat='identity') + 
  scale_fill_manual(values=rep(pal,4)) + 
  labs(x="Variable Loop", 
       y="Frequency", title="Deletion Lengths") +
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = 1.3, r = 1, b = 0.7, l = 1.5, unit = "cm"),
        axis.line = element_line(colour = "black"), 
        axis.title=element_text(size=18,margin=margin(t = 0, r = 3, b = 0, l = 12)),
        axis.text = element_text(size=16, colour="black"),
        plot.title = element_text(size=22, hjust = 0.5),
        legend.text=element_text(size=16), 
        legend.background=element_rect(colour="black"),
        legend.title=element_text(size=18))
dplot




# TEST : ALLUVIAL PLOT
# ----------------------------------------
iplot <- ggplot(idf, aes(y=count, axis1=)) + 
  geom_bar(aes(x=vloop, y=count, fill=bin), data=idf, stat='identity') + 







mosaic(~vloop + bin, data=ins,
       shade=T, main=NULL, direction="v",
       spacing=spacing_equal(sp = unit(0.7, "lines")),
       residuals_type="Pearson",
       margins=c(2,2,6,2),
       labeling_args = list(tl_labels = c(F,T), 
                            tl_varnames=c(F,T),
                            gp_labels=gpar(fontsize=20),
                            gp_varnames=gpar(fontsize=26),
                            set_varnames = c(vloop="Variable Loop", 
                                             bin="Insertion Length (nt)"),
                            offset_labels=c(0,0,0,0),rot_labels=c(0,0,0,0), just_labels=c("center","center","center","center")),
       legend=legend_resbased(fontsize = 20, fontfamily = "",
                              x = unit(0.5, "lines"), y = unit(2,"lines"),
                              height = unit(0.8, "npc"),
                              width = unit(1, "lines"), range=c(-10,10)),
       set_labels=list(vloop=c("V1","V2","V4","V5")))


mosaic(~vloop + bin, data=del,
       shade=T, main=NULL, direction="v",
       spacing=spacing_equal(sp = unit(0.7, "lines")),
       residuals_type="Pearson",
       margins=c(2,2,6,2),
       labeling_args = list(tl_labels = c(F,T), 
                            tl_varnames=c(F,T),
                            gp_labels=gpar(fontsize=20),
                            gp_varnames=gpar(fontsize=26),
                            set_varnames = c(vloop="Variable Loop", 
                                             bin="Deletion Length (nt)"),
                            offset_labels=c(0,0,0,0),rot_labels=c(0,0,0,0), just_labels=c("center","center","center","center")),
       legend=legend_resbased(fontsize = 20, fontfamily = "",
                              x = unit(0.5, "lines"), y = unit(2,"lines"),
                              height = unit(0.8, "npc"),
                              width = unit(1, "lines"), range=c(-10,10)),
       set_labels=list(vloop=c("V1","V2","V4","V5")))



ggplot()
ggplot(all.df, aes(x=Date, y=count, group=count)) + geom_density_ridges(colour="white", fill="blue", scale=1, bandwidth=5)
