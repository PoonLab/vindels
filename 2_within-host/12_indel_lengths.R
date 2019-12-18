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
    }else if(x > 9){
      ">9"
    }
  })))
  category
}
getLength <- function(seq){
  
}
path <- "~/PycharmProjects/hiv-withinhost/"
iLength <- read.csv(paste0(path,"12_lengths/ins-full.csv"), row.names=1, stringsAsFactors = F)
dLength <- read.csv(paste0(path,"12_lengths/del-full.csv"), row.names=1, stringsAsFactors = F)
 
iLength <- iLength[iLength$Count>0,]
dLength <- dLength[dLength$Count>0,]

iLength$Len <- sapply(iLength$Seq, nchar)
dLength$Len <- sapply(dLength$Seq, nchar)

iLength$Bin <- sapply(iLength$Seq,categorize)
dLength$Bin <- sapply(dLength$Seq,categorize)

# verify that  no commas are found within the iLength and dLength data frames 
"," %in% iLength$Seq
"," %in% dLength$Seq

# order the iLength and dLength dataframes 
iLength$Bin <- factor(iLength$Bin,levels=c(">9","9","7-8", "6", "4-5","3", "1-2"))
dLength$Bin <- factor(dLength$Bin,levels=c(">9","9","7-8", "6", "4-5","3", "1-2"))

# table manipulation for data display
itab <- table(iLength$Bin, iLength$Vloop)
dtab <- table(dLength$Bin, dLength$Vloop)

idf <- as.data.frame(itab)
ddf <- as.data.frame(dtab)

colnames(idf) <- c("Bin", "Vloop", "Count")
colnames(ddf) <- colnames(idf)

# STACK BAR PLOT 
# -------------------------------------------
require(RColorBrewer)
pal <- c("gray28", "blue4",  'tomato', 'dodgerblue',  'red',  "skyblue", 'darkred' )
require(ggplot2)
iplot <- ggplot() + 
  geom_bar(aes(x=Vloop, y=Count, fill=Bin), data=idf, stat='identity') + 
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
  geom_bar(aes(x=Vloop, y=Count, fill=Bin), data=ddf, stat='identity') + 
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
iplot <- ggplot(idf, aes(y=Count, axis1=)) + 
  geom_bar(aes(x=Vloop, y=Count, fill=Bin), data=idf, stat='identity') + 







mosaic(~Vloop + Bin, data=iLength,
       shade=T, main=NULL, direction="v",
       spacing=spacing_equal(sp = unit(0.7, "lines")),
       residuals_type="Pearson",
       margins=c(2,2,6,2),
       labeling_args = list(tl_labels = c(F,T), 
                            tl_varnames=c(F,T),
                            gp_labels=gpar(fontsize=20),
                            gp_varnames=gpar(fontsize=26),
                            set_varnames = c(Vloop="Variable Loop", 
                                             Bin="Insertion Length (nt)"),
                            offset_labels=c(0,0,0,0),rot_labels=c(0,0,0,0), just_labels=c("center","center","center","center")),
       legend=legend_resbased(fontsize = 20, fontfamily = "",
                              x = unit(0.5, "lines"), y = unit(2,"lines"),
                              height = unit(0.8, "npc"),
                              width = unit(1, "lines"), range=c(-10,10)),
       set_labels=list(Vloop=c("V1","V2","V4","V5")))


mosaic(~Vloop + Bin, data=dLength,
       shade=T, main=NULL, direction="v",
       spacing=spacing_equal(sp = unit(0.7, "lines")),
       residuals_type="Pearson",
       margins=c(2,2,6,2),
       labeling_args = list(tl_labels = c(F,T), 
                            tl_varnames=c(F,T),
                            gp_labels=gpar(fontsize=20),
                            gp_varnames=gpar(fontsize=26),
                            set_varnames = c(Vloop="Variable Loop", 
                                             Bin="Deletion Length (nt)"),
                            offset_labels=c(0,0,0,0),rot_labels=c(0,0,0,0), just_labels=c("center","center","center","center")),
       legend=legend_resbased(fontsize = 20, fontfamily = "",
                              x = unit(0.5, "lines"), y = unit(2,"lines"),
                              height = unit(0.8, "npc"),
                              width = unit(1, "lines"), range=c(-10,10)),
       set_labels=list(Vloop=c("V1","V2","V4","V5")))



ggplot()
ggplot(all.df, aes(x=Date, y=Count, group=Count)) + geom_density_ridges(colour="white", fill="blue", scale=1, bandwidth=5)
