
library(ggplot2)
library(ggpubr)
library(VennDiagram)
library(wesanderson)
library(data.table)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggvenn)
library(ggsci)
library(ggrepel)
library(eulerr)

epigResults$change <- factor(epigResults$change, levels=c("Gain of H3K27ac", "Gain of H3K9me3", "Gain of both", "Loss of H3K27ac", "Loss of H3K9me3"))
epigA<-ggplot(epigResults[epigResults$change=="Loss of H3K27ac" | epigResults$change=="Loss of H3K9me3" | epigResults$change=="Gain of H3K27ac" | epigResults$change=="Gain of H3K9me3" | epigResults$change=="Gain of both",], aes(x=TEpresence, y =log2(avgExpr))) + geom_boxplot(aes(fill=TEpresence)) + facet_grid(.~change) + coord_flip() +
  theme_Publication() + labs(y=expression(log[2](TMM)),x="",fill="") +  theme(panel.border = element_blank(), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_fill_manual(values=c("gray80", "firebrick3"), labels = c("Transcripts TE(-)","Transcripts TE(+)"))  + theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) +   scale_x_discrete(limits=rev)  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16))
brks <- c(0, 0.25, 0.5, 0.75, 1)
epiGB<-ggplot(dfFC, aes(x=change,y=value,fill=variable)) + geom_col(position = "fill") + theme_Publication() + scale_fill_manual(values=c("cadetblue3","brown3","gray80"),breaks = c("FC_low", "FC_high"), labels = c("FC<1", "FC>1"),na.value="gray80") +labs(x="", y ="Percentage", fill="Fold-change") +   scale_y_continuous(breaks = brks, labels = scales::percent(brks))  + theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()) + geom_text(data=subset(dfFC, value!=0), aes( label = value), vjust = 1.5, position="fill")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16))
fig6<-plot_grid(epigA,epiGB,ncol=1,rel_heights = c(0.5,0.75), align = 'v', axis = 'l', labels="AUTO")
ggsave(filename = "FigureS1.png", fig6, width = 10, height = 6)
ggsave(filename = "FigureS1.pdf", fig6, width = 10, height = 6)
