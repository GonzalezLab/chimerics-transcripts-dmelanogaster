library(ggplot2)
library(ggpubr)
library(VennDiagram)
library(wesanderson)
library(data.table)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggvenn)

freqtable <- data.frame(tissue=rep(c("Head", "Gut", "Ovary"),each = 9), group = rep("All",27), frequency=rep(c("Strain-specific", "Shared across 5 strains", "Shared across 2-4 strains"), n=3), position=rep(c("3'UTR", "5'UTR", "Internal\nexon"), each = 3), n = c(298,19,124,211,9,68,181,20,119,202,20,80,124,2,49,165,18,125,173,16,100,93,3,32,103,20,122))
freqtable$tissue <- factor(freqtable$tissue, levels = c("Head", "Gut", "Ovary"))
freqtable$frequency <- factor(freqtable$frequency, levels = c("Strain-specific", "Shared across 2-4 strains", "Shared across 5 strains"))
freqG1table <- data.frame(tissue=rep(c("Head", "Gut", "Ovary"),each = 9), group = rep("Overlap &\nAS insertions",27), frequency=rep(c("Strain-specific", "Shared across 5 strains", "Shared across 2-4 strains"), n=3), position=rep(c("3'UTR", "5'UTR", "Internal\nexon"), each = 3), n = c(214,5,58,155,1,26,11,1,7,131,3,37,120,0,22,26,1,4,140,9,41,77,0,18,14,0,10))
freqG2table <- data.frame(tissue=rep(c("Head", "Gut", "Ovary"),each = 9), group = rep("Internal\ninsertions",27), frequency=rep(c("Strain-specific", "Shared across 5 strains", "Shared across 2-4 strains"), n=3), position=rep(c("3'UTR", "5'UTR", "Internal\nexon"), each = 3), n = c(191,18,125,124,8,68,191,20,129,138,21,86,62,2,41,170,17,130,99,15,88,59,3,30,107,20,124))
freqTable<-rbind(freqtable,freqG1table,freqG2table)
freqTable$group <- factor(freqTable$group, levels = c("All", "Overlap &\nAS insertions", "Internal\ninsertions"))
freqTable$tissue <- factor(freqTable$tissue, levels = c("Head", "Gut", "Ovary"))

f3b<-ggplot(freqTable, aes(y=n, x = position)) + geom_boxplot(fill="gray90") + geom_point(aes(color=frequency)) + facet_grid(group~tissue ) + theme_Publication() + scale_colour_Publication()  + labs(color="Frequency", x = "", y = "Number of gene-TE chimeric transcripts") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16))


freq <- data.frame(position=rep(c("5'UTR", "Internal\nexon", "3'UTR"), each = 2), frequency=rep(c("Strain-specific", "Shared across\n2-5 strains"), n=3), n = c(282,100,226,248,392,211))
freq$position<-factor(freq$position,levels=c("5'UTR", "Internal\nexon", "3'UTR"))
freq$frequency<-factor(freq$frequency,levels=c("Strain-specific", "Shared across\n2-5 strains"))

f3a<-ggplot(freq, aes(y=n, x = position)) + geom_col(aes(fill=frequency))   + theme_Publication() + scale_fill_manual(values= c("#386cb0", "#a6cee3")) + labs(fill="Frequency", x = "", y = "Number of gene-TE chimeric transcripts")  +   scale_y_continuous(expand = c(0, 0))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())   + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16))

ggsave(plot_grid(f3a,f3b,  align="h", axis = "b", labels = c("A","B"), label_size = 18,  rel_widths = c(0.45, 1)),file="Figure3.png", width = 15, height = 6.5)
ggsave(plot_grid(f3a,f3b,  align="h", axis = "b", labels = c("A","B"), label_size = 18,  rel_widths = c(0.45, 1)),file="Figure3.pdf", width = 15, height = 6.5)
