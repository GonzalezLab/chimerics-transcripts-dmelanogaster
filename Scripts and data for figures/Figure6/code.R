
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

dataGO_BP[dataGO_BP$shortName=="Regulation, metabolic process",]$shortName <- "Regulation of metabolism"
dataGO_BP[dataGO_BP$shortName=="Metabolism" & dataGO_BP$tissue=="Ovary",]$type <- "Regulation"
dataGO_BP[dataGO_BP$shortName=="Metabolism" & dataGO_BP$tissue=="Ovary",]$shortName <- "Regulation of metabolism"
dataGO_BP[dataGO_BP$shortName=="Regulation of metabolism",]$type <- "Metabolism"
dataGO_BP[dataGO_BP$shortName=="Positive regulation of RNA transcription, metabolism",]$type <- "Metabolism"

dataGO_BP_g1 <- fread("../../resultsChimAnalysis/GO/GO_tissues/GO_tissues_byGroup/resultParseBP_g1.tsv")
dataGO_BP_g2 <- fread("../../resultsChimAnalysis/GO/GO_tissues/GO_tissues_byGroup/resultParseBP_g2.tsv")

map_color <- c(Development = "#E41A1C", Regulation = "#377EB8", Metabolism = "#4DAF4A", 
  Component = "#984EA3", "Transcription")

dataGO_BP$enrichmentScore<-as.numeric(as.character(dataGO_BP$enrichmentScore))
dataGO_BP$classification <- "All"
dataGO_BP$shortName<-factor(dataGO_BP$shortName, levels=rev(c("Regulation of metabolism", "Metabolism",   "Development, morphogenesis",   "Positive regulation of RNA transcription, metabolism",  "Organ development", "Regulation, response to stimulus, signaling", "Anatomical structure development", "Regulation, signaling, communication", "Cellular component organization")))
dataGO_BP<-dataGO_BP %>% arrange(factor(shortName, levels = rev(c("Regulation of metabolism", "Metabolism",   "Development, morphogenesis",   "Positive regulation of RNA transcription, metabolism",  "Organ development", "Regulation, response to stimulus, signaling", "Anatomical structure development", "Regulation, signaling, communication", "Cellular component organization"))))
dataGO_BP$shortName<-as.character(dataGO_BP$shortName)

dataGO_BP_g1$enrichmentScore<-as.numeric(as.character(dataGO_BP_g1$enrichmentScore))
dataGO_BP_g1$classification <- "Overlap & AS\ninsertions"
dataGO_BP_g1$shortName <- factor(dataGO_BP_g1$shortName, levels=rev(c("Cellular component organization","Nucleosome assembly and organization", "Cilium assembly and organization")))
dataGO_BP_g1<-dataGO_BP_g1 %>% arrange(factor(shortName, levels = rev(c("Cellular component organization","Nucleosome assembly and organization", "Cilium assembly and organization"))))

dataGO_BP_g2$enrichmentScore<-as.numeric(as.character(dataGO_BP_g2$enrichmentScore))
dataGO_BP_g2$classification <- "Internal\ninsertions"
dataGO_BP_g2$shortName<-factor(dataGO_BP_g2$shortName, levels=rev(c("Organ development", "Regulation of metabolism",  "Development, morphogenesis","Metabolism", "Positive regulation of metabolism", "Epithelium development", "Tissue development", "Cell differentiation and development",   "Transcription", "Regulation, development", "Regulation")))
dataGO_BP_g2<-dataGO_BP_g2 %>% arrange(factor(shortName, levels = rev(c("Organ development", "Regulation of metabolism",  "Development, morphogenesis","Metabolism", "Positive regulation of metabolism", "Epithelium development", "Tissue development", "Cell differentiation and development",   "Transcription", "Regulation, development", "Regulation"))))
dataGO_BP_g2[dataGO_BP_g2$shortName=="Transcription",]$type <- "Transcription"

dataGO_BP$classification<-"All chimeric\ngene-TE transcripts"  
map_color <- c("Development" = "#E41A1C", "Regulation" = "#377EB8", "Metabolism" = "#4DAF4A", "Component" = "#984EA3", "Transcription" = "#ff7f00")
GO_BP_plot <-ggplot(dataGO_BP, aes(y=enrichmentScore, x=shortName)) + geom_col(aes(fill=type),alpha=0.7) + coord_flip() + facet_grid(classification~ tissue, scales = "free", space = "free") + labs(x="", y = "Enrichment score", fill="Annotation cluster") + theme_Publication() + geom_text(aes(label = nGenes), hjust = 0, size=5)  + scale_fill_manual(values=map_color) +
  theme(axis.text = element_text( size = 16)) + scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2)) + guides(fill="none")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + theme (axis.line.x = element_line(color = "white", 
                                                                          size = 1, linetype = "solid"))
GO_BP_plot
dataGO_BP_g1$classification <- "Overlap\n & AS"

GO_BP_g1_plot <-ggplot(dataGO_BP_g1, aes(y=enrichmentScore, x=shortName)) + geom_col(aes(fill=type),alpha=0.7) + coord_flip() + facet_grid(classification~ tissue, scales = "free", space = "free") + labs(x="", y = "Enrichment score", fill="Annotation cluster") + theme_Publication() + geom_text(aes(label = nGenes), hjust = 0, size=5)  + scale_fill_manual(values=map_color) +
  theme(axis.text = element_text( size = 16)) + scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2)) + guides(fill="none")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16)) + 
  theme(strip.background.x = element_blank(),strip.text.x = element_blank()) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + theme (axis.line.x = element_line(color = "white", 
                                                                          size = 1, linetype = "solid"))
GO_BP_g1_plot
dataGO_BP_g2$tissue <- factor(dataGO_BP_g2$tissue, levels=c("Head", "Gut", "Ovary"))
dataGO_BP_g2[dataGO_BP_g2$shortName=="Regulation of metabolism",]$type <- "Metabolism"
dataGO_BP_g2[dataGO_BP_g2$shortName=="Positive regulation of metabolism",]$type <- "Metabolism"
dataGO_BP_g2[dataGO_BP_g2$shortName=="Transcription",]$type <- "Transcription"

GO_BP_g2_plot <-ggplot(dataGO_BP_g2, aes(y=enrichmentScore, x=shortName)) + geom_col(aes(fill=type),alpha=0.7) + coord_flip() + facet_grid(classification~ tissue, scales = "free", space = "free") + labs(x="", y = "Enrichment score", fill="Annotation cluster") + theme_Publication() + geom_text(aes(label = nGenes), hjust = 0, size=5)  + scale_fill_manual(values=map_color) +
  theme(axis.text = element_text( size = 16)) + scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2)) + guides(fill="none")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16))  + 
  theme(strip.background.x = element_blank(),strip.text.x = element_blank())
GO_BP_g2_plot

plot_A <- plot_grid(GO_BP_plot, GO_BP_g1_plot, GO_BP_g2_plot, nrow=3, align = 'v', axis = 'l', rel_heights = c(0.8,0.3,1))

color_map_mf<-c(`DNA binding` = "#fee08b", `Transmembrane transporter activity` = "#99d594", 
  `Enzyme activator activity` = "#A65628", Transcription = "#999999", 
  `Ion binding` = "#F781BF")

dataGO_MF$enrichmentScore<-as.numeric(as.character(dataGO_MF$enrichmentScore)) 
dataGO_MF$classification <- "All chimeric\ngene-TE transcripts"
dataGO_MF[dataGO_MF$shortName=="Ion channel activity",]$shortName <- "Transmembrane transporter activty"
dataGO_MF$shortName <- factor(dataGO_MF$shortName, levels=rev(c("Sequence-specific DNA binding", "RNA polymerase II transcription", "GTPase activator", "Nucleic acid binding" , "Ion binding", "Transcription regulatory region DNA binding","Transmembrane transporter activty"  , "Transcription factor activity"       )))
dataGO_MF_g1 <- fread("../../resultsChimAnalysis/GO/GO_tissues/GO_tissues_byGroup/resultParseMF_g1.tsv")
dataGO_MF_g2 <- fread("../../resultsChimAnalysis/GO/GO_tissues/GO_tissues_byGroup/resultParseMF_g2.tsv")

dataGO_MF_g1$enrichmentScore<-as.numeric(as.character(dataGO_MF_g1$enrichmentScore))
dataGO_MF_g1$classification <- "Overlap & AS\ninsertions"
dataGO_MF_g2$enrichmentScore<-as.numeric(as.character(dataGO_MF_g2$enrichmentScore))
dataGO_MF_g2$classification <- "Internal\ninsertions"

dataGO_MF_g2[dataGO_MF_g2$shortName=="Ion transport",]$shortName <- "Transmembrane transporter activty"
dataGO_MF_g2[dataGO_MF_g2$type =="Ion transport",]$type <- "Transmembrane transporter activity"

dataGO_MF_g2$shortName<-factor(dataGO_MF_g2$shortName, levels=rev(c("Nucleic acid binding", "Transcription factor activity" , "Activator activity" , "RNA polymerase II transcription", "Ion transport" , "Ion binding" )))
dataGO_MF_g2<-dataGO_MF_g2 %>% arrange(factor(shortName, levels=rev(c("Nucleic acid binding", "Transcription factor activity" , "Activator activity" , "RNA polymerase II transcription", "Ion transport" , "Ion binding" ))))

GO_MF_plot <-ggplot(dataGO_MF, aes(y=enrichmentScore, x=shortName)) + geom_col(aes(fill=type),alpha=0.7) + coord_flip() + facet_grid(classification~ tissue, scales = "free", space = "free") + labs(x="", y = "Enrichment score", fill="Annotation cluster") + theme_Publication() + geom_text(aes(label = nGenes), hjust = 0, size=5)  + scale_fill_manual(values=color_map_mf) +
  theme(axis.text = element_text( size = 16)) + scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2)) + guides(fill="none")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + theme (axis.line.x = element_line(color = "white", 
                                                                          size = 1, linetype = "solid"))
GO_MF_plot
dataGO_MF_g1$tissue <- factor(dataGO_MF_g1$tissue, levels=c("Head" ,"Gut","Ovary"))
dataGO_MF_g1$classification <-  "Overlap\n & AS"
GO_MF_g1_plot <-ggplot(dataGO_MF_g1, aes(y=enrichmentScore, x=shortName)) + geom_col(aes(fill=type),alpha=0.7) + coord_flip() + facet_grid(classification~ tissue, scales = "free", space = "free") + labs(x="", y = "Enrichment score", fill="Annotation cluster") + theme_Publication() + geom_text(aes(label = nGenes), hjust = 0, size=5)  + scale_fill_manual(values=color_map_mf) + 
  theme(axis.text = element_text( size = 16)) + scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2)) + guides(fill="none")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16)) + 
  theme(strip.background.x = element_blank(),strip.text.x = element_blank()) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + theme (axis.line.x = element_line(color = "white", 
                                                                          size = 1, linetype = "solid"))
GO_MF_g1_plot
dataGO_MF_g2$tissue <- factor(dataGO_MF_g2$tissue, levels=c("Head" ,"Gut","Ovary"))
GO_MF_g2_plot <-ggplot(dataGO_MF_g2, aes(y=enrichmentScore, x=shortName)) + geom_col(aes(fill=type),alpha=0.7) + coord_flip() + facet_grid(classification~ tissue, scales = "free", space = "free") + labs(x="", y = "Enrichment score", fill="Annotation cluster") + theme_Publication() + geom_text(aes(label = nGenes), hjust = 0, size=5)  + scale_fill_manual(values=color_map_mf) + 
  theme(axis.text = element_text( size = 16)) + scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2)) + guides(fill="none")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16))  + 
  theme(strip.background.x = element_blank(),strip.text.x = element_blank())
GO_MF_g2_plot

plot_grid(GO_BP_plot, GO_BP_g1_plot, GO_BP_g2_plot, GO_MF_plot, GO_MF_g1_plot, GO_MF_g2_plot, ncol=1, align = 'v', axis = 'l', labels=c("A","","","B","",""), label_size=18, rel_heights = c(1.2, 0.4,1.3,1,0.4,1))

ggsave("Figure6_3.png", height = 12, width = 13)
ggsave("Figure6_3.pdf", height = 12, width = 13)


dataGO_MF_g1$enrichmentScore<-as.numeric(as.character(dataGO_MF_g1$enrichmentScore)) 
dataGO_MF_g1$classification <- "Overlap & AS\ninsertions"
dataGO_MF_g2$enrichmentScore<-as.numeric(as.character(dataGO_MF_g2$enrichmentScore)) 
dataGO_MF_g2$classification <- "Internal insertions"

dataGO_MF_all <- rbind(dataGO_MF,dataGO_MF_g1,dataGO_MF_g2)


GO_MF_plot <-ggplot(dataGO_MF_all, aes(y=enrichmentScore, x=shortName)) + geom_col(aes(fill=type),alpha=0.7) + coord_flip() + facet_grid(classification~ tissue,scale="free",space="free") + labs(x="", y = "Enrichment score", fill="Annotation cluster") + theme_Publication() + geom_text(aes(label = nGenes), hjust = 0,size=5)  + scale_fill_manual(values=c( "#FF7F00","#FFFF33","#A65628", "#999999","#F781BF", "#F781BF", "#984EA3")) +
  theme(axis.text = element_text( size = 16)) +guides(fill=guide_legend(nrow=2,byrow=TRUE))  + scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2)) + guides(fill="none")   + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16))
GO_MF_plot

plotGO<-plot_grid(GO_BP_plot,GO_MF_plot,ncol=1,labels=c("A","B"),  align = 'v', axis = 'l',label_size = 18)
plotGO
ggsave(filename = "Figure6.png", plotGO, width = 14.3, height = 6.5)
ggsave(filename = "Figure6.pdf", plotGO, width = 14.3, height = 6.5)
