library(ggplot2)
library(ggpubr)
library(VennDiagram)
library(wesanderson)
library(data.table)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggvenn)

listTranscriptsTissueTotalCommonUniqueByCategories <- fread("../../resultsChimAnalysis/tableS3_fig1.tsv", header = T)

listTranscriptsTissueTotalCommonUniqueByCategories[listTranscriptsTissueTotalCommonUniqueByCategories$category=="all",]$category <- "All (2169)"
listTranscriptsTissueTotalCommonUniqueByCategories[listTranscriptsTissueTotalCommonUniqueByCategories$category=="g1",]$category <- "Overlap & AS\ninsertions (977)"
listTranscriptsTissueTotalCommonUniqueByCategories[listTranscriptsTissueTotalCommonUniqueByCategories$category=="g2",]$category <- "Internal\ninsertions (1587)"
listTranscriptsTissueTotalCommonUniqueByCategories[listTranscriptsTissueTotalCommonUniqueByCategories$analysis=="bodypart",]$analysis <- "Shared chimeric gene-TE transcripts across body parts"
listTranscriptsTissueTotalCommonUniqueByCategories[listTranscriptsTissueTotalCommonUniqueByCategories$analysis=="strain",]$analysis <- "Shared chimeric gene-TE transcripts across strains"

listTranscriptsTissueTotalCommonUniqueByCategories[listTranscriptsTissueTotalCommonUniqueByCategories$type=="specific",]$type <- "Body part-specific"
listTranscriptsTissueTotalCommonUniqueByCategories[listTranscriptsTissueTotalCommonUniqueByCategories$type=="unique",]$type <- "Strain-specific "
listTranscriptsTissueTotalCommonUniqueByCategories[listTranscriptsTissueTotalCommonUniqueByCategories$type=="shared2",]$type <- "Shared across 2 body parts"
listTranscriptsTissueTotalCommonUniqueByCategories[listTranscriptsTissueTotalCommonUniqueByCategories$type=="shared3",]$type <- "Shared across 3 body parts"
listTranscriptsTissueTotalCommonUniqueByCategories[listTranscriptsTissueTotalCommonUniqueByCategories$type=="pol",]$type <- "Shared across 2-4 strains"
listTranscriptsTissueTotalCommonUniqueByCategories[listTranscriptsTissueTotalCommonUniqueByCategories$type=="shared5",]$type <- "Shared across 5 strains"

listTranscriptsTissueTotalCommonUniqueByCategoriesPer<-listTranscriptsTissueTotalCommonUniqueByCategories%>%group_by(analysis, category) %>% 
  mutate(sum = sum(number)) %>% mutate(percentage = number/sum)

listTranscriptsTissueTotalCommonUniqueByCategoriesPer$category<-factor(listTranscriptsTissueTotalCommonUniqueByCategoriesPer$category,levels=c("All (2169)", "Overlap & AS\ninsertions (977)", "Internal\ninsertions (1587)"))
listTranscriptsTissueTotalCommonUniqueByCategoriesPer$type<-factor(listTranscriptsTissueTotalCommonUniqueByCategoriesPer$type,levels=rev(c("Body part-specific", "Strain-specific ", "Shared across 2 body parts", "Shared across 2-4 strains", "Shared across 3 body parts", "Shared across 5 strains")))


p1c<-ggplot(listTranscriptsTissueTotalCommonUniqueByCategoriesPer[listTranscriptsTissueTotalCommonUniqueByCategoriesPer$analysis=="Shared chimeric gene-TE transcripts across strains",] , aes(x=category, y = percentage, fill = type)) + geom_col(alpha=0.8) + theme_Publication() + scale_fill_Publication() + labs(y="Percentage of chimeric\ngene-TE transcript", x="", fill = "Type")  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16),plot.title = element_text(size = 16)) +
  geom_text(aes(label=number), position=position_stack(),vjust = 1.2, size = 4) + scale_y_continuous(labels=scales::percent,expand = c(0,0)) +
  guides(fill = guide_legend(reverse=TRUE)) + ggtitle("Chimeric gene-TE transcripts across strains")

p1d<-ggplot(listTranscriptsTissueTotalCommonUniqueByCategoriesPer[listTranscriptsTissueTotalCommonUniqueByCategoriesPer$analysis=="Shared chimeric gene-TE transcripts across body parts",] , aes(x=category, y = percentage, fill = type))  + geom_col(alpha=0.8) + theme_Publication() + scale_fill_manual(values=c("#ef3b2c","#662506","#a6cee3")) + labs(y="Percentage of chimeric\ngene-TE transcript", x="", fill = "Type")  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16)) +
  geom_text(aes(label=number), position=position_stack(),vjust = 1.2, size = 4) + scale_y_continuous(labels=scales::percent,expand = c(0,0)) +
  guides(fill = guide_legend(reverse=TRUE))  + ggtitle("Chimeric gene-TE transcripts across body parts")

ggarrange(p1c,p1d, ncol = 1, labels = c("C", "D"))
ggarrange(p1c,p1d, nrow = 1, labels = c("A", "B"))

ggsave(p1c, filename="Figure1C.png", width = 8, height = 5.25)
ggsave(p1d, filename="Figure1D.png", width = 8, height = 5.25)

ggsave(ggarrange(p1c,p1d, nrow = 1, labels = c("A", "B")), filename = "FigureS1.png",  width = 8, height = 5.25)
ggsave(ggarrange(p1c+guides(fill=guide_legend(ncol =1)),p1d+guides(fill=guide_legend(ncol =1)), nrow = 1, labels = c("A", "B"),  font.label = list(size = 18)), filename = "FigureS1.png",  width = 12.5, height = 5.25)

ggsave(ggarrange(p1c+guides(fill=guide_legend(ncol =1)),p1d+guides(fill=guide_legend(ncol =1)), nrow = 1, labels = c("A", "B"),  font.label = list(size = 18)), filename = "FigureS1.pdf",  width = 12.5, height = 5.25)


kk<-fread("tableClimates.tsv")
kkP<-ggplot(kk, aes(x = Climates, y=value, fill=Climates)) + geom_col() + scale_fill_manual(values=kk$col) + theme_Publication() +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
ggsave(kkP, filename="Figure1A_legend.png", width = 8, height = 5.25)
