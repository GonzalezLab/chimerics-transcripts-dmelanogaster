
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

df.expression.comparison.figure2_original <- df.expression.comparison.figure2
plotA<-ggplot(unique(df.expression.comparison.figure2), aes(x=status, y = log2(avgExpr))) + geom_boxplot(fill=c("#FF6A6A", "#FF6A6A", "firebrick3", "gray80")) + theme_Publication() + labs(y=expression(log[2](TMM)),x="")  + coord_flip() +  theme(panel.border = element_blank(), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))   + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16))

TE_chim_global_expr_v2$position <- ''
TE_chim_global_expr_v2[TE_chim_global_expr_v2$exonPosDescription=="First exon" | TE_chim_global_expr_v2$exonPosDescription=="First exon (TE inside)" ,]$position<-"5'UTR"
TE_chim_global_expr_v2[TE_chim_global_expr_v2$exonPosDescription=="Last exon" | TE_chim_global_expr_v2$exonPosDescription=="Last exon (TE inside)" ,]$position<-"3'UTR"
TE_chim_global_expr_v2[TE_chim_global_expr_v2$exonPosDescription=="Middle exon" | TE_chim_global_expr_v2$exonPosDescription=="Middle exon (TE inside)" ,]$position<-"Internal exons"
TE_chim_global_expr_v2_subset <- unique(TE_chim_global_expr_v2[,c("stringtieID", "FlyBaseGeneRef","position","assembly", "tissue","avgExpr")])
TE_chim_global_expr_v2_subset<-TE_chim_global_expr_v2_subset[TE_chim_global_expr_v2_subset$avgExpr>=1,]
colnames(TE_chim_global_expr_v2_subset) <- c("stringtieID", "FlyBaseGeneRef","status","assembly", "tissue","avgExpr")
df.expression.comparison.figure2<-rbind(df.expression.comparison.figure2,TE_chim_global_expr_v2_subset)
df.expression.comparison.figure2<-unique(df.expression.comparison.figure2)

df.expression.comparison.figure2$status <- factor(df.expression.comparison.figure2$status,levels=rev(c("Non-chimeric transcripts", "Chimeric transcripts", "Overlap and AS insertions", "Internal insertions", "5'UTR", "Internal exons", "3'UTR")))
df.expression.comparison.figure2_chim_noRoo_noChim$type<-NULL
df.expression.comparison.figure2 <- rbind(df.expression.comparison.figure2,df.expression.comparison.figure2_chim_noRoo_noChim[df.expression.comparison.figure2_chim_noRoo_noChim$status=="Chimeric transcripts (without roo)",])

df.expression.comparison.figure2$status <- factor(df.expression.comparison.figure2$status,levels=rev(c("Non-chimeric transcripts", "Chimeric transcripts", "Chimeric transcripts (without roo)", "Overlap and AS insertions", "Internal insertions", "5'UTR", "Internal exons", "3'UTR")))

plotA<-ggplot(unique(df.expression.comparison.figure2), aes(x=status, y = log2(avgExpr))) + geom_boxplot(fill=c("cadetblue2","cadetblue2","cadetblue2","#FF6A6A", "#FF6A6A","firebrick3",  "firebrick3", "gray80")) + theme_Publication() + labs(y=expression(log[2](TMM)),x="")  + coord_flip() +  theme(panel.border = element_blank(), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))   + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16))

newPlot <- df.expression.comparison.figure2[df.expression.comparison.figure2$status=="Non-chimeric transcripts" | df.expression.comparison.figure2$status=="5'UTR" | df.expression.comparison.figure2$status=="3'UTR" | df.expression.comparison.figure2$status=="Internal exons",]

newPlot$status <- factor(newPlot$status, levels=rev(c("Non-chimeric transcripts", "5'UTR", "Internal exons", "3'UTR")))

ggplot(newPlot, aes(x=status, y = log2(avgExpr))) + geom_boxplot(fill=c("#FF6A6A", "#FF6A6A", "#FF6A6A", "gray80")) + theme_Publication() + labs(y=expression(log[2](TMM)),x="")  + coord_flip() +  theme(panel.border = element_blank(), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))   + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16))

wilcox.test(df.expression.comparison.figure2[df.expression.comparison.figure2$status=="Non-chimeric transcripts",]$avgExpr,df.expression.comparison.figure2[df.expression.comparison.figure2$status=="Chimeric transcripts",]$avgExpr)
wilcox.test(df.expression.comparison.figure2[df.expression.comparison.figure2$status=="Non-chimeric transcripts",]$avgExpr,df.expression.comparison.figure2[df.expression.comparison.figure2$status=="Overlap and AS insertions",]$avgExpr)
wilcox.test(df.expression.comparison.figure2[df.expression.comparison.figure2$status=="Non-chimeric transcripts",]$avgExpr,df.expression.comparison.figure2[df.expression.comparison.figure2$status=="Internal insertions",]$avgExpr)
wilcox.test(df.expression.comparison.figure2[df.expression.comparison.figure2$status=="Non-chimeric transcripts",]$avgExpr,df.expression.comparison.figure2[df.expression.comparison.figure2$status=="5'UTR",]$avgExpr)
wilcox.test(df.expression.comparison.figure2[df.expression.comparison.figure2$status=="Non-chimeric transcripts",]$avgExpr,df.expression.comparison.figure2[df.expression.comparison.figure2$status=="3'UTR",]$avgExpr)
wilcox.test(df.expression.comparison.figure2[df.expression.comparison.figure2$status=="Non-chimeric transcripts",]$avgExpr,df.expression.comparison.figure2[df.expression.comparison.figure2$status=="Internal exons",]$avgExpr)
wilcox.test(df.expression.comparison.figure2[df.expression.comparison.figure2$status=="Internal exons",]$avgExpr,df.expression.comparison.figure2[df.expression.comparison.figure2$status=="3'UTR",]$avgExpr)
wilcox.test(df.expression.comparison.figure2[df.expression.comparison.figure2$status=="5'UTR",]$avgExpr,df.expression.comparison.figure2[df.expression.comparison.figure2$status=="3'UTR",]$avgExpr)

ggsave("FigS3.png", width = 12, height = 2)
ggsave("FigS3.pdf", width = 12, height = 2)

colors <- c(rep("#009ACD",29), "#FFA500")

dfType <- data.frame(matrix(ncol =2, nrow = 0))
x <- c("gene", "type")
colnames(dfType) <- x


for (gene in averageGeneExpressionContribution$gene) {
  result<-sum(unique(TE_chim_global[TE_chim_global$FlyBaseGeneRef==gene,]$group))
  
  if (result == 3) {
    df <- data.frame(gene=gene,type="Both groups")
    dfType <- rbind(dfType,df)
  } else if (result == 2) {
    df <- data.frame(gene=gene,type="Group 2")
    dfType <- rbind(dfType,df)
  } else if (result == 1) {
    df <- data.frame(gene=gene,type="Group 1")
    dfType <- rbind(dfType,df)
  }
}
View(dfType)

averageGeneExpressionContribution_original <- averageGeneExpressionContribution
averageGeneExpressionContribution <- merge(averageGeneExpressionContribution,dfType,by="gene")

averageGeneExpressionContribution$color<-""
averageGeneExpressionContribution[averageGeneExpressionContribution$type=="Both groups" & averageGeneExpressionContribution$Mean < 100, ]$color<-"#00A3D9"
averageGeneExpressionContribution[averageGeneExpressionContribution$type=="Group 2" & averageGeneExpressionContribution$Mean < 100, ]$color<-"#0086B3"
averageGeneExpressionContribution[averageGeneExpressionContribution$type=="Group 1" & averageGeneExpressionContribution$Mean < 100, ]$color<-"#00394D"
averageGeneExpressionContribution[averageGeneExpressionContribution$type=="Both groups" & averageGeneExpressionContribution$Mean == 100, ]$color<-"#FFA500"
averageGeneExpressionContribution[averageGeneExpressionContribution$type=="Group 2" & averageGeneExpressionContribution$Mean == 100, ]$color<-"#805300"
averageGeneExpressionContribution[averageGeneExpressionContribution$type=="Group 1" & averageGeneExpressionContribution$Mean == 100, ]$color<-"#402900"


averageGeneExpressionContribution$type<-factor(averageGeneExpressionContribution$type, levels=c("Group 1", "Group 2", "Both groups"))
averageGeneExpressionContribution$color<-factor(averageGeneExpressionContribution$color)

averageGeneExpressionContribution$type2<-""
averageGeneExpressionContribution[averageGeneExpressionContribution$type=="Both groups" & averageGeneExpressionContribution$Mean < 100, ]$type2<-"Both groups (variable)"
averageGeneExpressionContribution[averageGeneExpressionContribution$type=="Group 2" & averageGeneExpressionContribution$Mean < 100, ]$type2<-"Group 2 "
averageGeneExpressionContribution[averageGeneExpressionContribution$type=="Group 1" & averageGeneExpressionContribution$Mean < 100, ]$type2<-"Group 1 "
averageGeneExpressionContribution[averageGeneExpressionContribution$type=="Both groups" & averageGeneExpressionContribution$Mean == 100, ]$type2<-"Both groups (always chimeric)"
averageGeneExpressionContribution[averageGeneExpressionContribution$type=="Group 2" & averageGeneExpressionContribution$Mean == 100, ]$type2<-"Group 2"
averageGeneExpressionContribution[averageGeneExpressionContribution$type=="Group 1" & averageGeneExpressionContribution$Mean == 100, ]$type2<-"Group 1"
averageGeneExpressionContribution$type2<-factor(averageGeneExpressionContribution$type2)


plotB<-ggplot(averageGeneExpressionContribution,aes(x=Mean,fill=color)) + geom_histogram( boundary = 0, closed = "left") + theme_Publication() + scale_fill_identity(guide = 'legend', labels = c("Overlap & AS insertions", "Overlap & AS insertions", "Internal insertions", "Internal insertions", "Both groups (variable)", "Both groups (always chimeric)"), breaks = c("#00394D","#402900","#0086B3","#805300","#00A3D9","#FFA500"))    + geom_vline(xintercept=median(averageGeneExpressionContribution[averageGeneExpressionContribution$Mean!=100,]$Mean), linetype = "longdash") +  scale_y_continuous( expand = c(0, 0)) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + labs(x="Mean chimeric gene-TE transcript expression contribution", y = "Number of chimeric genes", fill="Type") + geom_text(aes(x=median(averageGeneExpressionContribution[averageGeneExpressionContribution$Mean!=100,]$Mean),  y=200), label=paste0("\nMedian = ",round(median(averageGeneExpressionContribution[averageGeneExpressionContribution$Mean!=100,]$Mean),1),"%"), angle=90, check_overlap = TRUE)   + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16))

plot_grid(plotA, plotB, align = "h", axis = "t", nrow=2,labels = c('A', 'B'), label_size = 18, rel_heights =  c(1.5,3))

ggsave(plot_grid(plotA, plotB, align = "h", axis = "t", nrow=2,labels = c('A', 'B'), label_size = 18, rel_heights =  c(1.5,3))
       ,file="Figure5_def2.png", width = 7.8, height = 8)

ggsave(plot_grid(plotA, plotB, align = "h", axis = "t", nrow=2,labels = c('A', 'B'), label_size = 18, rel_heights =  c(1.5,3))
       ,file="Figure5_def2.pdf", width = 7.8, height = 8)



# contribution by lengthTE to expression
### are TEs with shorter length TE the ones that have expression 100%?

TE_chim_length_TE_gene <- unique(TE_chim_global[c("FlyBaseGeneRef","lengthTE","TEconsensusID")])


dfType <- data.frame(matrix(ncol =2, nrow = 0))
x <- c("gene", "size")
colnames(dfType) <- x

for ( gene in unique(TE_chim_length_TE_gene$FlyBaseGeneRef) ) {
  n=length(unique(TE_chim_length_TE_gene[TE_chim_length_TE_gene$FlyBaseGeneRef==gene,]$lengthTE >= 120))
  if ( n == 2 ) {
    print(paste(gene,"mix"))
    df <- data.frame(gene=gene,size="mix")
    dfType <- rbind(dfType,df)
  }
  else if (n == 1) {
    type=unique(TE_chim_length_TE_gene[TE_chim_length_TE_gene$FlyBaseGeneRef==gene,]$lengthTE >= 120)
    if ( type == "TRUE" ) {
      print(paste(gene,"big"))
      df <- data.frame(gene=gene,size="big")
      dfType <- rbind(dfType,df)
    } else if (type == "FALSE") {
      print(paste(gene,"short"))
      df <- data.frame(gene=gene,size="short")
      dfType <- rbind(dfType,df)
    }
  }
}

averageGeneExpressionContribution_size <- merge(averageGeneExpressionContribution,dfType,by="gene")

df.expression.comparison.figure2

df.expression.chimeric <- df.expression.comparison.figure2[df.expression.comparison.figure2$status=="Chimeric transcripts",]
df.expression.chimeric<-merge(df.expression.chimeric,averageGeneExpressionContribution_size,by.x = "FlyBaseGeneRef",by.y="gene")

df.expression.chimeric[df.expression.chimeric$size=="short",]$size <- "Short insertions"
df.expression.chimeric[df.expression.chimeric$size=="big",]$size <- " ≥120bp insertions"

ggplot(df.expression.chimeric[df.expression.chimeric$size!="mix",], aes(x=size, y = log2(avgExpr))) + geom_boxplot(fill=c("#FF6A6A", "#FF6A6A")) + theme_Publication() + labs(y=expression(log[2](TMM)),x="")  + coord_flip() +  theme(panel.border = element_blank(), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))   + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16)) 

ggsave("FigS4.pdf", height = 2, width = 10)
ggsave("FigS4.png", height = 2, width = 10)


> length(unique(averageGeneExpressionContribution_size[averageGeneExpressionContribution_size$type=="Group 2" & averageGeneExpressionContribution_size$size=="big"  ,]$gene))

resultsRoo
df.expression.comparison.figure2_original
df.expression.comparison.figure2_chim <- df.expression.comparison.figure2[df.expression.comparison.figure2$status=="Chimeric transcripts",]
df.expression.comparison.figure2_noChim <- df.expression.comparison.figure2[df.expression.comparison.figure2$status=="Non-chimeric transcripts",]
df.expression.comparison.figure2_chim
resultsRoo <- unique(resultsRoo[,c("tissue","assembly","stringtieID","FlyBaseGeneRef","type")])
df.expression.comparison.figure2_chimRoo<-merge(df.expression.comparison.figure2_chim,resultsRoo,by=c("tissue","assembly","stringtieID","FlyBaseGeneRef"),all.x = TRUE)   

df.expression.comparison.figure2_chimRoo[is.na(df.expression.comparison.figure2_chimRoo$type),]$type<-"No roo"
df.expression.comparison.figure2_chimnoRoo<-df.expression.comparison.figure2_chimRoo[df.expression.comparison.figure2_chimRoo$type!="Repeat",]

df.expression.comparison.figure2_noChim$type <- "No roo"
df.expression.comparison.figure2_chim_noRoo_noChim <- rbind(df.expression.comparison.figure2_noChim,df.expression.comparison.figure2_chimnoRoo)
t.test(df.expression.comparison.figure2_chim_noRoo_noChim[df.expression.comparison.figure2_chim_noRoo_noChim$status=="Non-chimeric transcripts",]$avgExpr,df.expression.comparison.figure2_chim_noRoo_noChim[df.expression.comparison.figure2_chim_noRoo_noChim$status=="Chimeric transcripts",]$avgExpr)
wilcox.test(df.expression.comparison.figure2_chim_noRoo_noChim[df.expression.comparison.figure2_chim_noRoo_noChim$status=="Non-chimeric transcripts",]$avgExpr,df.expression.comparison.figure2_chim_noRoo_noChim[df.expression.comparison.figure2_chim_noRoo_noChim$status=="Chimeric transcripts",]$avgExpr)

df.expression.comparison.figure2_chim_noRoo_noChim[df.expression.comparison.figure2_chim_noRoo_noChim$status=="Chimeric transcripts",]$status <- "Chimeric transcripts (without roo)"
df.expression.comparison.figure2_chim_noRoo_noChim$status <- factor(df.expression.comparison.figure2_chim_noRoo_noChim$status, levels = rev(c("Non-chimeric transcripts", "Chimeric transcripts (without roo)")))
ggplot(df.expression.comparison.figure2_chim_noRoo_noChim, aes(x=status, y = log2(avgExpr))) + geom_boxplot(fill=c("#FF6A6A", "#FF6A6A")) + theme_Publication() + labs(y=expression(log[2](TMM)),x="")  + coord_flip() +  theme(panel.border = element_blank(), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))   + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16)) 

ggsave("figure.png", width = 12, height = 2)
ggsave("figure.pdf", width = 12, height = 2)

wilcox.test(df.expression.comparison.figure2_chim_noRoo_noChim[df.expression.comparison.figure2_chim_noRoo_noChim$status=="Non-chimeric transcripts",]$avgExpr,df.expression.comparison.figure2_chim_noRoo_noChim[df.expression.comparison.figure2_chim_noRoo_noChim$status=="Chimeric transcripts (without roo)",]$avgExpr)


averageGeneExpressionContribution_size_noRoo<-averageGeneExpressionContribution_size[averageGeneExpressionContribution_size$gene%in%df.expression.comparison.figure2_chim_noRoo_noChim[df.expression.comparison.figure2_chim_noRoo_noChim$status=="Chimeric transcripts (without roo)",]$FlyBaseGeneRef,]
median(averageGeneExpressionContribution_size_noRoo$Mean)
median(averageGeneExpressionContribution_size_noRoo[averageGeneExpressionContribution_size_noRoo$Mean<100,]$Mean)
