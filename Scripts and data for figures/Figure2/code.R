library(ggplot2)
library(ggpubr)
library(VennDiagram)
library(wesanderson)
library(data.table)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggvenn)
library(eulerr)

TE_chim_global_fig1_data <- fread("../../resultsChimAnalysis/tableS3_fig1.tsv", header = T)

f2b<-ggplot(TE_chim_global_fig1_data, aes(x=tissue, y = stringtieID)) + 
  geom_boxplot(fill="gray80") + 
  geom_point(aes(color=assembly)) + 
  facet_grid(.~group ) + 
  theme_Publication() + 
  scale_colour_Publication() +
  scale_y_continuous(breaks = pretty_breaks(8)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16)) + 
  labs(color="Strain", x = "", y = "Number of chimeric gene-TE transcripts")

listTissues = list(`Head`=unique(TE_chim_global[TE_chim_global$tissue=="head",]$stringtieID),
                   `Gut`=unique(TE_chim_global[TE_chim_global$tissue=="gut",]$stringtieID),
                   `Ovary`=unique(TE_chim_global[TE_chim_global$tissue=="ovary",]$stringtieID))


createVD_3g <- function(dataList, filename, category.names, mainTitle, ngroups = 5) {
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  myCol <-  c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506")
  cat.pos = c(-27, 27, 135, 135, 135)
  cat.dist = c(0.055, 0.055, 0.055, 0.055, 0.055)
  col = myCol
  fill = c(alpha(myCol[1],1), alpha(myCol[2],1),
           alpha(myCol[3],1), alpha(myCol[4],1),
           alpha(myCol[5],1))
  
  venn.diagram(
    x = dataList,
    category.names = category.names,
    filename = filename,
    output=TRUE,
    
    # Output features
    imagetype="png" ,
    height = 950, 
    width = 1080, 
    resolution = 300,
    
    # Circles
    lwd = 1,
    # Numbers
    cex = .8,
    fontfamily = "sans",
    # Set names
    cat.cex = 0.8,
    cat.default.pos = "outer",
    cat.fontfamily = "sans",
    cat.pos = cat.pos[1:ngroups],
    cat.dist = cat.dist[1:ngroups],
    col = "black",
    main.fontface = "bold",
    fill = fill[1:ngroups],
    cat.col = "black",  ext.length = 0.5,
    main = mainTitle, main.cex = 0.8, main.fontfamily = "sans"
  )
  print(paste("Shared elements:",length(Reduce(intersect, dataList))))
  print(paste("Total:",length(unique(unlist(dataList)))))
}

createVD_3g(dataList = listTissues, filename = "VD_by_tissue.png", mainTitle = "Shared chimeric gene-TE transcripts across body parts", ngroups = 3, category.names =names(listTissues))

listTissues = list(`Head (1,459)`=unique(TE_chim_global[TE_chim_global$tissue=="head",]$stringtieID),
                   `Gut (1,068)`=unique(TE_chim_global[TE_chim_global$tissue=="gut",]$stringtieID),
                   `Ovary (884)`=unique(TE_chim_global[TE_chim_global$tissue=="ovary",]$stringtieID))

f2a<-plot(euler(listTissues),
          quantities = list(type = c("counts"), cex = 1.2), labels = list(cex = 1.2),
          edges = list(col = "black"),
          fills =list(fill = c("#CD96CD","#D9D9D9","#FFFACD"), alpha = 0.7)  )

plot_grid(f2a,NULL, f2b, ncol = 3, labels = c("A", "B"), label_size = 18,  rel_widths = c(1, 0.2, 1.5))

ggsave(plot = plot_grid(f2a,NULL, f2b, ncol = 3, labels = c("A", "B"), label_size = 18,  rel_widths = c(1, 0.2, 1.5))
,
       filename="Figure2_2.png", width = 13, height = 6)
ggsave(plot = plot_grid(f2a,NULL, f2b, ncol = 3, labels = c("A", "B"), label_size = 18,  rel_widths = c(1, 0.1, 1.5))
,
       filename="Figure2.pdf", width = 13, height = 6)
