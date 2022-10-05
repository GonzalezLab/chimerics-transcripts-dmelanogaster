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

# Figure 6A

# All
ListStrains = list(`Head`=unique(TE_chim_global[TE_chim_global$tissue=="head",]$TEfamily),
                   `Gut`=unique(TE_chim_global[TE_chim_global$tissue=="gut",]$TEfamily),
                   `Ovary`=unique(TE_chim_global[TE_chim_global$tissue=="ovary",]$TEfamily))
ggvenn(ListStrains,
       stroke_size = 0.5, set_name_size = 4,fill_alpha=0.4,stroke_alpha = 1,stroke_color=NA ) + ggtitle("Shared TE families") + scale_fill_viridis_d() 

f4a<-plot(euler(ListStrains[c("Head","Gut","Ovary")]),
          quantities = list(type = c("counts"), cex = 1.2), labels = list(cex = 1.2),
          edges = list(col = NA),
          fills =list(fill = c("#CD96CD","#D9D9D9","#FFFACD"), alpha = 0.7) )

#createVD_3g(dataList = listTissues, filename = "~/Desktop/VD_by_tissue.png", mainTitle = "Shared chimeric gene-TE transcripts across body parts", ngroups = 3, category.names =names(listTissues))
# Group 1
ListStrains = list(`Head`=unique(TE_chim_global[TE_chim_global$tissue=="head" & TE_chim_global$group==1,]$TEfamily),
                   `Gut`=unique(TE_chim_global[TE_chim_global$tissue=="gut" & TE_chim_global$group==1,]$TEfamily),
                   `Ovary`=unique(TE_chim_global[TE_chim_global$tissue=="ovary" & TE_chim_global$group==1,]$TEfamily))
ggvenn(ListStrains,
       stroke_size = 0.5, set_name_size = 4,fill_alpha=0.4,stroke_alpha = 1,stroke_color=NA ) + ggtitle("Shared TE families") + scale_fill_viridis_d() 

f4a_g1<-plot(euler(ListStrains[c("Head","Gut","Ovary")]),
          quantities = list(type = c("counts"), cex = 1.2), labels = list(cex = 1.2),
          edges = list(col = NA),
          fills =list(fill = c("#CD96CD","#D9D9D9","#FFFACD"), alpha = 0.7) )
# Group 2
ListStrains = list(`Head`=unique(TE_chim_global[TE_chim_global$tissue=="head" & TE_chim_global$group==2,]$TEfamily),
                   `Gut`=unique(TE_chim_global[TE_chim_global$tissue=="gut" & TE_chim_global$group==2,]$TEfamily),
                   `Ovary`=unique(TE_chim_global[TE_chim_global$tissue=="ovary" & TE_chim_global$group==2,]$TEfamily))
ggvenn(ListStrains,
       stroke_size = 0.5, set_name_size = 4,fill_alpha=0.4,stroke_alpha = 1,stroke_color=NA ) + ggtitle("Shared TE families") + scale_fill_viridis_d() 

f4a_g2<-plot(euler(ListStrains[c("Head","Gut","Ovary")]),
             quantities = list(type = c("counts"), cex = 1.2), labels = list(cex = 1.2),
             edges = list(col = NA),
             fills =list(fill = c("#CD96CD","#D9D9D9","#FFFACD"), alpha = 0.7) )

# Figure 6B
#  pal_rickandmorty("schwifty", alpha = 1)(12)
color_map_families <- c(
  "roo"="#FAFD7CFF",
  "INE-1"="#82491EFF",
  "LARD"="#24325FFF",
  "P-element"="#B7E4F9FF",
  "1360"="#FB6467FF",
  "Gypsy-2_Dsim"="#526E2DFF",
  "S-element"="#E762D7FF",
  "jockey"="#E89242FF",
  "hopper2"="#FAE48BFF",
  "HB"="#A6EEE6FF",
  "297"= "#917C5DFF",
  "Others"="#69C8ECFF"
)

TE_chim_global_TE_commonTE <- unique(TE_chim_global_TE[,c("FlyBaseGeneRef", "TEfamily")])

occurrencesTE<-as.data.frame(table(TE_chim_global_TE_commonTE$TEfamily))
ggplot(occurrencesTE, aes(x=Var1,y=Freq)) + geom_col()
occurrencesTE<-as.data.frame(table(TE_chim_global_TE_commonTE$TEfamily))
ggplot(occurrencesTE, aes(x=Var1,y=Freq)) + geom_col()
occurrencesTE$Var1<-as.character(occurrencesTE$Var1)
occurrencesTE[occurrencesTE$Freq<15,]$Var1 <- "Others"
n<-sum(occurrencesTE$Var1 == "Others")
occurrencesTE$Var2 <- occurrencesTE$Var1
occurrencesTE[occurrencesTE$Var1=="Others",]$Var2 <- paste0("Others (",n,")")
occurrencesTE<-occurrencesTE %>% group_by(Var1,Var2) %>% summarise(Freq=sum(Freq))
occurrencesTE$Var1 <- factor(occurrencesTE$Var1, levels=c(occurrencesTE[occurrencesTE$Var1!="Others",] %>% arrange(.,Freq) %>% pull(Var1) %>%rev() %>% as.character(), "Others"))
occurrencesTE<-occurrencesTE %>% 
  group_by(Var1,Var2) %>% 
  mutate(per=Freq/sum(occurrencesTE$Freq)) %>% 
  arrange(desc(Freq))
occurrencesTE$label <- scales::percent(occurrencesTE$per)

occurrencesTE$label <- scales::percent(occurrencesTE$per, accuracy = 0.1)
f4b<-ggplot(occurrencesTE) + geom_col(width=1, aes(x="",y=per,fill=Var1)) + coord_polar("y", start=0) +
  theme_void() + scale_fill_manual(values=color_map_families) + geom_label_repel(aes(x="",y= cumsum(per) - per/2,label = label), size=4, show.legend = F, max.overlaps = 15, nudge_x = 0.1)  +  theme(plot.title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) + labs(fill="TE family")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16))
f4b



TE_chim_global_TE_commonTE <- unique(TE_chim_global_TE[,c("FlyBaseGeneRef", "TEfamily")])
occurrencesTE<-as.data.frame(table(TE_chim_global_TE_commonTE$TEfamily))
ggplot(occurrencesTE, aes(x=Var1,y=Freq)) + geom_col()
occurrencesTE<-as.data.frame(table(TE_chim_global_TE_commonTE$TEfamily))
ggplot(occurrencesTE, aes(x=Var1,y=Freq)) + geom_col()
occurrencesTE$Var1<-as.character(occurrencesTE$Var1)
occurrencesTE[occurrencesTE$Freq<15,]$Var1 <- "Others"
occurrencesTE<-aggregate(occurrencesTE$Freq, by=list(Category=occurrencesTE$Var1), FUN=sum)
occurrencesTE$Category <- factor(occurrencesTE$Category, levels=c("roo", "INE-1", "LARD", "P-element", "1360", "Gypsy-2_Dsim", "hopper2", "jockey", "S-element", "Others"))
occurrencesTE<-occurrencesTE %>% 
  group_by(Category) %>% 
  mutate(per=x/sum(occurrencesTE$x)) %>% 
  arrange(desc(Category))
occurrencesTE$label <- scales::percent(occurrencesTE$per)
occurrencesTE$label <- scales::percent(occurrencesTE$per, accuracy = 0.1)

p1<-ggplot(occurrencesTE) + geom_col(width=1, aes(x="",y=per,fill=Category)) + coord_polar("y", start=0) +
  theme_void() + scale_fill_manual(values=color_map_families) + geom_label_repel(aes(x="",y= cumsum(per) - per/2,label = label), size=4, show.legend = F, max.overlaps = 15, nudge_x = 0.1)  +  theme(plot.title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) + labs(fill="TE family")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16))
p1


TE_chim_global_TE_g1 <- TE_chim_global_TE[TE_chim_global_TE$group==1,]
TE_chim_global_TE_commonTE <- unique(TE_chim_global_TE_g1[,c("FlyBaseGeneRef", "TEfamily")])
occurrencesTE<-as.data.frame(table(TE_chim_global_TE_commonTE$TEfamily))
ggplot(occurrencesTE, aes(x=Var1,y=Freq)) + geom_col()
occurrencesTE<-as.data.frame(table(TE_chim_global_TE_commonTE$TEfamily))
ggplot(occurrencesTE, aes(x=Var1,y=Freq)) + geom_col()
occurrencesTE$Var1<-as.character(occurrencesTE$Var1)
occurrencesTE[occurrencesTE$Freq<15,]$Var1 <- "Others"
occurrencesTE<-aggregate(occurrencesTE$Freq, by=list(Category=occurrencesTE$Var1), FUN=sum)
occurrencesTE$Category <- factor(occurrencesTE$Category, levels=c( "INE-1", "roo",  "P-element", "LARD", "1360",   "S-element", "Gypsy-2_Dsim", "Others"))
occurrencesTE<-occurrencesTE %>% 
  group_by(Category) %>% 
  mutate(per=x/sum(occurrencesTE$x)) %>% 
  arrange(desc(Category))
occurrencesTE$label <- scales::percent(occurrencesTE$per)
occurrencesTE$label <- scales::percent(occurrencesTE$per, accuracy = 0.1)

p2<-ggplot(occurrencesTE) + geom_col(width=1, aes(x="",y=per,fill=Category)) + coord_polar("y", start=0) +
  theme_void() + scale_fill_manual(values=color_map_families) + geom_label_repel(aes(x="",y= cumsum(per) - per/2,label = label), size=4, show.legend = F, max.overlaps = 15, nudge_x = 0.1)  +  theme(plot.title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) + labs(fill="TE family")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16))
p2


TE_chim_global_TE_g2 <- TE_chim_global_TE[TE_chim_global_TE$group==2,]
TE_chim_global_TE_commonTE <- unique(TE_chim_global_TE_g2[,c("FlyBaseGeneRef", "TEfamily")])
occurrencesTE<-as.data.frame(table(TE_chim_global_TE_commonTE$TEfamily))
ggplot(occurrencesTE, aes(x=Var1,y=Freq)) + geom_col()
occurrencesTE<-as.data.frame(table(TE_chim_global_TE_commonTE$TEfamily))
ggplot(occurrencesTE, aes(x=Var1,y=Freq)) + geom_col()
occurrencesTE$Var1<-as.character(occurrencesTE$Var1)
occurrencesTE[occurrencesTE$Freq<15,]$Var1 <- "Others"
occurrencesTE<-aggregate(occurrencesTE$Freq, by=list(Category=occurrencesTE$Var1), FUN=sum)
occurrencesTE$Category <- factor(occurrencesTE$Category, levels=c(  "roo",  "INE-1","LARD", "hopper2", "Others"))
occurrencesTE<-occurrencesTE %>% 
  group_by(Category) %>% 
  mutate(per=x/sum(occurrencesTE$x)) %>% 
  arrange(desc(Category))
occurrencesTE$label <- scales::percent(occurrencesTE$per)
occurrencesTE$label <- scales::percent(occurrencesTE$per, accuracy = 0.1)

p3<-ggplot(occurrencesTE) + geom_col(width=1, aes(x="",y=per,fill=Category)) + coord_polar("y", start=0) +
  theme_void() + scale_fill_manual(values=color_map_families) + geom_label_repel(aes(x="",y= cumsum(per) - per/2,label = label), size=4, show.legend = F, max.overlaps = 15, nudge_x = 0.1)  +  theme(plot.title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) + labs(fill="TE family")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16))
p3

ggarrange(p1 + ggtitle("All"), p2+ ggtitle("Overlap & AS insertions"), p3+ ggtitle("Internal insertions"),   nrow=1, common.legend = TRUE, legend="bottom")
ggsave("Figure4_3.png", width = 12, height = 5)
ggsave("Figure4_3.pdf", width = 12, height = 5)

TE_chim_global_TE_120 <- TE_chim_global_TE[TE_chim_global_TE$lengthTE>=120,]
TE_chim_global_TE_commonTE <- unique(TE_chim_global_TE_120[,c("FlyBaseGeneRef", "TEfamily")])
occurrencesTE<-as.data.frame(table(TE_chim_global_TE_commonTE$TEfamily))
ggplot(occurrencesTE, aes(x=Var1,y=Freq)) + geom_col()
occurrencesTE<-as.data.frame(table(TE_chim_global_TE_commonTE$TEfamily))
ggplot(occurrencesTE, aes(x=Var1,y=Freq)) + geom_col()
occurrencesTE$Var1<-as.character(occurrencesTE$Var1)
occurrencesTE[occurrencesTE$Freq<15,]$Var1 <- "Others"
occurrencesTE<-aggregate(occurrencesTE$Freq, by=list(Category=occurrencesTE$Var1), FUN=sum)
occurrencesTE$Category <- factor(occurrencesTE$Category, levels=c(   "INE-1", "roo",  "LARD", "P-element", "1360", "Gypsy-2_Dsim", "S-element","jockey", "Others"))
occurrencesTE<-occurrencesTE %>% 
  group_by(Category) %>% 
  mutate(per=x/sum(occurrencesTE$x)) %>% 
  arrange(desc(Category))
occurrencesTE$label <- scales::percent(occurrencesTE$per)
occurrencesTE$label <- scales::percent(occurrencesTE$per, accuracy = 0.1)

p1_120<-ggplot(occurrencesTE) + geom_col(width=1, aes(x="",y=per,fill=Category)) + coord_polar("y", start=0) +
  theme_void() + scale_fill_manual(values=color_map_families) + geom_label_repel(aes(x="",y= cumsum(per) - per/2,label = label), size=4, show.legend = F, max.overlaps = 15, nudge_x = 0.1)  +  theme(plot.title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) + labs(fill="TE family")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16))
p1_120

TE_chim_global_TE_g1_120 <- TE_chim_global_TE[TE_chim_global_TE$group==1 & TE_chim_global_TE$lengthTE>=120,]
TE_chim_global_TE_commonTE <- unique(TE_chim_global_TE_g1_120[,c("FlyBaseGeneRef", "TEfamily")])
occurrencesTE<-as.data.frame(table(TE_chim_global_TE_commonTE$TEfamily))
ggplot(occurrencesTE, aes(x=Var1,y=Freq)) + geom_col()
occurrencesTE<-as.data.frame(table(TE_chim_global_TE_commonTE$TEfamily))
ggplot(occurrencesTE, aes(x=Var1,y=Freq)) + geom_col()
occurrencesTE$Var1<-as.character(occurrencesTE$Var1)
occurrencesTE[occurrencesTE$Freq<15,]$Var1 <- "Others"
occurrencesTE<-aggregate(occurrencesTE$Freq, by=list(Category=occurrencesTE$Var1), FUN=sum)
occurrencesTE$Category <- factor(occurrencesTE$Category, levels=c(   "INE-1", "roo",  "P-element",  "LARD", "1360", "Gypsy-2_Dsim", "S-element", "Others"))
occurrencesTE<-occurrencesTE %>% 
  group_by(Category) %>% 
  mutate(per=x/sum(occurrencesTE$x)) %>% 
  arrange(desc(Category))
occurrencesTE$label <- scales::percent(occurrencesTE$per)
occurrencesTE$label <- scales::percent(occurrencesTE$per, accuracy = 0.1)

p2_120<-ggplot(occurrencesTE) + geom_col(width=1, aes(x="",y=per,fill=Category)) + coord_polar("y", start=0) +
  theme_void() + scale_fill_manual(values=color_map_families) + geom_label_repel(aes(x="",y= cumsum(per) - per/2,label = label), size=4, show.legend = F, max.overlaps = 15, nudge_x = 0.1)  +  theme(plot.title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) + labs(fill="TE family")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16))
p2_120

TE_chim_global_TE_g2_120 <- TE_chim_global_TE[TE_chim_global_TE$group==2 & TE_chim_global_TE$lengthTE>=120,]
TE_chim_global_TE_commonTE <- unique(TE_chim_global_TE_g2_120[,c("FlyBaseGeneRef", "TEfamily")])
occurrencesTE<-as.data.frame(table(TE_chim_global_TE_commonTE$TEfamily))
ggplot(occurrencesTE, aes(x=Var1,y=Freq)) + geom_col()
occurrencesTE<-as.data.frame(table(TE_chim_global_TE_commonTE$TEfamily))
ggplot(occurrencesTE, aes(x=Var1,y=Freq)) + geom_col()
occurrencesTE$Var1<-as.character(occurrencesTE$Var1)
occurrencesTE[occurrencesTE$Freq<15,]$Var1 <- "Others"
occurrencesTE<-aggregate(occurrencesTE$Freq, by=list(Category=occurrencesTE$Var1), FUN=sum)
occurrencesTE$Category <- factor(occurrencesTE$Category, levels=c(   "roo",  "INE-1",  "LARD", "Others"))
occurrencesTE<-occurrencesTE %>% 
  group_by(Category) %>% 
  mutate(per=x/sum(occurrencesTE$x)) %>% 
  arrange(desc(Category))
occurrencesTE$label <- scales::percent(occurrencesTE$per)
occurrencesTE$label <- scales::percent(occurrencesTE$per, accuracy = 0.1)

p3_120<-ggplot(occurrencesTE) + geom_col(width=1, aes(x="",y=per,fill=Category)) + coord_polar("y", start=0) +
  theme_void() + scale_fill_manual(values=color_map_families) + geom_label_repel(aes(x="",y= cumsum(per) - per/2,label = label), size=4, show.legend = F, max.overlaps = 15, nudge_x = 0.1)  +  theme(plot.title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) + labs(fill="TE family")  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16))
p3_120

ggarrange(p1_120 + ggtitle("All"), p2_120+ ggtitle("Overlap & AS insertions"), p3_120+ ggtitle("Internal insertions"),   nrow=1, common.legend = TRUE, legend="bottom")
ggsave("FigureS_120.png", width = 12, height = 5)
ggsave("FigureS_120.pdf", width = 12, height = 5)


# Figure 6C
dataplotTE[dataplotTE$exonPosDescription=="Internal exon",]$exonPosDescription <- "Internal\nexon"
dataplotTE$exonPosDescription<-factor(dataplotTE$exonPosDescription,level=c("5'UTR", "Internal\nexon","3'UTR"))
f4c<-ggplot(dataplotTE , aes(x=exonPosDescription, fill=Order)) + geom_bar(position="fill") + facet_grid(.~tissue) + theme_Publication() + scale_colour_Publication()  + labs(x = "", y = "Proportion of gene-TE\nchimeric genes")+scale_fill_manual(values=c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948")) +#+scale_fill_tableau() 
  scale_y_continuous(labels=scales::percent) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text(size=16), plot.title = element_text(size = 16))

titleA <- ggdraw() + draw_label("Shared TE families across body parts", fontface='bold')
titleB <- ggdraw() + draw_label("TE superfamily distribution\nacross chimeric gene-TE transcripts", fontface='bold')
titles<-plot_grid(titleA, titleB, align = "h", axis = "t", labels = c('A', 'B'), label_size = 18)

plots <- align_plots(f4a, f4c, align = 'v', axis = 'l')
# then build the bottom row
top_row <- plot_grid(f4a, f4b, align = "h", axis = "t", rel_widths = c(1,1.5))
top_row_title <- plot_grid(titles, top_row,  rel_heights = c(0.1, 1), nrow = 2, ncol=1)

# then combine with the top row for final plot
plot_grid(top_row_title, NULL, f4c, labels = c('', '', 'C'), label_size = 18, ncol = 1, rel_heights = c(1,0.1, 1.3), nrow=3)
   
ggsave(plot_grid(top_row_title, NULL, f4c, labels = c('', '', 'C'), label_size = 18, ncol = 1, rel_heights = c(1,0.1, 1.3), nrow=3) ,file="Figure4.png", width = 9, height = 9)
ggsave(plot_grid(top_row_title, NULL, f4c, labels = c('', '', 'C'), label_size = 18, ncol = 1, rel_heights = c(1,0.1, 1.3), nrow=3) ,file="Figure4.pdf", width = 9, height = 9)

ggsave( plot_grid(titles, top_row,  rel_heights = c(0.1, 1), nrow = 2, ncol=1),  file="Figure4_2.pdf", width = 8.5, height = 4.5)
ggsave( plot_grid(titles, top_row,  rel_heights = c(0.1, 1), nrow = 2, ncol=1),  file="Figure4_2.png", width = 8.5, height = 4.5)

