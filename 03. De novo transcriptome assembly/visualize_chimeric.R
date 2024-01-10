suppressPackageStartupMessages({
  
  library(gggenomes)
  library(ggplot2)
  library(reshape2)
  library(data.table)
  library(stringr)
  
})

args <- commandArgs(TRUE)


g0 <- read_feats(args[1])
transcript <- args[2]
transcript_flybase <- args[3]
out <- args[4]

g1<-g0[g0$type=="repeat_region",]

TE<-unique(g1$name)

#g0[g0$geom_id==transcript,]$name <- transcript

plotIns<-gggenomes(g0) %>% add_feats(TE=g1) + geom_gene(aes(fill=name),size = 5, position="strandpile") + geom_bin_label(size=6)+
  scale_fill_brewer("", palette="Dark2", na.value="cornsilk3") + theme_minimal(base_size = 16) + theme(panel.grid = element_blank(),axis.line.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(), axis.line.x=element_line(),axis.text.x = element_text(color = "black"),legend.position = "bottom") + labs(x="")
plotIns


ggsave(file=paste0(out,transcript,".TE.pdf"),plot = plotIns, width = 8, height = 4, bg="white")

