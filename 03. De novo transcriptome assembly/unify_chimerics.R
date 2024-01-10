library(data.table)
library(dplyr)

chimerics_tmp <- fread("clean_set_minimap2/chimerics.tmp.tab", header = F)
colnames(chimerics_tmp) <- c("tissue", "strain","transcript_trinity", "transcript_stringtie", "class","transcript", "gene", "lengthTranscript", "totalExons", "posExon", "exon", "coordExon", "coordTE", "type", "TE_consensus","TE_family", "TE_superfamily","TE_order","TE_class","TE_length_incorporated","TE_length_total", "CP", "geneStatus", "expr", "SS", "result", "group", "roo_type")

chimerics<-chimerics_tmp %>% group_by(tissue, strain, transcript_trinity, transcript_stringtie, coordTE) %>% mutate(TE_length_incorporated_ok=sum(TE_length_incorporated))

chimerics$TE_length_incorporated <- chimerics$TE_length_incorporated_ok
chimerics$TE_length_incorporated_ok <- NULL

write.table(chimerics, "clean_set_minimap2/chimerics.tab", col.names = T, row.names = F, quote = F, sep = "\t")
