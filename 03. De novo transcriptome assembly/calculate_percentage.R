suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(rstatix)
})

args <- commandArgs(TRUE)
file <- args[1]
tissue <- args[2]
strain <- args[3]

blastn_result <- read.table(file)

colnames(blastn_result) <- c("qseqid","qlen","sseqid","slen", "qcovs","pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

blastn_result <- blastn_result[blastn_result$length > 100,]
blastn_result$scov <- ifelse(blastn_result$length > blastn_result$slen, blastn_result$slen/blastn_result$slen, blastn_result$length/blastn_result$slen)
blastn_result <- blastn_result[blastn_result$scov>0.8,]

reference_transcripts <- read.table("reference_transcripts.lst", header = T)
blastn_result_merge <- merge(blastn_result, reference_transcripts[,c("ref_gene_id", "qry_id")], by.x="sseqid", by.y="qry_id")

blastn_result_merge_filter <- blastn_result_merge %>% group_by(qseqid,ref_gene_id) %>%  filter(bitscore == max(bitscore))

blastn_result_merge_filter<-blastn_result_merge_filter[,c("qseqid","qlen","sseqid","ref_gene_id","slen", "qcovs","pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore","scov")]

write.table(blastn_result_merge_filter,file=paste0("check_fussion/output/",strain,"/",tissue,"/chimerics.filterLengthMin100.filterScovMin80.tmp.blastn"), quote = F, sep = "\t", col.names = F, row.names = F)
