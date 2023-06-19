#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name=roo
#SBATCH --output=logs/roo.out
#SBATCH --error=logs/roo.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate chimerics
module load foss/2021b
module load R/4.1.2

DIR="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/chimerics"
DIRDATA="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/DATA"

#tissues="head=Heads gut=Gut ovary=Ovary"
tissues="head gut ovary"
assembly="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"

#wget http://ftp.flybase.net/releases/FB2019_06/dmel_r6.31/fasta/dmel-all-gene-r6.31.fasta.gz

#gunzip dmel-all-gene-r6.31.fasta.gz

#cat roo.fasta dmel-all-gene-r6.31.fasta > dmel-genes-roo.fasta

#makeblastdb -in  dmel-genes-roo.fasta -dbtype nucl
#makeblastdb -in  dmel-chr.fasta -dbtype nucl

cd ${DIR}/roo_analysis

while read insertionInfo
do
tissue=$(echo "$insertionInfo" | cut -f 1)
strain=$(echo "$insertionInfo" | cut -f 2)
transcript_stringtie=$(echo "$insertionInfo" | cut -f 3)
transcript_trinity=$(echo "$insertionInfo" | cut -f 4)
transcript=$(echo "$insertionInfo" | cut -f 5)
gene=$(echo "$insertionInfo" | cut -f 6)

startRoo=$(echo "$insertionInfo" | cut -f 7)
startRoo=$(( startRoo - 1  ))
endRoo=$(echo "$insertionInfo" | cut -f 8)
cons=$(echo "$insertionInfo" | cut -f 9)
#echo -e "$tissue\t$strain\t$transcript\t$transcript_trinity"
i=$(grep -Pl "^(?=.*$startRoo)(?=.*$endRoo)"  ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.roo.*.fa | grep -o "roo\.[0-9]*" | sed "s/roo\.//" | head -n1 )
start=NA
end=NA
roo=$(blastn -query ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.roo.$i.fa -db ${DIR}/roo/roo.fasta -dust no -soft_masking false -word_size 7 -outfmt 6 -max_target_seqs 1 -evalue 0.05 -gapopen 5 -gapextend 2 | head -n1 | cut -f 9,10 )

start=$(echo "$roo" | cut -f1)
end=$(echo "$roo" | cut -f2)

if [ $start -gt $end ]; then
echo -e "$tissue\t$strain\t$transcript\t$transcript_trinity\t$roo_type\t$end\t$start" 
else
echo -e "$tissue\t$strain\t$transcript\t$transcript_trinity\t$roo_type\t$start\t$end" 
fi

done <  <(awk -F'\t' ' $28 == "repeat" ' ${DIR}/clean_set_minimap2/chimerics_v2.tab | cut -f 1,2,3,4,6,7,13,15 |tr ':' '\t' |sort -u) > roo_match_position.tab