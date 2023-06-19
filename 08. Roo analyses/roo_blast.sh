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

> resultsBLAST.strict.txt
while read insertionInfo
do
tissue=$(echo "$insertionInfo" | cut -f 1)
strain=$(echo "$insertionInfo" | cut -f 2)
transcript_trinity=$(echo "$insertionInfo" | cut -f 4)
transcript_stringtie=$(echo "$insertionInfo" | cut -f 3)
transcript=$(echo "$insertionInfo" | cut -f 5)
gene=$(echo "$insertionInfo" | cut -f 6)

startRoo=$(echo "$insertionInfo" | cut -f 7)
startRoo=$(( startRoo - 1  ))
endRoo=$(echo "$insertionInfo" | cut -f 8)
cons=$(echo "$insertionInfo" | cut -f 9)
echo -e "$tissue\t$strain\t$transcript\t$transcript_trinity"
i=$(grep -Pl "^(?=.*$startRoo)(?=.*$endRoo)"  ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.roo.*.fa | grep -o "roo\.[0-9]*" | sed "s/roo\.//" | head -n1 )
file=$(ls ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.roo.$i.fa)

cp $file tmp/$transcript.$transcript_trinity.roo.$i.fa
sed -i "1 s/^.*$/>${transcript_trinity}:${tissue}:${strain}:${startRoo}:${endRoo}/" tmp/$transcript.$transcript_trinity.roo.$i.fa

blastn -query tmp/$transcript.$transcript_trinity.roo.$i.fa -db dmel-chr.fasta -dust no -soft_masking false -outfmt 6 -word_size 7 -evalue 0.05 -gapopen 5 -gapextend 2  -qcov_hsp_perc 85 -perc_identity 75 | head -n 20 >> resultsBLAST.strict.txt

echo -e "$tissue\t$assembly\t$transcript"

done <  <(awk -F'\t' ' $28 == "repeat" ' ${DIR}/clean_set_minimap2/chimerics_v2.tab | cut -f 1,2,3,4,6,7,13,15 |tr ':' '\t' |sort -u)

# I recovered some other to repeat (recover_roo_other), so I add them:
comm -13 <(cat ../roo_analysis/resultsBLAST.strict.txt | cut -f1 -d':' | sort -u)  <(awk -F'\t' ' $30 == "repeat" ' ${DIR}/clean_set_minimap2/chimerics_v2.tab | cut -f 1,2,3,4,6,7,13,15 |tr ':' '\t' |sort -u  | cut -f 4 | sort -u) > list_missing

while read insertionInfo
do
tissue=$(echo "$insertionInfo" | cut -f 1)
strain=$(echo "$insertionInfo" | cut -f 2)
transcript_trinity=$(echo "$insertionInfo" | cut -f 4)
transcript_stringtie=$(echo "$insertionInfo" | cut -f 7)
transcript=$(echo "$insertionInfo" | cut -f 3)
gene=$(echo "$insertionInfo" | cut -f 8)

startRoo=$(echo "$insertionInfo" | cut -f 5)
startRoo=$(( startRoo - 1  ))
endRoo=$(echo "$insertionInfo" | cut -f 6)
cons=$(echo "$insertionInfo" | cut -f 9)
echo -e "$tissue\t$strain\t$transcript\t$transcript_trinity"
i=$(grep -Pl "^(?=.*$startRoo)(?=.*$endRoo)"  ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.roo.*.fa | grep -o "roo\.[0-9]*" | sed "s/roo\.//" | head -n1 )
file=$(ls ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.roo.$i.fa)

cp $file tmp/$transcript.$transcript_trinity.roo.$i.fa
sed -i "1 s/^.*$/>${transcript_trinity}:${tissue}:${strain}:${startRoo}:${endRoo}/" tmp/$transcript.$transcript_trinity.roo.$i.fa

blastn -query tmp/$transcript.$transcript_trinity.roo.$i.fa -db dmel-chr.fasta -dust no -soft_masking false -outfmt 6 -word_size 7 -evalue 0.05 -gapopen 5 -gapextend 2  -qcov_hsp_perc 85 -perc_identity 75 | head -n 20 >> resultsBLAST.strict.txt

echo -e "$tissue\t$assembly\t$transcript"

done <  <(awk -F'\t' ' $30 == "repeat" ' ${DIR}/clean_set_minimap2/chimerics_v2.tab | grep -f list_missing | cut -f 1,2,3,4,6,7,9,16 |tr ':' '\t' |sort -u)




#awk -v OFS='\t' '$3 == "gene"' ../../../trinity-repeatMasker/resultsChimAnalysis/roo_analysis/dmel-chr-r6.31.gff | cut -f1-9 > ../../../trinity-repeatMasker/resultsChimAnalysis/roo_analysis/dmel-chr-genes-r6.31.gff 
#cat Natural\ TE-* > Natural_TEs.gff3
#awk -v OFS='\t' '{ print $2,"blast","region",$9,$10,".",".",".",$1 }'o resultsBLAST.strict.txt | cut -f1 -d':' > resultsBLAST.gff

awk -v OFS='\t' '{
if ( $13 > $14 )
print $6,"blast","region",$14,$13,".",".",".",$1";"$2";"$3";"$4+1":"$5";evalue="$15";bit-score="$16
else
 print $6, "blast","region",$13,$14,".",".",".",$1";"$2";"$3";"$4+1":"$5";evalue="$15";bit-score="$16
}' <(cat resultsBLAST.strict.txt | tr ':' '\t' )  > resultsBLAST.strict.gff


##muy lento
# > resultsBLAST_top20.txt 
# while read transcript
# do
# awk -v transcript="$transcript" ' $1 == transcript ' resultsBLAST.strict.txt | head -n 20 >> resultsBLAST_top20.txt 
# done < <(cut -f1 resultsBLAST.strict.txt | sort -u )


# awk -v OFS='\t' '{
# if ( $9 > $10 )
# print $2,"blast","region",$10,$9,".",".",".",$1";evalue="$11";bit-score="$12
# else
#  print $2,"blast","region",$9,$10,".",".",".",$1";evalue="$11";bit-score="$12
# }' resultsBLAST_top20.txt  | cut -f1 -d':' > resultsBLAST_top20.gff
##
 

#cat resultsBLAST.gff | sort -u > resultsBLAST.strict.unique2.gff
#opcion ordenada
cat resultsBLAST.strict.gff | awk '!x[$0]++' > resultsBLAST.strict.unique2.gff


#awk -v OFS='\t' '{
#if ( $9 > $10 )
#print $2,$10-1,$9 
#else
# print $2,$9-1,$10
#}' resultsBLAST.strict.txt  > resultsBLAST.bed

#cat dmel-chr-genes-r6.31.gff Natural_TEs.gff3 > dmel-chr-genes-TE-r6.31.tmp.gff 
#gff3sort.pl --precise dmel-chr-genes-TE-r6.31.tmp.gff > dmel-chr-genes-TE-r6.31.gff

#bedtools intersect  -a resultsBLAST.strict.unique2.gff -b dmel-chr-genes-TE-r6.31.gff  -wa -wb -v | cut -f 1,4,5 | sort -u | wc -l
#bedtools intersect  -a resultsBLAST.strict.unique2.gff -b dmel-chr-genes-TE-r6.31.gff -wb > resultsBLAST.strict.unique.bedtoolsGenomeTes.tab
#bedtools intersect  -a resultsBLAST.strict.unique2.gff -b dmel-chr-genes-TE-r6.31.gff -wb -v > resultsBLAST.strict.unique.bedtoolsInterGenic.tab

#cat resultsBLAST.strict.unique.bedtoolsGenomeTes.tab resultsBLAST.strict.unique.bedtoolsInterGenic.tab > resultsBLAST.strict.unique.bedtoolsAll.tab

#grep "roo" resultsBLAST.strict.unique.bedtoolsAll.tab | grep "transposable" | cut -f 9| cut -f1 -d';' | sort -u | sed "s/_head_/\thead\t/g" | sed "s/_gut_/\tgut\t/g" | sed "s/_ovary_/\tovary\t/g" > resultsBLAST.strict.unique.bedtoolsAll.roo.lst

# > resultsBLAST.strict.unique2_top20.gff 
# while read transcript
# do
# grep "$transcript" resultsBLAST.strict.unique2.gff | head -n 20 >> resultsBLAST.strict.unique2_top20.gff 
# done < <(cut -f9 resultsBLAST.strict.unique2.gff | cut -f1 -d ';' | sort -u)

bedtools intersect  -a resultsBLAST.strict.unique2.gff -b dmel-chr-genes-TE-r6.31.gff -wb > resultsBLAST.strict.unique2_top20.bedtoolsGenomeTes.tab
bedtools intersect  -a resultsBLAST.strict.unique2.gff -b dmel-chr-genes-TE-r6.31.gff -wb -v > resultsBLAST.strict.unique2_top20.bedtoolsInterGenic.tab

cat resultsBLAST.strict.unique2_top20.bedtoolsGenomeTes.tab resultsBLAST.strict.unique2_top20.bedtoolsInterGenic.tab > resultsBLAST.strict.unique2_top20.bedtoolsAll.tab

>  resultsBLAST.strict.unique2_top20.cross.gff
while read transcriptInfo
do
start=$(echo "$transcriptInfo" | cut -f4)
end=$(echo "$transcriptInfo" | cut -f5)
transcript=$(echo "$transcriptInfo"  | cut -f9 | cut -f1 -d ';' )
result=$(grep $transcript resultsBLAST.strict.unique2_top20.bedtoolsAll.tab | grep $start | grep $end |  cut -f 18 | sed -r '/^\s*$/d' | cut -f1 -d';' | tr '\n' ';' )
if [ -z "$result" ];then
echo -e "$transcriptInfo\tintergenic" >>  resultsBLAST.strict.unique2_top20.cross.gff
else
echo -e "$transcriptInfo\t$result"  >>  resultsBLAST.strict.unique2_top20.cross.gff
fi
done <  resultsBLAST.strict.unique2.gff 

> resultsBLAST.strict.unique2_top20.cross.result2.gff
while read transcriptInfo
do
transcript=$(echo "$transcriptInfo" | cut -f1 -d ';' )
tissue=$(echo "$transcriptInfo"  | cut -f2 -d ';' )
strain=$(echo "$transcriptInfo" | cut -f3 -d ';')
pos=$(echo "$transcriptInfo" | cut -f4 -d ';')
nMatches=$(grep "$transcript" resultsBLAST.strict.unique2_top20.cross.gff | grep $tissue | grep $strain | grep $pos | wc -l)
firstGene=$(grep "$transcript" resultsBLAST.strict.unique2_top20.cross.gff | grep $tissue | grep $strain | grep $pos | head -n1 | cut -f 10| sed 's/ID=//g' | tr -d ';' | sort -u)
nroo=$(grep "$transcript" resultsBLAST.strict.unique2_top20.cross.gff | grep $tissue | grep $strain | grep $pos | grep -n -m 1 "roo" | cut -f1 -d':')
ngene=$(grep "$transcript" resultsBLAST.strict.unique2_top20.cross.gff | grep $tissue | grep $strain | grep $pos  | grep -n -m 1 "FBgn" | cut -f1 -d':')
nintergenic=$(grep "$transcript" resultsBLAST.strict.unique2_top20.cross.gff | grep $tissue | grep $strain | grep $pos  | grep -n -m 1 "intergenic" | cut -f1 -d':')
gene=$(grep -w $transcript ${DIR}/clean_set_minimap2/chimerics_v2.tab | grep -w $tissue | grep -w $strain | cut -f9 | sort -u)
if [[ $firstGene == $gene ]];then
geneCheck="YES"
else
geneCheck="NO"
fi

if [[ -z "$nroo" ]];then
nroo="NA"
fi

if [[ -z "$nintergenic" ]];then
nintergenic="NA"
fi

echo -e "$transcriptInfo"
echo -e "$transcript\t$tissue\t$strain\tpos\t$pos\tmatches\t$nMatches\troo\t$nroo\tgene\t$ngene\tintergenic\t$nintergenic\t$geneCheck" >> resultsBLAST.strict.unique2_top20.cross.result2.gff

done <   <(cut -f9 resultsBLAST.strict.unique2_top20.cross.gff | cut -f1-4 -d';' | sort -u )  





