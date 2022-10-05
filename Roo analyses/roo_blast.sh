#!/bin/bash

# set the partition where the job will run

#SBATCH --partition=normal

# set the number of nodes

#SBATCH --cpus-per-task=12

#SBATCH --mem=78G

# mail alert at start, end and abortion of execution

#SBATCH --mail-type=ALL

#SBATCH --job-name=roo

#SBATCH --output=roo.output

#SBATCH --error=roo.error

# send mail to this address

#SBATCH --mail-user=marta.coronado@ibe.upf-csic.es

module load Miniconda3/4.7.10

source activate /homes/users/mcoronado/.conda/envs/trinity

# Run

DIR=/homes/users/mcoronado/scratch/5GenomesProject/ANALYSES

#tissues="head=Heads gut=Gut ovary=Ovary"
tissues="head gut ovary"
assembly="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"

#wget http://ftp.flybase.net/releases/FB2019_06/dmel_r6.31/fasta/dmel-all-gene-r6.31.fasta.gz

#gunzip dmel-all-gene-r6.31.fasta.gz

cat roo.fasta dmel-all-gene-r6.31.fasta > dmel-genes-roo.fasta

makeblastdb -in  dmel-genes-roo.fasta -dbtype nucl
makeblastdb -in  dmel-chr.fasta -dbtype nucl

> resultsBLAST.txt


while read transcriptInfo
do
tissue=$(echo "$transcriptInfo" | cut -f1 )
assembly=$(echo "$transcriptInfo" | cut -f2 )
transcript=$(echo "$transcriptInfo" | cut -f3 )
TEconsensusID=$(echo "$transcriptInfo" | cut -f6 )
TEfamily=$(echo "$transcriptInfo" | cut -f7 )
exonPosDescription=$(echo "$transcriptInfo" | cut -f8 )
length=$(echo "$transcriptInfo" | cut -f9 )
group=$(echo "$transcriptInfo" | cut -f10 )

awk -v transcript="$transcript" '$1 == transcript ' ${DIR}/trinity-repeatMasker/Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/Trinity.hits_bin90.fasta.out.gff | grep con48_roo | bedtools merge > ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/tmp/${transcript}_${tissue}_${assembly}_con48_roo.gff

bedtools getfasta \
-fi ${DIR}/trinity-repeatMasker/Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/Trinity.hits_bin90.fasta.masked \
-bed ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/tmp/${transcript}_${tissue}_${assembly}_con48_roo.gff > \
${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/tmp/${transcript}_${tissue}_${assembly}_con48_roo.fasta

sed -i "s/^>.*$/>${transcript}_${tissue}_${assembly}/g" ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/tmp/${transcript}_${tissue}_${assembly}_con48_roo.fasta

blastn -query ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/tmp/${transcript}_${tissue}_${assembly}_con48_roo.fasta -db dmel-chr.fasta -dust no -soft_masking false -outfmt 6 -word_size 7 -evalue 0.05 -gapopen 5 -gapextend 2 >> resultsBLAST.txt

echo -e "$tissue\t$assembly\t$transcript"

done <  ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/transcripts_TEfamily_group.tab

#awk -v OFS='\t' '$3 == "gene"' ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/dmel-chr-r6.31.gff | cut -f1-9 > ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/dmel-chr-genes-r6.31.gff 
#cat Natural\ TE-* > Natural_TEs.gff3
#awk -v OFS='\t' '{ print $2,"blast","region",$9,$10,".",".",".",$1 }'o resultsBLAST.txt | cut -f1 -d':' > resultsBLAST.gff

awk -v OFS='\t' '{
if ( $9 > $10 )
print $2,"blast","region",$10,$9,".",".",".",$1";evalue="$11";bit-score="$12
else
 print $2,"blast","region",$9,$10,".",".",".",$1";evalue="$11";bit-score="$12
}' resultsBLAST.txt | cut -f1 -d':' > resultsBLAST.gff


##muy lento
> resultsBLAST_top20.txt 
while read transcript
do
awk -v transcript="$transcript" ' $1 == transcript ' resultsBLAST.txt | head -n 20 >> resultsBLAST_top20.txt 
done < <(cut -f1 resultsBLAST.txt | sort -u )


awk -v OFS='\t' '{
if ( $9 > $10 )
print $2,"blast","region",$10,$9,".",".",".",$1";evalue="$11";bit-score="$12
else
 print $2,"blast","region",$9,$10,".",".",".",$1";evalue="$11";bit-score="$12
}' resultsBLAST_top20.txt  | cut -f1 -d':' > resultsBLAST_top20.gff
##
 

cat resultsBLAST.gff | sort -u > resultsBLAST.unique.gff
#opcion ordenada
cat resultsBLAST.gff | awk '!x[$0]++' > resultsBLAST.unique2.gff


#awk -v OFS='\t' '{
#if ( $9 > $10 )
#print $2,$10-1,$9 
#else
# print $2,$9-1,$10
#}' resultsBLAST.txt  > resultsBLAST.bed

#cat dmel-chr-genes-r6.31.gff Natural_TEs.gff3 > dmel-chr-genes-TE-r6.31.tmp.gff 
#gff3sort.pl --precise dmel-chr-genes-TE-r6.31.tmp.gff > dmel-chr-genes-TE-r6.31.gff

#bedtools intersect  -a resultsBLAST.unique.gff -b dmel-chr-genes-TE-r6.31.gff  -wa -wb -v | cut -f 1,4,5 | sort -u | wc -l
bedtools intersect  -a resultsBLAST.unique.gff -b dmel-chr-genes-TE-r6.31.gff -wb > resultsBLAST.unique.bedtoolsGenomeTes.tab
bedtools intersect  -a resultsBLAST.unique.gff -b dmel-chr-genes-TE-r6.31.gff -wb -v > resultsBLAST.unique.bedtoolsInterGenic.tab

cat resultsBLAST.unique.bedtoolsGenomeTes.tab resultsBLAST.unique.bedtoolsInterGenic.tab > resultsBLAST.unique.bedtoolsAll.tab

> resultsBLAST.unique2_top20.gff 
while read transcript
do
grep "$transcript" resultsBLAST.unique2.gff | head -n 20 >> resultsBLAST.unique2_top20.gff 
done < <(cut -f9 resultsBLAST.unique2.gff | cut -f1 -d ';' | sort -u)

bedtools intersect  -a resultsBLAST.unique2_top20.gff -b dmel-chr-genes-TE-r6.31.gff -wb > resultsBLAST.unique2_top20.bedtoolsGenomeTes.tab
bedtools intersect  -a resultsBLAST.unique2_top20.gff -b dmel-chr-genes-TE-r6.31.gff -wb -v > resultsBLAST.unique2_top20.bedtoolsInterGenic.tab

cat resultsBLAST.unique2_top20.bedtoolsGenomeTes.tab resultsBLAST.unique2_top20.bedtoolsInterGenic.tab > resultsBLAST.unique2_top20.bedtoolsAll.tab

while read transcriptInfo
do
start=$(echo "$transcriptInf"o | cut -f4)
end=$(echo "$transcriptInfo" | cut -f5)
transcript=$(echo "$transcriptInfo"  | cut -f9 | cut -f1 -d ';' )
result=$(grep $transcript resultsBLAST.unique2_top20.bedtoolsAll.tab | grep $start | grep $end |  cut -f 18 | sed -r '/^\s*$/d' | cut -f1 -d';' | tr '\n' ';' )
if [ -z "$result" ];then
echo -e "$transcriptInfo\tintergenic"
else
echo -e "$transcriptInfo\t$result"
fi
done <  resultsBLAST.unique2_top20.gff | cut -f 1,9,10 >  resultsBLAST.unique2_top20.cross.gff

> resultsBLAST.unique2_top20.cross.result.gff
while read transcriptInfo
do
firstGene=$(grep -m1 "$transcriptInfo" resultsBLAST.unique2_top20.cross.gff | cut -f 3| sed 's/ID=//g' | tr -d ';')
assembly=$(echo "$transcriptInfo" | rev | cut -f1 -d'_' | rev)
tissue=$(echo "$transcriptInfo" | rev | cut -f2 -d'_' | rev)
transcript=$(echo "$transcriptInfo" | rev | cut -f3- -d'_' | rev)
nroo=$(grep "$transcriptInfo" resultsBLAST.unique2_top20.cross.gff | grep -n -m 1 "roo" | cut -f1 -d':')
ngene=$(grep "$transcriptInfo" resultsBLAST.unique2_top20.cross.gff | grep -n -m 1 "FBgn" | cut -f1 -d':')
nintergenic=$(grep "$transcriptInfo" resultsBLAST.unique2_top20.cross.gff | grep -n -m 1 "intergenic" | cut -f1 -d':')
gene=$(grep -w $transcript ../TE_chimeric_global_REVISED_v2.tab | grep -w $tissue | grep -w $assembly | cut -f8 | sort -u)
if [ $firstGene == $gene ];then
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
echo -e "$transcript\t$tissue\t$assembly\troo\t$nroo\tgene\t$ngene\tintergenic\t$nintergenic\t$geneCheck">> resultsBLAST.unique2_top20.cross.result.gff

done <   <(cut -f2 resultsBLAST.unique2_top20.cross.gff | cut -f1 -d';' | sort -u )  

