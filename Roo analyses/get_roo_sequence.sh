#!/bin/bash

# set the partition where the job will run

#SBATCH --partition=normal

# set the number of nodes

#SBATCH --cpus-per-task=2

#SBATCH --mem=16G

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


# keep transcripts with a roo
mkdir ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/tmp/

cat ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/TE_chimeric_global_REVISED_v2.tab | cut -f1,2,4,6,8,13,14,17,21,24 | awk '$7 == "roo" '  |  sort -u > ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/transcripts_TEfamily_group.tab

makeblastdb -in roo.fasta -dbtype nucl

fimo --oc fimo_repetitive_region --verbosity 1 --thresh 1.0E-4 ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/20220928102932_JASPAR2022_combined_matrices_21733_meme.txt ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/repetitive_region_roo.fasta

fimo --oc fimo_repetitive_region --verbosity 1 --thresh 1.0E-4 ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/20220928102932_JASPAR2022_combined_matrices_21733_meme.txt ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/repetitive_region_roo.fasta

fimo --oc fimo_whole_roo --verbosity 1 --thresh 1.0E-4 ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/20220928102932_JASPAR2022_combined_matrices_21733_meme.txt ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/roo.fasta


#> ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/results.tab
#> "${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/roo_group_1.fasta"
#> "${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/roo_group_2.fasta"
#> "${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/roo_group_1_Middle exon.fasta"
#> "${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/roo_group_1_First exon.fasta"
#> "${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/roo_group_1_Last exon.fasta"
#> "${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/roo_group_2_Middle exon (TE inside).fasta"
#> "${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/roo_group_2_First exon (TE inside).fasta"
#> "${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/roo_group_2_Last exon (TE inside).fasta"
> results.fimo.tab

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

fimo --oc fimo --verbosity 1 --thresh 1.0E-4 ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/20220928102932_JASPAR2022_combined_matrices_21733_meme.txt ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/tmp/${transcript}_${tissue}_${assembly}_con48_roo.fasta


while read motifD
do
motif=$(echo "$motifD" | cut -f1,7,9,10)
echo -e "$tissue\t$assembly\t$transcript\t$motif" >> results.fimo.tab
#echo -e "$tissue\t$assembly\t$transcript\t$length\t$group\t$motif" >> results.fimo.tab
done < <(tail -n+2 fimo/fimo.tsv | grep -v '#' | grep -v "^$") 

cat ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/tmp/${transcript}_${tissue}_${assembly}_con48_roo.fasta #>> "${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/roo_group_${group}_${exonPosDescription}.fasta"

cat ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/tmp/${transcript}_${tissue}_${assembly}_con48_roo.fasta #>> "${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/roo_group_${group}.fasta"

#start=$(blastn -query ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/tmp/${transcript}_${tissue}_${assembly}_con48_roo.fasta -db roo.fasta -dust no -soft_masking false -word_size 7 -outfmt 6 -max_target_seqs 1 -evalue 0.05 -gapopen 5 -gapextend 2 | head -n1 | cut -f 9,10 | tr '\t' '\n' | sort | head -n1)

#end=$(blastn -query ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/tmp/${transcript}_${tissue}_${assembly}_con48_roo.fasta -db roo.fasta -dust no -soft_masking false -word_size 7 -outfmt 6 -max_target_seqs 1 -evalue 0.05 -gapopen 5 -gapextend 2 | head -n1 | cut -f 9,10 | tr '\t' '\n' | sort | tail -n1)

echo -e "$tissue\t$assembly\t$transcript\t$length\t$group\t$start\t$end"

done < ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/transcripts_TEfamily_group.tab #>> ${DIR}/trinity-repeatMasker/resultsChimAnalysis/roo_analysis/results.tab


cat results.fimo.tab | sort -u > results.fimo.sort.unique.tab


