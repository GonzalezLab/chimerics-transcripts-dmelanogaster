#!/bin/bash

#SBATCH -J fastp
#SBATCH --cpus-per-task=12
#SBATCH --mem 16G 
#SBATCH -o logs/fastp.out # Standard output 
#SBATCH -e logs/fastp.err # Standard error 

source ~/miniconda3/etc/profile.d/conda.sh
conda activate chimerics

DIR="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/chimerics"
DIRDATA="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/DATA"

cd $DIR

# tissues="gut head ovary"

# for tissue in $tissues
# 	do
# 		while read sample
# 		do
# 			case "$sample" in
# 				S*)
# 					echo -e "$sample\tSLA-001\t$tissue"
# 				;;
# 				T*)
# 					echo -e "$sample\tTOM-007\t$tissue"
# 				;;
# 				A*)
# 					echo -e "$sample\tAKA-017\t$tissue"
# 				;;
# 				M*)
# 					echo -e "$sample\tMUN-016\t$tissue"
# 				;;
# 				J*)
# 					echo -e "$sample\tJUT-011\t$tissue"
# 			esac
# 	done < <(cat $DIR/${tissue}Data.txt)
# done > $DIR/samples.txt

# cut -f2,3 $DIR/samples.txt | sort -u  > $DIR/data.txt

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/samples.txt) | cut -f 1)
tissue=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/samples.txt) | cut -f 3)

read1=$(echo "${DIRDATA}/RNA-seq/${tissue}/${sample}_read1.fastq.gz")
read2=$(echo "${DIRDATA}/RNA-seq/${tissue}/${sample}_read2.fastq.gz")

mkdir -p ${DIRDATA}/report_fastp/${sample}

fastp -w 12 -h ${DIRDATA}/report_fastp/${sample}/${sample}.html \
-j ${DIRDATA}/report_fastp/${sample}/${sample}.json \
-i ${read1} \
-I ${read2} \
-o ${DIRDATA}/fastp/${sample}_trim_1.fq.gz \
-O ${DIRDATA}/fastp/${sample}_trim_2.fq.gz \
-l 20 -q 20
