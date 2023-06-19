#!/bin/bash

# set the partition where the job will run

#SBATCH --partition=normal

# set the number of nodes

#SBATCH --cpus-per-task=12

#SBATCH --mem=64G

# mail alert at start, end and abortion of execution

#SBATCH --mail-type=ALL

#SBATCH --job-name=mappingRepeatMaskerImproved/-017

#SBATCH --output=mappingJUT-011.output

#SBATCH --error=mappingJUT-011.error

# send mail to this address

#SBATCH --mail-user=marta.coronado@ibe.upf-csic.es

module load R/3.5.1-foss-2018b
module load csem/2.3-Perl-5.26.1-foss-2018a
module load SRA-Toolkit/2.8.2-1-centos_linux64
module load FastQC/0.11.7-Java-1.8.0_74
module load BBMap/38.00
module load SAMtools/1.6-foss-2016b
module load Bowtie/1.2.2-foss-2016b
module load Trimmomatic/0.36-Java-1.8.0_92


histones="H3K27ac H3K9me3 input"
tissues="head gut ovary"
#assemblies="JUT-011 JUT-011 MUN-016 SLA-001 TOM-007"

for tissue in $tissues
	do
		mkdir mapping/JUT-011
		mkdir mapping/JUT-011/${tissue}
		for histone in $histones
			do
				mkdir mapping/JUT-011/${tissue}/${histone}
				Rscript permSeq.R "JUT-011" $tissue $histone
			done
	done
