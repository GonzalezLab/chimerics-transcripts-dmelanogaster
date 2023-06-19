#!/bin/bash

# set the partition where the job will run

#SBATCH --partition=normal

# set the number of nodes

#SBATCH --cpus-per-task=12

#SBATCH --mem=64G

# mail alert at start, end and abortion of execution

#SBATCH --mail-type=ALL

#SBATCH --job-name=bowtie2Index

#SBATCH --output=bowtie2Index.output

#SBATCH --error=bowtie2Index.error

# send mail to this address

#SBATCH --mail-user=marta.coronado@ibe.upf-csic.es

module load Miniconda3/4.7.10

source activate /homes/users/mcoronado/.conda/envs/encode-chip-seq-pipeline

# Run

for i in `less genomeAssembly.txt`
    do
    	cd $i
    	faidx ${i}.fasta -i chromsizes > ${i}.chrom.sizes
    	cat ${i}.chrom.sizes | awk '{s+=$2}END{print s}' > ${i}.total.chrom.sizes
    	# mkdir bowtieIndex_${i}
    	# cd bowtieIndex_${i}
    	# bowtie-build ../${i}.fasta ${i}
    	mkdir bowtie2Index_${i}
    	cd bowtie2Index_${i}
    	bowtie2-build ../${i}.fasta ${i}
    	tar -czvf ${i}_bowtie2_index.tar.gz .
    	cd ../..
done
