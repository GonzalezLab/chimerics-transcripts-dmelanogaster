#!/bin/bash

# set the partition where the job will run

#SBATCH --partition=normal

# set the number of nodes

#SBATCH --cpus-per-task=12

#SBATCH --mem=78G

# mail alert at start, end and abortion of execution

#SBATCH --mail-type=ALL

#SBATCH --job-name=chipSeq_JUT-011_gut_H3K9me3

#SBATCH --output=chipSeq_JUT-011_gut_H3K9me3.output

#SBATCH --error=chipSeq_JUT-011_gut_H3K9me3.error

# send mail to this address

#SBATCH --mail-user=marta.coronado@ibe.upf-csic.es

module load Miniconda3/4.7.10
module load Python/3.7.2-GCCcore-8.2.0
module load Java/12.0.2

source activate /homes/users/mcoronado/.conda/envs/encode-chip-seq-pipeline


INPUT=JUT-011_gut_H3K9me3.json

sbatch -p normal -J chipSeq_JUT-011_gut_H3K9me3 \
	--error=chipSeq_JUT-011_gut_H3K9me3.error \
	--output=chipSeq_JUT-011_gut_H3K9me3.output \
	--export=ALL \
	--cpus-per-task=12 --mem 78G -t 4-0 \
	--wrap "caper run /homes/users/mcoronado/bin/chip-seq-pipeline2/chip.wdl -i ${INPUT}"
