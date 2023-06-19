#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=clean_set_transcripts
#SBATCH --output=logs/clean_set_transcripts_%a.out
#SBATCH --error=logs/clean_set_transcripts_%a.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate chimerics
module load foss/2021b
module load R/4.1.2

DIR="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/chimerics"
DIRDATA="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/DATA"

cd ${DIR}

cat ${DIR}/clean_set_minimap2/*/chimerics.tmp.tab > ${DIR}/clean_set_minimap2/chimerics.tmp.tab 

Rscript unify_chimerics.R 
