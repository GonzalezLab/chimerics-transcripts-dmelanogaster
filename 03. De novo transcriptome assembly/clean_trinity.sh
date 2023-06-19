#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --job-name=clean
#SBATCH --output=logs/clean_%a.out
#SBATCH --error=logs/clean_%a.err

DIR="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/chimerics"
DIRDATA="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/DATA"

cd ${DIR}

strain=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/data.txt) | cut -f 1)
tissue=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/data.txt) | cut -f 2)

rm -rf ${DIR}/trinity/trinity_${strain}_${tissue}/scaffolding_entries.sam
rm -rf ${DIR}/trinity/trinity_${strain}_${tissue}/*fa
rm -rf ${DIR}/trinity/trinity_${strain}_${tissue}/insilico_read_normalization
rm -rf ${DIR}/trinity/trinity_${strain}_${tissue}/chrysalis
rm -rf ${DIR}/trinity/trinity_${strain}_${tissue}/read_partitions
rm -rf ${DIR}/trinity/trinity_${strain}_${tissue}/__salmon_filt.chkpts
rm -rf ${DIR}/trinity/trinity_${strain}_${tissue}/Trinity.tmp.fasta.salmon.idx

