#!/bin/bash

# set the partition where the job will run

#SBATCH --partition=normal

# set the number of nodes

# mail alert at start, end and abortion of execution

#SBATCH --mail-type=ALL

#SBATCH --cpus-per-task=2

#SBATCH --mem=8G

#SBATCH --job-name=TE_length

#SBATCH --output=TE_length.output

#SBATCH --error=TE_length.error

# send mail to this address

#SBATCH --mail-user=marta.coronado@ibe.upf-csic.es

module load Miniconda3/4.7.10

source activate /homes/users/mcoronado/.conda/envs/trinity

# Run

DIR="/homes/users/mcoronado/scratch/5GenomesProject"
assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"
tissues="head gut ovary"
#assemblies="AKA-017"
#tissues="head"

while read TEfamily
do
tissue=$(awk -v TEfamily="$TEfamily" ' $2 == TEfamily ' tissue_family.tab | cut -f1 | tr '\n' ';' | sed 's/.$//')
echo -e "$TEfamily\t$tissue"
done < <(tail -n +2 ${DIR}/DATA/TEs/TE_library_family.csv| cut -f 8 | sort -u)
