#!/bin/bash

# set the partition where the job will run

#SBATCH --partition=normal

# set the number of nodes

# mail alert at start, end and abortion of execution

#SBATCH --mail-type=ALL

#SBATCH --cpus-per-task=1

#SBATCH --mem=8G

#SBATCH --job-name=blast	

#SBATCH --output=blast	.output

#SBATCH --error=blast	.error

# send mail to this address

#SBATCH --mail-user=marta.coronado@ibe.upf-csic.es


module load Miniconda3/4.7.10

source activate /homes/users/mcoronado/.conda/envs/trinity


makeblastdb -in final_library_v3_Dsim.fa  -dbtype nucl
makeblastdb -in roo_consensus.fasta  -dbtype nucl

blastn -query roo_consensus.fasta -db final_library_v3_Dsim.fa -outfmt 6 -qcov_hsp_perc 80 -perc_identity 80 > result_blast_Dsim_library.blastoutput

# get name: Dsim-B-P72.17-Map10#LTR-BelPao

cut -f2 result_blast_Dsim_library.blastoutput > id_name.txt

seqtk subseq final_library_v3_Dsim.fa id_name.txt > roo_Dsim.fa

makeblastdb -in roo_Dsim.fa -dbtype nucl

blastn -query repetitive_region_roo.fasta -db roo_Dsim.fa  -outfmt 6 -qcov_hsp_perc 80 -perc_identity 80  -dust no -soft_masking false -word_size 7 -max_target_seqs 1 -evalue 0.05 -gapopen 5 -gapextend 2   > result_blast_Dsim_library.blastoutput
