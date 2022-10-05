#!/bin/bash

# set the partition where the job will run

#SBATCH --partition=normal

# set the number of nodes

#SBATCH --cpus-per-task=12

#SBATCH --mem=32G

# mail alert at start, end and abortion of execution

#SBATCH --mail-type=ALL

#SBATCH --job-name=minimap

#SBATCH --output=minimap.output

#SBATCH --error=minimap.error

# send mail to this address

#SBATCH --mail-user=marta.coronado@ibe.upf-csic.es

module load Miniconda3/4.7.10

source activate /homes/users/mcoronado/.conda/envs/trinity

# Run

DIR=/homes/users/mcoronado/scratch/5GenomesProject

#tissues="head=Heads gut=Gut ovary=Ovary"
tissues="head gut ovary"

for assembly in `less ${DIR}/DATA/genomeAssembly/genomeAssembly.txt`
	do
		for tissueData in $tissues
			do
				tissue=$(echo $tissueData |cut -f1 -d"=")
				mkdir Trinity_${assembly}_${tissue}
				seqtk subseq Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/Trinity.hits_bin90.fasta \
				Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/transcripts_with_TE.txt > \
				Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/Trinity.hits_bin90.TE.fasta
				
				mkdir Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly

				minimap2 -ax splice ${DIR}/DATA/FlyBase_r6.31/dmel-chr.fasta \
				Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/Trinity.hits_bin90.TE.fasta >  \
				Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2/alignment.sam

				minimap2 -ax splice --secondary=no --sam-hit-only -C5 -t4 ${DIR}/DATA/genomeAssembly/${assembly}/${assembly}.fasta \
				Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/Trinity.hits_bin90.TE.fasta >  \
				Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/alignment.sam

				samtools view -h -Sq 1 Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/alignment.sam > \
				Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/alignment.nomapq0.sam


				bedtools bamtobed -split -i Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/alignment.nomapq0.sam > \
				Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/annotation.tmp.bed

				bedtools bamtobed -i Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/alignment.nomapq0.sam > \
				Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/annotationTranscripts.tmp.bed

				Rscript correctCoordinates.R $assembly $tissue

				cut -f1 Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/alignment.nomapq0.sam | grep -v "^@" | sort |uniq -c | grep " 1 " | sed 's/  */ /g'  | sed "s/^ //g" | tr ' ' '\t' | cut -f2 |sort -u  > \
				Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/uniqueID.txt

				grep -w -f Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/uniqueID.txt Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/annotation.tmp.bed > \
				Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/annotation.bed


				grep -w -f Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/uniqueID.txt Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/annotationTranscripts.tmp2.bed > \
				Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/annotationTranscripts.bed
			
				gffread -w Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/transcripts.fasta -g ${DIR}/DATA/genomeAssembly/${assembly}/${assembly}.fasta Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/annotationTranscripts.bed
			done
	done

