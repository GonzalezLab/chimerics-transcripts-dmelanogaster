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

				#minimap2 -ax splice ${DIR}/DATA/FlyBase_r6.31/dmel-chr.fasta \
				#Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/Trinity.hits_bin90.TE.fasta >  \
				#Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2/alignment.sam

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
				#RepeatMasker -pa 12 -e rmblast -lib consensuses_curated_v4.fasta \
				#-dir Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts \
				#-norna -nolow -s -cutoff 250 -xsmall -no_is -gff  Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/transcripts.fasta 
				#grep Unspecified Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/transcripts.fasta.out | sed 's/  */ /g'  | sed "s/^ //g" | tr ' ' '\t' | sort -nr -k1  > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/transcripts.fasta.out.TE
				#bedtools intersect -a Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/annotation.bed \
				#-b ${DIR}/DATA/TEs/${assembly}.bed -wa -wb > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/intersectTEannot.noStrand.bed
				
				#bedtools window -a Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/annotationTranscripts.bed -b <(cut -f1-7 ${DIR}/DATA/TEs/${assembly}.bed ) -sm > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/intersectTEtranscript.bed
				#> Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/transcriptsFilter.fasta.out.TE
				#> Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/transcriptsFilter.fasta.out.gff
				#mkdir Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSam
				#samtools faidx Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/transcripts.fasta
				#samtools faidx Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/Trinity.hits_bin90.fasta
			# 	> Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/cleanIntersectTranscripts/totalIntersectTranscripts.gff
			# 	> Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/cleanIntersectTranscripts/exonOverlapTE.tmp.gff
			# 	> Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/cleanIntersectTranscripts/TEoverlapExon.tmp.gff

			# 	#mkdir Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/cleanIntersectTranscripts
			# 	echo "Analyzing transcripts of ${assembly}_${tissue}"
			# 	for transcript in `less  Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/uniqueID.txt`
			# 	do
			# 	#TEs=$(grep -w $transcript Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/Trinity.hits_bin90.fasta.out.TE | cut -f 10)
			# 	#grep -w $transcript Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/transcripts.fasta.out.TE | grep -w "$TEs" >> Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/transcriptsFilter.fasta.out.TE
			# 	#grep -w $transcript Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/transcripts.fasta.out.gff | grep -w "$TEs" >> Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/transcriptsFilter.fasta.out.gff
			# 	#samtools faidx Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/transcripts.fasta $transcript > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/$transcript.tmp.fa
			# 	#samtools faidx Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/Trinity.hits_bin90.fasta $transcript > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/${transcript}.expressed.tmp.fa
			# 	#mkdir Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSam/${transcript}
			# 	#minimap2 -ax splice --secondary=no --sam-hit-only -C5 -t4 Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/$transcript.tmp.fa \
			# 	#Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/${transcript}.expressed.tmp.fa >  \
			# 	#Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSam/${transcript}/${transcript}_alignment.sam
			# 	#bedtools bamtobed -split -i  Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSam/${transcript}/${transcript}_alignment.sam > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSam/${transcript}/${transcript}_alignment.bed
			# 	nAlign=$(samtools view -c Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSam/${transcript}/${transcript}_alignment.sam) 
			# 	if [[ $nAlign == 1 ]]; then
			# 		echo $transcript >> Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/uniqueID_after2nMapping.txt
			# 	fi
			# 	done

			# for transcript in `less  Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/uniqueID_after2nMapping.txt`
			# 	do
			# 	bed_to_gff3_converter.py Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSam/${transcript}/${transcript}_alignment.bed Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSam/${transcript}/${transcript}_alignment.tmp.gff
			# 	cat Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSam/${transcript}/${transcript}_alignment.tmp.gff | awk '!a[$0]++' > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSam/${transcript}/${transcript}_alignment.gff
			# 	#rm Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/$transcript.tmp.fa
			# 	#rm Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/${transcript}.expressed.tmp.fa
			# 	bedtools intersect -a Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSam/${transcript}/${transcript}_alignment.gff \
			# 	 -b Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/transcriptsFilter.fasta.out.gff -wa -wb >> \
			# 	  Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/cleanIntersectTranscripts/totalIntersectTranscripts.gff
				

			# 	bedtools intersect -a Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSam/${transcript}/${transcript}_alignment.gff \
			# 	 -b Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/transcriptsFilter.fasta.out.gff -wa -wb -f 1 >> \
			# 	  Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/cleanIntersectTranscripts/exonOverlapTE.tmp.gff


			# 	bedtools intersect -a Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSam/${transcript}/${transcript}_alignment.gff \
			# 	 -b Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/transcriptsFilter.fasta.out.gff -wa -wb -F 1 >> \
			# 	  Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/cleanIntersectTranscripts/TEoverlapExon.tmp.gff

				
				
			# 	done
			# 	echo "Finishing transcripts of ${assembly}_${tissue}"
			# 	echo "Creating partialOverlap.gff of ${assembly}_${tissue}"
			# 	grep -v -f Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/cleanIntersectTranscripts/exonOverlapTE.tmp.gff  Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/cleanIntersectTranscripts/totalIntersectTranscripts.gff | grep  -v -f Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/cleanIntersectTranscripts/TEoverlapExon.tmp.gff | sort -u > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/cleanIntersectTranscripts/partialOverlap.gff 
			# 	echo "Creating totalOverlapTEexon.gff of ${assembly}_${tissue}"
			# 	comm -12 <(sort Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/cleanIntersectTranscripts/exonOverlapTE.tmp.gff) <(sort Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/cleanIntersectTranscripts/TEoverlapExon.tmp.gff) | sort -u > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/cleanIntersectTranscripts/totalOverlapTEexon.gff 
			# 	echo "Creating final TEoverlapExon.gff"
			# 	grep -v -f Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/cleanIntersectTranscripts/totalOverlapTEexon.gff Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/cleanIntersectTranscripts/exonOverlapTE.tmp.gff | sort -u > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/cleanIntersectTranscripts/exonOverlapTE.gff
			# 	echo "Creating final exonOverlapTE.gff"
			# 	grep -v -f Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/cleanIntersectTranscripts/totalOverlapTEexon.gff Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/cleanIntersectTranscripts/TEoverlapExon.tmp.gff  | sort -u > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/cleanIntersectTranscripts/TEoverlapExon.gff



			# 	bedtools intersect -a Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSam/${transcript}/${transcript}_alignment.gff -b Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/transcriptsFilter.fasta.out.gff -wa -wb -s >  Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscripts/intersectTranscripts/${transcript}_strand.gff

			# 	bedtools intersect -a Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/annotation.tmp.bed \
			# 	-b ${DIR}/DATA/TEs/${assembly}.bed -s -wa -wb > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/intersectTEannot.tmp.bed

			# 	bedtools intersect -a Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/annotation.bed \
			# 	-b ${DIR}/DATA/TEs/${assembly}.bed -s -wa -wb -f 0.8 -r  > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/intersectTEannotFull.bed

			done
	done


# for assembly in `less ${DIR}/DATA/genomeAssembly/genomeAssembly.txt`
# do
# for tissueData in $tissues
# do
# tissue=$(echo $tissueData |cut -f1 -d"=")
# nT=$(cut -f4 Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/intersectTEannot.bed | sort -u |wc -l)
# echo -e "$assembly\t$tissue\t$nT"
# done
# done
