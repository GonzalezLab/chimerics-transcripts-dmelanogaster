#!/bin/bash

# set the partition where the job will run

#SBATCH --partition=normal

# set the number of nodes

#SBATCH --cpus-per-task=8

#SBATCH --mem=32G

# mail alert at start, end and abortion of execution

#SBATCH --mail-type=ALL

#SBATCH --job-name=minimap2

#SBATCH --output=minimap2.output

#SBATCH --error=minimap2.error

# send mail to this address

#SBATCH --mail-user=marta.coronado@ibe.upf-csic.es

module load Miniconda3/4.7.10

source activate /homes/users/mcoronado/.conda/envs/repeatmasker

# Run

DIR=/homes/users/mcoronado/scratch/5GenomesProject

#tissues="head=Heads gut=Gut ovary=Ovary"
tissues="head gut ovary"
assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"

#tissues="head"

#transcripts="TRINITY_DN1946_c0_g1_i2 TRINITY_DN1032_c2_g1_i13 TRINITY_DN1112_c0_g3_i3 TRINITY_DN1182_c5_g1_i5 TRINITY_DN5372_c0_g1_i4 TRINITY_DN1554_c0_g1_i15 TRINITY_DN1020_c1_g1_i1 TRINITY_DN1067_c0_g1_i7 TRINITY_DN145_c0_g2_i9"
for assembly in $assemblies
	do

		for tissue in $tissues
			do
				samtools faidx Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/transcripts.fasta
					rm -r Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSamCorrectOrientation/
					mkdir Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSamCorrectOrientation/
					> Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/transcriptsCorrectOrientation.fa
					> Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/transcriptsNotAnalyzed.tsv 
					for transcript in `less  Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/uniqueID_after2nMapping.txt`
						do
							mkdir  Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSamCorrectOrientation/${transcript}
							orientationSam=$(grep -w $transcript Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/annotationTranscripts.bed | cut -f 6 )
							
							transcriptGtf=$(grep -w $transcript Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/transcriptsTE.tsv | cut -f 4)

							transcriptClass=$(grep -w $transcript Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/transcriptsTE.tsv | cut -f 5)

							orientationGtf=$(grep -w $transcriptGtf ../GeneTransference/${assembly}/${assembly}.gtf | head -n1 | cut -f 7)
							
							echo $transcript $orientationSam $transcriptGtf $orientationGtf $transcriptClass

							if [ $transcriptClass  == "u"  ]; then
								echo not analyzed $transcript $transcriptGtf $transcriptClass
								echo $transcript $transcriptGtf $transcriptClass >> Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/transcriptsNotAnalyzed.tsv 
							elif [ $orientationSam == $orientationGtf -a $transcriptClass  != "u"  ]; then
								echo equal $transcript $transcriptClass 
								samtools faidx Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/transcripts.fasta $transcript >> Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/transcriptsCorrectOrientation.fa
							elif [ $orientationSam != $orientationGtf -a $transcriptClass  != "u" ]; then
								echo not equal $transcript $transcriptClass 
								samtools faidx Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/transcripts.fasta $transcript | seqkit seq -r -p -t dna -v  >> Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/transcriptsCorrectOrientation.fa
							else
								echo $transcript $transcriptGtf $transcriptClass >> Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/transcriptsNotAnalyzed.tsv
							fi
								echo "###############"
							samtools faidx Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/Trinity.hits_bin90.fasta $transcript > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/${transcript}.expressed.tmp.fa
							samtools faidx Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/transcriptsCorrectOrientation.fa
							samtools faidx Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/transcriptsCorrectOrientation.fa $transcript > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/$transcript.tmp.fa

							minimap2 -ax splice --secondary=no --sam-hit-only -C5 -t4 Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/$transcript.tmp.fa \
							Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/${transcript}.expressed.tmp.fa >  \
							Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSamCorrectOrientation/${transcript}/${transcript}_alignment.sam
							bedtools bamtobed -split -i  Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSamCorrectOrientation/${transcript}/${transcript}_alignment.sam > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSamCorrectOrientation/${transcript}/${transcript}_alignment.bed
							python3 ~/bin/bed_to_gff3/bed_to_gff3_converter.py Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSamCorrectOrientation/${transcript}/${transcript}_alignment.bed Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSamCorrectOrientation/${transcript}/${transcript}_alignment.tmp.gff
							cat Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSamCorrectOrientation/${transcript}/${transcript}_alignment.tmp.gff | awk '!a[$0]++' > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSamCorrectOrientation/${transcript}/${transcript}_alignment.gff
							rm Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/${transcript}.expressed.tmp.fa
							rm Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/$transcript.tmp.fa
					done

				 RepeatMasker -pa 8 -e rmblast -lib consensuses_curated_v4.fasta \
				 -dir Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation \
				 -norna -nolow -s -cutoff 250 -xsmall -no_is -gff  Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/transcriptsCorrectOrientation.fa
				 grep Unspecified Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/transcriptsCorrectOrientation.fa.out | sed 's/  */ /g'  | sed "s/^ //g" | tr ' ' '\t' | sort -nr -k1  > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/transcriptsCorrectOrientation.fa.out.TE

				 > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/transcriptsCorrectOrientationFilter.fasta.out.TE
				 > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/transcriptsCorrectOrientationFilter.fasta.out.gff
				 for transcript in `less  Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/uniqueID_after2nMapping.txt`
						do

							TEs=$(grep -w $transcript Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/Trinity.hits_bin90.fasta.out.TE | cut -f 10)
							grep -w $transcript Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/transcriptsCorrectOrientation.fa.out.TE | grep -w "$TEs" >> Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/transcriptsCorrectOrientationFilter.fasta.out.TE
							grep -w $transcript Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/transcriptsCorrectOrientation.fa.out.gff | grep -w "$TEs" >> Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/transcriptsCorrectOrientationFilter.fasta.out.gff

					done

				mkdir Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/cleanIntersectTranscripts
				cut -f 5 Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/transcriptsCorrectOrientationFilter.fasta.out.TE | sort -u > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/transcriptsAnalyzedFilter.txt
				 > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/cleanIntersectTranscripts/totalIntersectTranscripts.gff
				 > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/cleanIntersectTranscripts/exonOverlapTE.tmp.gff
				 > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/cleanIntersectTranscripts/TEoverlapExon.tmp.gff
				 

				

					for transcript in `less  Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/transcriptsAnalyzedFilter.txt`
						do

							bedtools intersect -a Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSamCorrectOrientation/${transcript}/${transcript}_alignment.gff \
						 -b Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/transcriptsCorrectOrientationFilter.fasta.out.gff -wa -wb >> \
						  Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/cleanIntersectTranscripts/totalIntersectTranscripts.gff
						

						bedtools intersect -a Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSamCorrectOrientation/${transcript}/${transcript}_alignment.gff \
						 -b Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/transcriptsCorrectOrientationFilter.fasta.out.gff -wa -wb -f 1 >> \
						  Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/cleanIntersectTranscripts/exonOverlapTE.tmp.gff


						bedtools intersect -a Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/tmpSamCorrectOrientation/${transcript}/${transcript}_alignment.gff \
						 -b Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/transcriptsCorrectOrientationFilter.fasta.out.gff -wa -wb -F 1 >> \
						  Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/cleanIntersectTranscripts/TEoverlapExon.tmp.gff



					done

				echo "Finishing transcripts of ${assembly}_${tissue}"
				echo "Creating partialOverlap.gff of ${assembly}_${tissue}"
				grep -v -f Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/cleanIntersectTranscripts/exonOverlapTE.tmp.gff  Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/cleanIntersectTranscripts/totalIntersectTranscripts.gff | grep  -v -f Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/cleanIntersectTranscripts/TEoverlapExon.tmp.gff | sort -u > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/cleanIntersectTranscripts/partialOverlap.gff 
				echo "Creating totalOverlapTEexon.gff of ${assembly}_${tissue}"
				comm -12 <(sort Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/cleanIntersectTranscripts/exonOverlapTE.tmp.gff) <(sort Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/cleanIntersectTranscripts/TEoverlapExon.tmp.gff) | sort -u > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/cleanIntersectTranscripts/totalOverlapTEexon.gff 
				echo "Creating final TEoverlapExon.gff"
				grep -v -f Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/cleanIntersectTranscripts/totalOverlapTEexon.gff Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/cleanIntersectTranscripts/exonOverlapTE.tmp.gff | sort -u > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/cleanIntersectTranscripts/exonOverlapTE.gff
				echo "Creating final exonOverlapTE.gff"
				grep -v -f Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/cleanIntersectTranscripts/totalOverlapTEexon.gff Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/cleanIntersectTranscripts/TEoverlapExon.tmp.gff  | sort -u > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/cleanIntersectTranscripts/TEoverlapExon.gff

			done
	done