#!/bin/bash

# set the partition where the job will run

#SBATCH --partition=normal

# set the number of nodes

# mail alert at start, end and abortion of execution

#SBATCH --mail-type=ALL

#SBATCH --cpus-per-task=1

#SBATCH --mem=8G

#SBATCH --job-name=PFAMScan

#SBATCH --output=PFAMScan.output

#SBATCH --error=PFAMScan.error

# send mail to this address

#SBATCH --mail-user=marta.coronado@ibe.upf-csic.es

module load Miniconda3/4.7.10

source activate /homes/users/mcoronado/.conda/envs/trinity

# Run

DIR=/homes/users/mcoronado/scratch/5GenomesProject
assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"
tissues="head gut ovary"
#assemblies="AKA-017"
#tissues="head"

> globalResult.tsv

for assembly in $assemblies
	do
		for tissue in $tissues
			do
				#rm -r $assembly/$tissue/tmp/gff
				#rm -r $assembly/$tissue/tmp/fasta
				#rm $assembly/$tissue/*output
				mkdir -p $assembly/$tissue/tmp/gff
				mkdir -p $assembly/$tissue/tmp/fasta
				mkdir -p $assembly/$tissue/pfamOutput
				insertions=$(grep "Middle exon (TE inside)" ../TE_chimeric_global_REVISED_v2.tab | grep ${assembly} | grep ${tissue})
				while IFS= read -r insertion
					do
						transcript=$(echo "$insertion" | cut -f4)
						echo $tissue $assembly $transcript
						stringtieID=$(echo "$insertion" | cut -f6)
						geneID=$(echo "$insertion" | cut -f8)
						startTE=$(echo "$insertion" | cut -f10 | cut -f1 -d'-')
						endTE=$(echo "$insertion" | cut -f10 | cut -f2 -d'-')
						TEconsensusID=$(echo "$insertion" | cut -f13 | cut -f1 -d'(')
						strand=$(echo "$insertion" | cut -f13 | cut -f2 -d'(' |  head -c 1)
						TEfam=$(echo "$insertion" | cut -f14)
						#echo -e "$transcript\t$stringtieID\t$geneID\t$startTE\t$endTE\t$TEconsensusID\t$strand\t$TEfam"
						echo -e "$transcript\tfeature\tTE\t$startTE\t$endTE\t.\t$strand\t.\tname=$transcript" > $assembly/$tissue/tmp/gff/${transcript}_${TEfam}_${startTE}-${endTE}.gff
						bedtools getfasta -fi ../../Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/minimap2assembly/RepeatMaskerTranscriptsCorrectOrientation/transcriptsCorrectOrientation.fa.masked -bed $assembly/$tissue/tmp/gff/${transcript}_${TEfam}_${startTE}-${endTE}.gff -s  > $assembly/$tissue/tmp/fasta/${transcript}_${TEfam}_${startTE}-${endTE}.fasta 
						getorf -sequence $assembly/$tissue/tmp/fasta/${transcript}_${TEfam}_${startTE}-${endTE}.fasta  -outseq $assembly/$tissue/tmp/fasta/${transcript}_${TEfam}_${startTE}-${endTE}.ORF.fasta
						faidx  --transform chromsizes $assembly/$tissue/tmp/fasta/${transcript}_${TEfam}_${startTE}-${endTE}.ORF.fasta | sort -k2,2nr | head -n1 |cut -f1 > $assembly/$tissue/tmp/fasta/longestIsoform_${transcript}_${TEfam}_${startTE}-${endTE}.lst
						seqtk subseq $assembly/$tissue/tmp/fasta/${transcript}_${TEfam}_${startTE}-${endTE}.ORF.fasta $assembly/$tissue/tmp/fasta/longestIsoform_${transcript}_${TEfam}_${startTE}-${endTE}.lst > $assembly/$tissue/tmp/fasta/${transcript}_${TEfam}_${startTE}-${endTE}.longORF.fasta
						sed -i "1 s/./>${transcript}_${tissue}_${assembly}_${stringtieID}_${geneID}_${TEfam}_g2_/" $assembly/$tissue/tmp/fasta/${transcript}_${TEfam}_${startTE}-${endTE}.longORF.fasta
						pfam_scan.pl -fasta $assembly/$tissue/tmp/fasta/${transcript}_${TEfam}_${startTE}-${endTE}.longORF.fasta  -dir pfamFiles -outfile $assembly/$tissue/pfamOutput/${transcript}_${TEfam}_${startTE}-${endTE}.g2.pfam.output
				done <<< "$insertions"

					
			done
	done
