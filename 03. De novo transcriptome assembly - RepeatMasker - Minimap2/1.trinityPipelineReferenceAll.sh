#!/bin/bash

# set the partition where the job will run

#SBATCH --partition=normal

# set the number of nodes

#SBATCH --cpus-per-task=12

#SBATCH --mem=78G

# mail alert at start, end and abortion of execution

#SBATCH --mail-type=ALL

#SBATCH --job-name=trinityPipelineAll

#SBATCH --output=trinityPipelineAll.output

#SBATCH --error=trinityPipelineAll.error

# send mail to this address

#SBATCH --mail-user=marta.coronado@ibe.upf-csic.es

module load Miniconda3/4.7.10

source activate /homes/users/mcoronado/.conda/envs/trinity

# Run

DIR=/homes/users/mcoronado/scratch/5GenomesProject

#tissues="head=Heads gut=Gut ovary=Ovary"
tissues="head gut ovary"

mkdir DataAssemblyMergedReference

cp ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/assemblyImproved/stringtie_merged.fasta DataAssemblyMergedReference/
cp ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/assemblyImproved/mergedSamples.stringtie_merged.gtf.tmap DataAssemblyMergedReference/

cd DataAssemblyMergedReference
makeblastdb -in stringtie_merged.fasta -dbtype nucl
cd ..

for assembly in `less ${DIR}/DATA/genomeAssembly/genomeAssembly.txt`
	do
		for tissueData in $tissues
			do
				tissue=$(echo $tissueData |cut -f1 -d"=")
				mkdir Trinity_${assembly}_${tissue}
				cd Trinity_${assembly}_${tissue}
				Trinity --seqType fq --max_memory 78G --trimmomatic --output output_Trinity_${assembly}_${tissue}  --samples_file ../samples.${assembly}.${tissue}.txt --CPU 12
				TrinityStats.pl output_Trinity_${assembly}_${tissue}/Trinity.fasta > Trinity.stats
				
				mkdir RepeatMasker_Blast_Merged_Reference
				
				cd RepeatMasker_Blast_Merged_Reference 
				

				blastn -query ../output_Trinity_${assembly}_${tissue}/Trinity.fasta \
				      -db ../../DataAssemblyMergedReference/stringtie_merged.fasta \
				      -out Trinity_${assembly}_${tissue}_refTrans.blastn \
				      -evalue 1e-20 -dust no -task megablast -num_threads 12 \
				      -max_target_seqs 1 -outfmt 6

				analyze_blastPlus_topHit_coverage.pl \
				      Trinity_${assembly}_${tissue}_refTrans.blastn \
				      ../output_Trinity_${assembly}_${tissue}/Trinity.fasta  \
				      ../../DataAssemblyMergedReference/stringtie_merged.fasta

				 cut -f2,3 Trinity_${assembly}_${tissue}_refTrans.blastn.hist.list > Trinity_${assembly}_${tissue}_refTrans.blastn.hist.list.names
				 #grep "Bin_100" Trinity_${assembly}_${tissue}_refTrans.blastn.hist.list | cut -f2,3  > Trinity_${assembly}_${tissue}_refTrans.blastn.hist_bin100.list.names
				 #grep "Bin_100" Trinity_${assembly}_${tissue}_refTrans.blastn.hist.list | cut -f2 > names_bin100_Trinity.txt
  				 #seqtk subseq ../output_Trinity_${assembly}_${tissue}/Trinity.fasta names_bin100_Trinity.txt > Trinity.hits_bin100.fasta
				 #rm -rf RepeatMaskerImproved
				 #RepeatMasker -pa 12 -e rmblast -lib ../../consensuses_curated_v4.fasta -dir RepeatMaskerImproved -norna -nolow -s -cutoff 250 -xsmall -no_is -gff Trinity.hits_bin100.fasta
				 #grep Unspecified RepeatMaskerImproved/Trinity.hits_bin100.fasta.out | sed 's/  */ /g'  | sed "s/^ //g" | tr ' ' '\t' | sort -nr -k1  > RepeatMaskerImproved/Trinity.hits_bin100.fasta.out.TE
      			 #grep "Bin_100" Trinity_${assembly}_${tissue}_refTrans.blastn.hist.list | cut -f2,3  > Trinity_${assembly}_${tissue}_refTrans.blastn.hist_bin100.list.names
				 grep -E "Bin_90|Bin_100" Trinity_${assembly}_${tissue}_refTrans.blastn.hist.list | cut -f2,3  > Trinity_${assembly}_${tissue}_refTrans.blastn.hist_bin90.list.names
				 grep -E "Bin_90|Bin_100" Trinity_${assembly}_${tissue}_refTrans.blastn.hist.list | cut -f2 > names_bin90_Trinity.txt
  				 seqtk subseq ../output_Trinity_${assembly}_${tissue}/Trinity.fasta names_bin90_Trinity.txt > Trinity.hits_bin90.fasta
				 mkdir RepeatMaskerImproved90
				 RepeatMasker -pa 12 -e rmblast -lib ../../consensuses_curated_v4.fasta -dir RepeatMaskerImproved90 -norna -nolow -s -cutoff 250 -xsmall -no_is -gff Trinity.hits_bin90.fasta
				 grep Unspecified RepeatMaskerImproved90/Trinity.hits_bin90.fasta.out | sed 's/  */ /g'  | sed "s/^ //g" | tr ' ' '\t' | sort -nr -k1  > RepeatMaskerImproved90/Trinity.hits_bin90.fasta.out.TE

				cd ../..

			done
	done
