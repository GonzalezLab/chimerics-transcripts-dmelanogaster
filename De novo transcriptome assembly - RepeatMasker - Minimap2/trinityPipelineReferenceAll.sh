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

#mkdir DataAssemblyMergedReference

#cp ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/assemblyImproved/stringtie_merged.fasta DataAssemblyMergedReference/
#cp ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/assemblyImproved/mergedSamples.stringtie_merged.gtf.tmap DataAssemblyMergedReference/

#cd DataAssemblyMergedReference
#makeblastdb -in stringtie_merged.fasta -dbtype nucl
#cd ..

for assembly in `less ${DIR}/DATA/genomeAssembly/genomeAssembly.txt`
	do
		for tissueData in $tissues
			do
				tissue=$(echo $tissueData |cut -f1 -d"=")
				#mkdir Trinity_${assembly}_${tissue}
				cd Trinity_${assembly}_${tissue}
				#Trinity --seqType fq --max_memory 78G --trimmomatic --output output_Trinity_${assembly}_${tissue}  --samples_file ../samples.${assembly}.${tissue}.txt --CPU 12
				#TrinityStats.pl output_Trinity_${assembly}_${tissue}/Trinity.fasta > Trinity.stats
				
				#mkdir RepeatMasker_Blast_Merged_Reference
				
				cd RepeatMasker_Blast_Merged_Reference 
				

				#blastn -query ../output_Trinity_${assembly}_${tissue}/Trinity.fasta \
				#       -db ../../DataAssemblyMergedReference/stringtie_merged.fasta \
				#       -out Trinity_${assembly}_${tissue}_refTrans.blastn \
				#       -evalue 1e-20 -dust no -task megablast -num_threads 12 \
				#       -max_target_seqs 1 -outfmt 6

				#analyze_blastPlus_topHit_coverage.pl \
				#       Trinity_${assembly}_${tissue}_refTrans.blastn \
				#       ../output_Trinity_${assembly}_${tissue}/Trinity.fasta  \
				#       ../../DataAssemblyMergedReference/stringtie_merged.fasta

				 #cut -f2,3 Trinity_${assembly}_${tissue}_refTrans.blastn.hist.list > Trinity_${assembly}_${tissue}_refTrans.blastn.hist.list.names
				 # grep "Bin_100" Trinity_${assembly}_${tissue}_refTrans.blastn.hist.list | cut -f2,3  > Trinity_${assembly}_${tissue}_refTrans.blastn.hist_bin100.list.names
				 # grep "Bin_100" Trinity_${assembly}_${tissue}_refTrans.blastn.hist.list | cut -f2 > names_bin100_Trinity.txt
  			# 	 seqtk subseq ../output_Trinity_${assembly}_${tissue}/Trinity.fasta names_bin100_Trinity.txt > Trinity.hits_bin100.fasta
				 # rm -rf RepeatMaskerImproved
				 # RepeatMasker -pa 12 -e rmblast -lib ../../consensuses_curated_v4.fasta -dir RepeatMaskerImproved -norna -nolow -s -cutoff 250 -xsmall -no_is -gff Trinity.hits_bin100.fasta
				 # grep Unspecified RepeatMaskerImproved/Trinity.hits_bin100.fasta.out | sed 's/  */ /g'  | sed "s/^ //g" | tr ' ' '\t' | sort -nr -k1  > RepeatMaskerImproved/Trinity.hits_bin100.fasta.out.TE
     #  			 grep "Bin_100" Trinity_${assembly}_${tissue}_refTrans.blastn.hist.list | cut -f2,3  > Trinity_${assembly}_${tissue}_refTrans.blastn.hist_bin100.list.names
				 grep -E "Bin_90|Bin_100" Trinity_${assembly}_${tissue}_refTrans.blastn.hist.list | cut -f2,3  > Trinity_${assembly}_${tissue}_refTrans.blastn.hist_bin90.list.names
				 grep -E "Bin_90|Bin_100" Trinity_${assembly}_${tissue}_refTrans.blastn.hist.list | cut -f2 > names_bin90_Trinity.txt
  				 seqtk subseq ../output_Trinity_${assembly}_${tissue}/Trinity.fasta names_bin90_Trinity.txt > Trinity.hits_bin90.fasta
				 mkdir RepeatMaskerImproved90
				 RepeatMasker -pa 12 -e rmblast -lib ../../consensuses_curated_v4.fasta -dir RepeatMaskerImproved90 -norna -nolow -s -cutoff 250 -xsmall -no_is -gff Trinity.hits_bin90.fasta
				 grep Unspecified RepeatMaskerImproved90/Trinity.hits_bin90.fasta.out | sed 's/  */ /g'  | sed "s/^ //g" | tr ' ' '\t' | sort -nr -k1  > RepeatMaskerImproved90/Trinity.hits_bin90.fasta.out.TE

				cd ../..

			done
	done

assemblies="AKA-017 SLA-001 JUT-011 MUN-016 TOM-007"
tissues="head gut ovary"

# # for assembly in `less ${DIR}/DATA/genomeAssembly/genomeAssembly.txt`
# # do
# # for tissue in $tissues
# # do
# # cut -f1,2,3,4 -d'_' Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/names_Trinity.txt |sort |uniq -c | sed 's/^ *//' |tr -d '>' > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/transcripts_by_gene_analyzed.txt
# # done
# # done

# # for assembly in `less ${DIR}/DATA/genomeAssembly/genomeAssembly.txt`
# # do
# # for tissue in $tissues
# # do
# # 	echo $assembly $tissue
# #  cut -f5,10 Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMasker/Trinity.hits.fasta.out.TE |sed "s/_i.*\t/\t/g" | sort -u | cut -f2 | sort |uniq -c | sed 's/^ *//'  |sort -k1hr | head -n5
# # echo
# # done
# # done

# # for assembly in $assemblies
# # do
# # for tissue in $tissues
# # do
# # cut -f5 Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMasker/Trinity.hits.fasta.out.TE | cut -f1,2,3,4 -d'_' | sort -u > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMasker/genes_with_TE.txt
# # done

# # for assembly in $assemblies
# # do
# # for tissue in $tissues
# # do
# # cut -f5 Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMasker/Trinity.hits.fasta.out.TE | sort -u > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMasker/transcripts_with_TE.txt
# # done

# # for assembly in $assemblies
# # do
# # for tissue in $tissues
# # do
# # 	echo $assembly $tissue
	
# # 	cat Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/names_Trinity.txt | cut -f1,2,3,4 -d'_' | sort -u > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/names_gene_Trinity.txt
# # 	done
# # done

# # for assembly in $assemblies
# # 	do
# # 		for tissue in $tissues
# # 			do
# # 				echo $assembly $tissue
# # 				grep Unclassified Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMasker/Trinity.hits.fasta.tbl
# # 			done
# # 	done

# # for assembly in $assemblies
# # do
# # for tissue in $tissues
# # do
# # echo $assembly $tissue
# # echo -n "Number of transcripts analyzed: "
# # cat Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/names_Trinity.txt |wc -l
# # echo -n "Number of genes analyzed: "
# # #cat Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/names_Trinity.txt |cut -f1,2,3,4 -d'_' | sort -u| wc -l
# # cat Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/names_gene_Trinity.txt |wc -l
# # done
# # done


# # for assembly in `less ${DIR}/DATA/genomeAssembly/genomeAssembly.txt`
# # do
# # for tissue in $tissues
# # do
# # echo $assembly $tissue
# # echo -n "Number of genes: "
# # cat Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMasker/genes_with_TE.txt | wc -l 
# # echo -n "Number of transcripts: "
# # cut -f5 Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMasker/Trinity.hits.fasta.out.TE |sort -u | wc -l
# # join -1 1 -2 2 Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMasker/genes_with_TE.txt <(sort -k2 Trinity_${assembly}_${tissue}/transcripts_by_gene_Trinity.txt)  | awk '{ total += $2; count++ } END { print total/count }'
# # join -1 1 -2 2 Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/names_gene_Trinity.txt <(sort -k2 Trinity_${assembly}_${tissue}/transcripts_by_gene_Trinity.txt)  | awk '{ total += $2; count++ } END { print total/count }'

# # done
# # done

# # # he de sacar de names_gene_Trinity los k no tengan TE



# # for assembly in $assemblies
# # do
# # for tissue in $tissues
# # do
# # 	echo $assembly $tissue
	
# # 	cat Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/names_Trinity.txt > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/names_gene_Trinity.txt
# # 	done
# # done

for assembly in $assemblies
do
for tissue in $tissues
do
echo $assembly $tissue
echo -n "Number of transcripts with TE: "
cut -f5 Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/Trinity.hits_bin90.fasta.out.TE |sort -u | wc -l
	
transcripts_TE=$(cut -f5 Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/Trinity.hits_bin90.fasta.out.TE |sort -u | wc -l)
echo -n "Number of transcripts analyzed: "
cat Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/names_bin90_Trinity.txt |wc -l
transcripts_analyzed=$(cat Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/names_bin90_Trinity.txt |wc -l )
proportionTranscripts=$(echo "scale=6; $transcripts_TE/$transcripts_analyzed*100" | bc)
echo -n "Proportion transcripts: "
echo $proportionTranscripts

echo -n "Number of genes with TE: "
cut -f5 Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/Trinity.hits_bin90.fasta.out.TE | cut -f1,2,3,4 -d'_' |sort -u | wc -l
genes_TE=$(cut -f5 Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/Trinity.hits_bin90.fasta.out.TE | cut -f1,2,3,4 -d'_' |sort -u | wc -l)
	
echo -n "Number of genes analyzed: "
cat Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/names_bin90_Trinity.txt | cut -f1,2,3,4 -d'_' | sort -u > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/names_gene_bin90_Trinity.txt
cat Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/names_gene_bin90_Trinity.txt |wc -l
genes_analyzed=$(cat Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/names_gene_bin90_Trinity.txt |wc -l )
proportionGenes=$(echo "scale=6; $genes_TE/$genes_analyzed*100" | bc)
echo -n "Proportion genes: "
echo $proportionGenes
echo ""
done
done

for assembly in $assemblies
do
for tissue in $tissues
do
echo $assembly $tissue
grep Unclassified Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/Trinity.hits_bin90.fasta.tbl
done
done


for assembly in $assemblies
do
for tissue in $tissues
do
cut -f5 Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/Trinity.hits_bin90.fasta.out.TE | sort -u > Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/RepeatMaskerImproved90/transcripts_with_TE.txt
done
done