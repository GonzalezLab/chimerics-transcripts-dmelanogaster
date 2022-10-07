#!/bin/bash

# set the partition where the job will run

#SBATCH --partition=normal

# set the number of nodes

# mail alert at start, end and abortion of execution

#SBATCH --mail-type=ALL

#SBATCH --cpus-per-task=1

#SBATCH --mem=8G

#SBATCH --job-name=salmonAnalysis

#SBATCH --output=salmonAnalysis.output

#SBATCH --error=salmonAnalysis.error

# send mail to this address

#SBATCH --mail-user=marta.coronado@ibe.upf-csic.es

module load Miniconda3/4.7.10

source activate /homes/users/mcoronado/.conda/envs/trinity

mkdir -p DATA/assembly
mkdir DATA/all
mkdir RESULTS

tissues="head gut ovary"
assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"
#tissues="head"
#assemblies="TOM-007"
DIR=/homes/users/mcoronado/scratch/5GenomesProject/ANALYSES/trinity-repeatMasker

# copy data
cp  ${DIR}/DataAssemblyMergedReference/mergedSamples.stringtie_merged.gtf.tmap DATA/assembly/mergedSamples.stringtie_merged.gtf.tmap
cut -f1,5  DATA/assembly/mergedSamples.stringtie_merged.gtf.tmap | sort -u > DATA/assembly/genesTranscript.gtf.tmap
> DATA/all/Trinity.all.TMM.EXPR.matrix

# for assembly in $assemblies
# do
# 	for tissue in $tissues
# 	do
# 	cat DATA/$assembly/$tissue/Trinity_${assembly}_${tissue}_refTrans.blastn.hist_bin90.list.names | sed "s/$/\t${tissue}\t$assembly/"   >> DATA/assembly/Trinity_refTrans.blastn.hist_bin90.list.names
# 	done

# 	done


for assembly in $assemblies
do
	for tissue in $tissues
	do
		echo $assembly - $tissue
		mkdir -p DATA/$assembly/$tissue
		cp  ${DIR}/Trinity_${assembly}_${tissue}/salmon_outdir/Trinity_trans.isoform.TMM.EXPR.matrix DATA/$assembly/$tissue
		cp  ${DIR}/Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/Trinity_${assembly}_${tissue}_refTrans.blastn.hist_bin90.list.names DATA/$assembly/$tissue			
		grep $assembly ${DIR}/resultsChimAnalysis/TE_chimeric_global_REVISED_v2.tab | grep $tissue > DATA/$assembly/$tissue/TE_chimeric_global_REVISED_${assembly}_${tissue}.tab
		> DATA/$assembly/$tissue/Trinity_${assembly}_${tissue}_listTranscriptGene.names
		while IFS= read -r info
		do
			transcript=$(echo "$info" | cut -f2)
			gene=$(grep -wF "${transcript}" DATA/assembly/genesTranscript.gtf.tmap | cut -f1)
			echo -e "$info\t$gene" >> DATA/$assembly/$tissue/Trinity_${assembly}_${tissue}_listTranscriptGene.names
		done < DATA/$assembly/$tissue/Trinity_${assembly}_${tissue}_refTrans.blastn.hist_bin90.list.names
		cut -f8 DATA/$assembly/$tissue/TE_chimeric_global_REVISED_${assembly}_${tissue}.tab | sort -u > DATA/$assembly/$tissue/genesTE.tab
		grep -w -f DATA/$assembly/$tissue/genesTE.tab DATA/$assembly/$tissue/Trinity_${assembly}_${tissue}_listTranscriptGene.names > DATA/$assembly/$tissue/transcriptsGeneTE.lst
		cut -f 4 DATA/$assembly/$tissue/TE_chimeric_global_REVISED_${assembly}_${tissue}.tab | sort -u >  DATA/$assembly/$tissue/transcriptTE.lst
		> DATA/$assembly/$tissue/transcriptsGeneTEStatus.lst
		while IFS= read -r info
		do
			transcript=$(echo "$info" | cut -f1)
			if grep -qwF "${transcript}" DATA/$assembly/$tissue/transcriptTE.lst
			then
				TE=present
			else
				TE=notPresent
			fi

			echo -e "$info\t$TE" >> DATA/$assembly/$tissue/transcriptsGeneTEStatus.lst
		done < DATA/$assembly/$tissue/transcriptsGeneTE.lst
		

		nGenes=$(cut -f3 DATA/$assembly/$tissue/transcriptsGeneTEStatus.lst | sort -u | wc -l)
		nTranscripts=$(cut -f1 DATA/$assembly/$tissue/transcriptsGeneTEStatus.lst |sort -u |wc -l )
		avgTotal=$( echo "scale=3;$nTranscripts/$nGenes" | bc -l |  awk '{printf "%.3f\n", $0}')
		nGenesTE=$(grep "notPresent" DATA/$assembly/$tissue/transcriptsGeneTEStatus.lst | cut -f3 |sort -u |wc -l)
		GenesTE=$(grep "notPresent" DATA/$assembly/$tissue/transcriptsGeneTEStatus.lst | cut -f3 |sort -u)
		ntranscriptsTE=$(grep -wF "${GenesTE}" DATA/$assembly/$tissue/transcriptsGeneTEStatus.lst | cut -f1 |sort -u|wc -l)
		transcriptsTE=$(grep -wF "${GenesTE}" DATA/$assembly/$tissue/transcriptsGeneTEStatus.lst | cut -f1 |sort -u)
		avgTE=$( echo "scale=3;$ntranscriptsTE/$nGenesTE" | bc -l |  awk '{printf "%.3f\n", $0}')
		nGenesAlwaysTE=$(grep -v -f <(grep notP DATA/$assembly/$tissue/transcriptsGeneTEStatus.lst | cut -f3 | sort -u) DATA/$assembly/$tissue/transcriptsGeneTEStatus.lst | cut -f3 |sort -u |wc -l)
		GenesAlwaysTE=$(grep -v -f <(grep notP DATA/$assembly/$tissue/transcriptsGeneTEStatus.lst | cut -f3 | sort -u) DATA/$assembly/$tissue/transcriptsGeneTEStatus.lst | cut -f3 |sort -u )
		nTranscriptsAlwaysTE=$(grep -wF "${GenesAlwaysTE}"  DATA/$assembly/$tissue/transcriptsGeneTEStatus.lst | cut -f1 |sort -u|wc -l)
		TranscriptsAlwaysTE=$(grep -wF "${GenesAlwaysTE}"  DATA/$assembly/$tissue/transcriptsGeneTEStatus.lst | cut -f1 |sort -u)
		avgAlwaysTE=$( echo "scale=3;$nTranscriptsAlwaysTE/$nGenesAlwaysTE" | bc -l |  awk '{printf "%.3f\n", $0}')
		
				# Expression matrix genes that produce transcripts with and without TE
				grep -wF "${transcriptsTE}"  DATA/$assembly/$tissue/Trinity_trans.isoform.TMM.EXPR.matrix > DATA/$assembly/$tissue/Trinity.TEisoforms.TMM.EXPR.tmp.matrix 
				join -j 1 <(sort -f -k 1 DATA/$assembly/$tissue/Trinity.TEisoforms.TMM.EXPR.tmp.matrix)  <(sort -f -k 1 DATA/$assembly/$tissue/transcriptsGeneTEStatus.lst) > DATA/$assembly/$tissue/Trinity.TEisoforms.TMM.EXPR.tmp2.matrix
				while IFS= read -r transcriptLine
				do
					transcript=$(echo "$transcriptLine" | cut -f1 -d' ' )
					group=$(grep -w "$transcript" ${DIR}/resultsChimAnalysis/TE_chimeric_global_REVISED_v2.tab | grep $tissue | grep $assembly | cut -f 24 | sort -u | tr '\n' ',' | sed 's/.$//')
					if [ -z "$group" ]
					then
						group="NA"
					fi
					transcriptID=$(echo "$transcriptLine" | cut -f5 -d' ' )
						#nTissues=$(grep -w "$transcriptID" ${DIR}/resultsChimAnalysis/TE_chimeric_global_REVISED_v2.tab | cut -f1 | sort -u | wc -l )
						#nStrains=$(grep -w "$transcriptID" ${DIR}/resultsChimAnalysis/TE_chimeric_global_REVISED_v2.tab | cut -f2 | sort -u | wc -l )
						#if [ "$nTissues" -eq 0 -a "$nStrains" -eq 0 ]
						#	then
						#		nTissues="NA"
						#		nStrains="NA"
						#fi
						echo $transcript $group #$nTissues $nStrains
					done < <(cat DATA/$assembly/$tissue/Trinity.TEisoforms.TMM.EXPR.tmp2.matrix) > DATA/$assembly/$tissue/Trinity.TEisoforms.TMM.EXPR.matrix.group
					join -j 1 <(sort -f -k 1 DATA/$assembly/$tissue/Trinity.TEisoforms.TMM.EXPR.tmp2.matrix)  <(sort -f -k 1 DATA/$assembly/$tissue/Trinity.TEisoforms.TMM.EXPR.matrix.group) > DATA/$assembly/$tissue/Trinity.TEisoforms.TMM.EXPR.matrix	
 				# Expression matrix genes that produce transcripts always with TE
 				grep -wF "${TranscriptsAlwaysTE}"  DATA/$assembly/$tissue/Trinity_trans.isoform.TMM.EXPR.matrix > DATA/$assembly/$tissue/Trinity.alwaysTEisoforms.TMM.EXPR.tmp.matrix 
 				join -j 1 <(sort -f -k 1 DATA/$assembly/$tissue/Trinity.alwaysTEisoforms.TMM.EXPR.tmp.matrix)  <(sort -f -k 1 DATA/$assembly/$tissue/transcriptsGeneTEStatus.lst) > DATA/$assembly/$tissue/Trinity.alwaysTEisoforms.TMM.EXPR.tmp2.matrix
 				sed -i  's/present/alwaysPresent/g' DATA/$assembly/$tissue/Trinity.alwaysTEisoforms.TMM.EXPR.tmp2.matrix
 				while IFS= read -r transcriptLine
 				do
 					transcript=$(echo "$transcriptLine" | cut -f1 -d' ')
 					group=$(grep -w "$transcript" ${DIR}/resultsChimAnalysis/TE_chimeric_global_REVISED_v2.tab | grep $tissue | grep $assembly | cut -f 24 | sort -u | tr '\n' ',' | sed 's/.$//')
 					if [ -z "$group" ]
 					then
 						group="NA"
 					fi
 					transcriptID=$(echo "$transcriptLine" | cut -f5 -d' ' )
						#nTissues=$(grep -w "$transcriptID" ${DIR}/resultsChimAnalysis/TE_chimeric_global_REVISED_v2.tab | cut -f1 | sort -u | wc -l )
						#nStrains=$(grep -w "$transcriptID" ${DIR}/resultsChimAnalysis/TE_chimeric_global_REVISED_v2.tab | cut -f2 | sort -u | wc -l )
						echo $transcript $group #$nTissues $nStrains
					done < <(cat DATA/$assembly/$tissue/Trinity.alwaysTEisoforms.TMM.EXPR.tmp2.matrix) > DATA/$assembly/$tissue/Trinity.alwaysTEisoforms.TMM.EXPR.matrix.group
					join -j 1 <(sort -f -k 1 DATA/$assembly/$tissue/Trinity.alwaysTEisoforms.TMM.EXPR.tmp2.matrix)  <(sort -f -k 1 DATA/$assembly/$tissue/Trinity.alwaysTEisoforms.TMM.EXPR.matrix.group) > DATA/$assembly/$tissue/Trinity.alwaysTEisoforms.TMM.EXPR.matrix

				# Expression matrix genes that produce transcripts never with TE
				genesNeverTE=$(comm -13 <(sort DATA/$assembly/$tissue/genesTE.tab) <(cut -f 3 DATA/$assembly/$tissue/Trinity_${assembly}_${tissue}_listTranscriptGene.names | grep -v '-' | sort -u))
				nGenesNeverTE=$(comm -13 <(sort DATA/$assembly/$tissue/genesTE.tab) <(cut -f 3 DATA/$assembly/$tissue/Trinity_${assembly}_${tissue}_listTranscriptGene.names | grep -v '-' | sort -u) | wc -l)
				transcriptsNeverTE=$(grep -wF "$genesNeverTE" DATA/$assembly/$tissue/Trinity_${assembly}_${tissue}_listTranscriptGene.names | cut -f1 |sort -u )
				ntranscriptsNeverTE=$(grep -wF "$genesNeverTE" DATA/$assembly/$tissue/Trinity_${assembly}_${tissue}_listTranscriptGene.names | cut -f1 |sort -u |wc -l )
				avNeverTE=$( echo "scale=3;$ntranscriptsNeverTE/$nGenesNeverTE" | bc -l |  awk '{printf "%.3f\n", $0}')

				echo "$transcriptsNeverTE" > DATA/$assembly/$tissue/transcriptsNeverTE.lst
				grep -wF -f "DATA/$assembly/$tissue/transcriptsNeverTE.lst" DATA/$assembly/$tissue/Trinity_trans.isoform.TMM.EXPR.matrix > DATA/$assembly/$tissue/Trinity.neverTEisoforms.TMM.EXPR.tmp.matrix 
				join -j 1 <(sort -f -k 1 DATA/$assembly/$tissue/Trinity.neverTEisoforms.TMM.EXPR.tmp.matrix)  <(sort -f -k 1 DATA/$assembly/$tissue/Trinity_${assembly}_${tissue}_listTranscriptGene.names) > DATA/$assembly/$tissue/Trinity.neverTEisoforms.TMM.EXPR.matrix
				sed -i "s/$/\tnever/" DATA/$assembly/$tissue/Trinity.neverTEisoforms.TMM.EXPR.matrix
				sed -i "s/$/\tNA/" DATA/$assembly/$tissue/Trinity.neverTEisoforms.TMM.EXPR.matrix
				sed -i "s/ /\t/g" DATA/$assembly/$tissue/Trinity.neverTEisoforms.TMM.EXPR.matrix
				# create file with three types
				cat DATA/$assembly/$tissue/Trinity.TEisoforms.TMM.EXPR.matrix DATA/$assembly/$tissue/Trinity.alwaysTEisoforms.TMM.EXPR.matrix DATA/$assembly/$tissue/Trinity.neverTEisoforms.TMM.EXPR.matrix > DATA/$assembly/$tissue/Trinity.all.TMM.EXPR.matrix
				sed -i "s/$/\t${assembly}/" DATA/$assembly/$tissue/Trinity.all.TMM.EXPR.matrix
				sed -i "s/$/\t${tissue}/" DATA/$assembly/$tissue/Trinity.all.TMM.EXPR.matrix

				echo -e "TISSUE ANALYZED: $tissue - ASSEMBLY: $assembly\nNUMBER OF TOTAL GENES: $nGenes - NUMBER OF TOTAL TRANSCRIPTS: $nTranscripts (avg: $avgTotal)\nNUMBER OF GENES WITH ISOFORMS WITH AND WITHOUT TE: $nGenesTE - NUMBER ISOFORMS FOR EXPRESSION TEST: $ntranscriptsTE (avg: $avgTE)\nNUMBER OF GENES WITH ISOFORMS ALWAYS TE: $nGenesAlwaysTE - NUMBER ISOFORMS FOR EXPRESSION TEST: $nTranscriptsAlwaysTE (avg: $avgAlwaysTE)\nNUMBER OF GENES WITH ISOFORMS NEVER TE: $nGenesNeverTE - NUMBER TRANSCRIPTS NEVER TE: $ntranscriptsNeverTE (avg: $avNeverTE)" > DATA/$assembly/$tissue/infoGenes.txt 
				echo -ne "Genes producing isoforms with/without TE\t$tissue\t$assembly\t$ntranscriptsTE ($nGenesTE)\nGenes producing isoforms always with TE\t$tissue\t$assembly\t$nTranscriptsAlwaysTE ($nGenesAlwaysTE)\nGenes producing isoforms without TE\t$tissue\t$assembly\t$ntranscriptsNeverTE ($nGenesNeverTE)\n"

				cat DATA/$assembly/$tissue/Trinity.all.TMM.EXPR.matrix >> DATA/all/Trinity.all.TMM.EXPR.matrix
				sed -i "s/ /\t/g" DATA/all/Trinity.all.TMM.EXPR.matrix
			done
		done

