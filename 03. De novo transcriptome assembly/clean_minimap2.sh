#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name=clean_set_transcripts
#SBATCH --output=logs/clean_set_transcripts_%a.out
#SBATCH --error=logs/clean_set_transcripts_%a.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate chimerics
module load foss/2021b
module load R/4.1.2

DIR="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/chimerics"
DIRDATA="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/DATA"

#wget http://ftp.flybase.net/releases/FB2019_06/precomputed_files/genes/fbgn_fbtr_fbpp_expanded_fb_2019_06.tsv.gz
#grep "^Dmel" fbgn_fbtr_fbpp_expanded_fb_2019_06.tsv | cut -f 3,7 | sort -u > genes_status.tab

# PREBUILT MODELS FOR DROSOPHILA
# https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/fly_Hexamer.tsv/download
# https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/Fly_logitModel.RData/download

cd ${DIR}
#SLURM_ARRAY_TASK_ID=2
strain=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/data.txt) | cut -f 1)
tissue=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/data.txt) | cut -f 2)

mkdir -p ${DIR}/clean_set_minimap2/${strain}_${tissue}
mkdir -p ${DIR}/minimap2/${strain}_${tissue}/tmp/CPAT/

> ${DIR}/clean_set_minimap2/${strain}_${tissue}/transcripts.lst
> ${DIR}/clean_set_minimap2/${strain}_${tissue}/chimerics.tmp.tab

while read transcriptInfo
	do
		transcript_trinity=$(echo "$transcriptInfo" | cut -f1 )
		lengthTranscript=$(grep -w "$transcript_trinity" ${DIR}/minimap2/${strain}_${tissue}/Trinity.hits_bin90.clean.fasta | cut -f2 -d' ' | cut -f2 -d'='  )
		gene=$(echo "$transcriptInfo" | cut -f3 )
		transcript=$(echo "$transcriptInfo" | cut -f2 )
		transcript_stringtie=$(echo "$transcriptInfo" | cut -f4 )
		class=$(echo "$transcriptInfo" | cut -f5)
		type=$(echo "$transcriptInfo" | cut -f6)

		if cat ${DIR}/trf/${strain}_${tissue}/transcripts_trf.tab | grep -vw "other_fam"  | grep -qw $transcript_trinity; then
			
			echo -e "$transcript_trinity\tTE_insertion" >> ${DIR}/clean_set_minimap2/${strain}_${tissue}/transcripts.lst

			# CPAT
			samtools faidx ${DIR}/minimap2/${strain}_${tissue}/Trinity.hits_bin90.clean.fasta $transcript_trinity > ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.fa

			cpat.py -x ${DIR}/flyPrebuild/fly_Hexamer.tsv -d ${DIR}/flyPrebuild/Fly_logitModel.RData \
			--top-orf=100 \
			--antisense \
			-g ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.fa \
			-o ${DIR}/minimap2/${strain}_${tissue}/tmp/CPAT/$transcript_trinity

			CP=$(cut -f 11 ${DIR}/minimap2/${strain}_${tissue}/tmp/CPAT/$transcript_trinity.ORF_prob.best.tsv | tail -n +2 )

			# Retrieve expression
			expr=$(grep -w "$transcript_trinity" ${DIR}/trinity/trinity_${strain}_${tissue}/salmon_outdir/quant.sf | cut -f 4)

			# Transcript status
			geneStatus=$(grep -w "$gene" ${DIRDATA}/FlyBase_r6.31/genes_status.tab | cut -f2 | sort -u | head -n1)

			if [[ $type == "analyzed" ]];then
				annotation_TEs="${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.TEs.gff"
			else
				annotation_TEs="${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.gff"
			fi

			totalExons=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.bed | wc -l)

			intersect=$(bedtools intersect -a ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.bed -b $annotation_TEs -wa -wb )

			i=1
			while read TE_insertion
				do
					TE_consensus=$(echo "$TE_insertion" | cut -f 15  | cut -f2 -d' ' | tr -d '"' | cut -f2 -d':')
					TE_order=$(grep -w "$TE_consensus" ${DIRDATA}/TEs/TE_library_family.csv | cut -f 6 )
					TE_superfamily=$(grep -w  "$TE_consensus" ${DIRDATA}/TEs/TE_library_family.csv | cut -f 7 )
					TE_family=$(grep -w  "$TE_consensus" ${DIRDATA}/TEs/TE_library_family.csv | cut -f 8 )
					TE_class=$(grep -w  "$TE_consensus" ${DIRDATA}/TEs/TE_library_family.csv | cut -f 5 )

					exon=$(echo "$TE_insertion" | cut -f1,2,3,4)
					cat ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.bed | grep "$exon" > ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.$i.bed
					posExon=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.bed | grep -n "$exon" | cut -f1 -d':')

					if [[ $totalExons -eq 1 ]]; then
						exon="Unique"
					else
						if [[ $posExon -eq 1 ]]; then
							exon="First exon"
						elif [[ $posExon -eq $totalExons ]]; then
							exon="Last exon"
						else
							exon="Middle exon"
						fi
					fi



					TE=$(echo "$TE_insertion" | cut -f7,8,9,10,11)
					cat $annotation_TEs | grep "$TE" > ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.TE.$i.gff

					#echo $i $transcript_trinity

					TE_within_exon=$(bedtools intersect -a ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.$i.bed -b ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.TE.$i.gff -wa -wb -F 1 | wc -l )
					
					exon_within_TE=$(bedtools intersect -a ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.$i.bed -b ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.TE.$i.gff -wa -wb -f 1  | wc -l )

					
					if [[ $TE_consensus == "con48_roo" ]];then

						bedtools getfasta \
						-fi ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.fa \
						-bed ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.TE.$i.gff > \
						${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.roo.$i.fa

						start_TE=$(blastn -query ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.roo.$i.fa -db ${DIR}/roo/roo.fasta -dust no -soft_masking false -word_size 7 -outfmt 6 -max_target_seqs 1 -evalue 0.05 -gapopen 5 -gapextend 2 | head -n1 | cut -f 9,10 | tr '\t' '\n' | sort | head -n1)

						if [ "$start_TE" -le 429 ] || [ "$start_TE" -ge 8660 ]; then
							roo_type=LTR
						elif [ "$start_TE" -gt 429 ] && [ "$start_TE" -lt 1275 ];then
							roo_type=repeat
						else
							roo_type=other
						fi
						
					else
						roo_type=NA
					fi

					SS=NA
					result=NA
					SS_AG=NA
					result_AG=NA
					SS_GR=NA
					result_GT=NA

					if [[ $TE_within_exon -eq 1 ]];then

						group=2

						TE_length_incorporated=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.TE.$i.gff | awk '{print $5 - $4 + 1}')
						TE_length_total=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.TE.$i.gff | awk '{print $5 - $4 + 1}')

						coordTE=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.TE.$i.gff | cut -f 4,5 | tr '\t' ':')
						coordExon=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.$i.bed | cut -f 2,3 | awk '{printf "%d %d\n", $1+1, $2}'  | tr ' ' ':')

						echo -e "$tissue\t$strain\t$transcript_trinity\t$transcript_stringtie\t$class\t$transcript\t$gene\t$lengthTranscript\t$totalExons\t$posExon\t$exon\t$coordExon\t$coordTE\tTE_within_exon\t$TE_consensus\t$TE_family\t$TE_superfamily\t$TE_order\t$TE_class\t$TE_length_incorporated\t$TE_length_total\t$CP\t$geneStatus\t$expr\t$SS\t$result\t$group\t$roo_type" >> ${DIR}/clean_set_minimap2/${strain}_${tissue}/chimerics.tmp.tab

					elif [[ $exon_within_TE -eq 1 ]];then

						group=1

						coordTE=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.TE.$i.gff | cut -f 4,5 | tr '\t' ':')
						coordExon=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.$i.bed | cut -f 2,3 | awk '{printf "%d %d\n", $1+1, $2}'  | tr ' ' ':')

						start_TE=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.TE.$i.gff | cut -f 4)
						end_TE=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.TE.$i.gff | cut -f 5)
						start_exon=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.$i.bed | cut -f2 )
						start_exon=$((start_exon + 1))
						end_exon=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.$i.bed | cut -f3 )


						TE_length_incorporated=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.$i.bed | awk '{print $3 - $2 }')
						TE_length_total=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.TE.$i.gff | awk '{print $5 - $4 + 1}')

						#echo "$transcript_trinity exon_within_TE $exon $TE_length_incorporated $TE_length_total"

						if [[ $exon == "Last exon" ]]; then
							SSstart=$(($start_exon-11))
							SSend=$(($start_exon+1))

							echo -e "$transcript_trinity\t$SSstart\t$SSend\t$transcript\t.\t+" > ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.bed
							sed "1s/.*/>$transcript_trinity/" ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.fa > ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.name.fa
							bedtools getfasta -fi ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.name.fa -bed ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.bed -name > ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.fasta 
							#mkdir totalOverlapTEexon/SS/fimo_${transcript}
							conda activate meme
							fimo --thresh 0.05 --verbosity 1 --oc ${DIR}/minimap2/${strain}_${tissue}/tmp/ ${DIR}/SSmotifs/AGmotif.txt ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.fasta 
							conda activate chimerics
							fimoResult=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/fimo.gff | wc -l )
							
							if [[ $fimoResult -gt 1 ]]; then
								result=$(tail -n +2 ${DIR}/minimap2/${strain}_${tissue}/tmp/fimo.gff | cut -f 4,6 -d';' | head -n1) 
								SS="Canonical AG"
								#echo "$transcript_trinity $exon exon_within_TE_5' $TE_length_incorporated $TE_length_total $SS"
							else
								result="No enrichment"
								sequence=$(tail -n +2 ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.fasta)
								SS="No real SS/Non-canonical AG ($sequence)"
								#echo "$transcript_trinity $exon exon_within_TE_5' $TE_length_incorporated $TE_length_total $SS"
							fi

							echo -e "$tissue\t$strain\t$transcript_trinity\t$transcript_stringtie\t$class\t$transcript\t$gene\t$lengthTranscript\t$totalExons\t$posExon\t$exon\t$coordExon\t$coordTE\texon_within_TE_5\t$TE_consensus\t$TE_family\t$TE_superfamily\t$TE_order\t$TE_class\t$TE_length_incorporated\t$TE_length_total\t$CP\t$geneStatus\t$expr\t$SS\t$result\t$group\t$roo_type" >> ${DIR}/clean_set_minimap2/${strain}_${tissue}/chimerics.tmp.tab
						
						elif [[ $exon == "First exon" ]]; then
							SSstart=$(($end_exon-3))
							SSend=$(($end_exon+6))

							echo -e "$transcript_trinity\t$SSstart\t$SSend\t$transcript\t.\t+" > ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.bed
							sed "1s/.*/>$transcript_trinity/" ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.fa > ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.name.fa
							bedtools getfasta -fi ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.name.fa -bed ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.bed -name > ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.fasta 
							#mkdir totalOverlapTEexon/SS/fimo_${transcript}
							conda activate meme
							fimo --thresh 0.05 --verbosity 1 --oc ${DIR}/minimap2/${strain}_${tissue}/tmp/ ${DIR}/SSmotifs/GTmotif.txt ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.fasta 
							conda activate chimerics
							fimoResult=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/fimo.gff | wc -l )
							
							if [[ $fimoResult -gt 1 ]]; then
								result=$(tail -n +2 ${DIR}/minimap2/${strain}_${tissue}/tmp/fimo.gff | cut -f 4,6 -d';' | head -n1) 
								SS="Canonical GT"
								#echo "$transcript_trinity $exon exon_within_TE_3' $TE_length_incorporated $TE_length_total $SS"

							else
								result="No enrichment"
								sequence=$(tail -n +2 ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.fasta)
								SS="No real SS/Non-canonical GT ($sequence)"
								#echo "$transcript_trinity $exon exon_within_TE_3' $TE_length_incorporated $TE_length_total $SS"
							fi 
							
							echo -e "$tissue\t$strain\t$transcript_trinity\t$transcript_stringtie\t$class\t$transcript\t$gene\t$lengthTranscript\t$totalExons\t$posExon\t$exon\t$coordExon\t$coordTE\texon_within_TE_3\t$TE_consensus\t$TE_family\t$TE_superfamily\t$TE_order\t$TE_class\t$TE_length_incorporated\t$TE_length_total\t$CP\t$geneStatus\t$expr\t$SS\t$result\t$group\t$roo_type" >> ${DIR}/clean_set_minimap2/${strain}_${tissue}/chimerics.tmp.tab
						else
							
							SSstart=$(($start_exon-11))
							SSend=$(($start_exon+1))

							echo -e "$transcript_trinity\t$SSstart\t$SSend\t$transcript\t.\t+" > ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.bed
							sed "1s/.*/>$transcript_trinity/" ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.fa > ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.name.fa
							bedtools getfasta -fi ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.name.fa -bed ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.bed -name > ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.fasta 
							#mkdir totalOverlapTEexon/SS/fimo_${transcript}
							conda activate meme
							fimo --thresh 0.05 --verbosity 1 --oc ${DIR}/minimap2/${strain}_${tissue}/tmp/ ${DIR}/SSmotifs/AGmotif.txt ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.fasta 
							conda activate chimerics
							fimoResult=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/fimo.gff | wc -l )
							
							if [[ $fimoResult -gt 1 ]]; then
								result_AG=$(tail -n +2 ${DIR}/minimap2/${strain}_${tissue}/tmp/fimo.gff | cut -f 4,6 -d';' | head -n1) 
								SS_AG="Canonical AG"
								#echo "$transcript_trinity $exon exon_within_TE_5' $TE_length_incorporated $TE_length_total $SS"
							else
								result_AG="No enrichment"
								sequence_AG=$(tail -n +2 ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.fasta)
								SS_AG="No real SS/Non-canonical AG ($sequence_AG)"
								#echo "$transcript_trinity $exon exon_within_TE_5' $TE_length_incorporated $TE_length_total $SS"
							fi

							SSstart=$(($end_exon-3))
							SSend=$(($end_exon+6))

							echo -e "$transcript_trinity\t$SSstart\t$SSend\t$transcript\t.\t+" > ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.bed
							sed "1s/.*/>$transcript_trinity/" ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.fa > ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.name.fa
							bedtools getfasta -fi ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.name.fa -bed ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.bed -name > ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.fasta 
							#mkdir totalOverlapTEexon/SS/fimo_${transcript}
							conda activate meme
							fimo --thresh 0.05 --verbosity 1 --oc ${DIR}/minimap2/${strain}_${tissue}/tmp/ ${DIR}/SSmotifs/GTmotif.txt ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.fasta 
							conda activate chimerics
							fimoResult=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/fimo.gff | wc -l )
							
							if [[ $fimoResult -gt 1 ]]; then
								result_GT=$(tail -n +2 ${DIR}/minimap2/${strain}_${tissue}/tmp/fimo.gff | cut -f 4,6 -d';' | head -n1) 
								SS_GT="Canonical GT"
								#echo "$transcript_trinity $exon exon_within_TE_3' $TE_length_incorporated $TE_length_total $SS"

							else
								result_GT="No enrichment"
								sequence_GT=$(tail -n +2 ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.fasta)
								SS_GT="No real SS/Non-canonical GT ($sequence_GT)"
								#echo "$transcript_trinity $exon exon_within_TE_3' $TE_length_incorporated $TE_length_total $SS"
							fi 
							
							echo -e "$tissue\t$strain\t$transcript_trinity\t$transcript_stringtie\t$class\t$transcript\t$gene\t$lengthTranscript\t$totalExons\t$posExon\t$exon\t$coordExon\t$coordTE\texon_within_TE\t$TE_consensus\t$TE_family\t$TE_superfamily\t$TE_order\t$TE_class\t$TE_length_incorporated\t$TE_length_total\t$CP\t$geneStatus\t$expr\t$SS_AG;$SS_GT\t$result_AG;$result_GT\t$group\t$roo_type" >> ${DIR}/clean_set_minimap2/${strain}_${tissue}/chimerics.tmp.tab
						fi

					else
						group=1

						coordTE=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.TE.$i.gff | cut -f 4,5 | tr '\t' ':')
						coordExon=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.$i.bed | cut -f 2,3 | awk '{printf "%d %d\n", $1+1, $2}'  | tr ' ' ':')

						start_TE=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.TE.$i.gff | cut -f 4)
						end_TE=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.TE.$i.gff | cut -f 5)
						start_exon=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.$i.bed | cut -f2 )
						start_exon=$((start_exon + 1))
						end_exon=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.$i.bed | cut -f3 )

						if [ $start_TE -lt $start_exon -a $end_TE -lt $end_exon ]; then
							TE_length_total=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.TE.$i.gff | awk '{print $5 - $4 + 1}')
							TE_length_incorporated=$(( $end_TE - $start_exon + 1 ))

							if [ "$exon" != "First exon" -a $totalExons -gt 1  ]; then
								#echo "TRINITY_DN543_c0_g1_i3"
								SSstart=$(($start_exon-11))
								SSend=$(($start_exon+1))

								echo -e "$transcript_trinity\t$SSstart\t$SSend\t$transcript\t.\t+" > ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.bed
								sed "1s/.*/>$transcript_trinity/" ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.fa > ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.name.fa
								bedtools getfasta -fi ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.name.fa -bed ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.bed -name > ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.fasta 
								#mkdir totalOverlapTEexon/SS/fimo_${transcript}
								conda activate meme
								fimo --thresh 0.05 --verbosity 1 --oc ${DIR}/minimap2/${strain}_${tissue}/tmp/ ${DIR}/SSmotifs/AGmotif.txt ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.fasta 
								conda activate chimerics
								fimoResult=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/fimo.gff | wc -l )
								
								if [[ $fimoResult -gt 1 ]]; then
									result=$(tail -n +2 ${DIR}/minimap2/${strain}_${tissue}/tmp/fimo.gff | cut -f 4,6 -d';' | head -n1) 
									SS="Canonical AG"
									#echo "$transcript_trinity $exon overlap_5' $TE_length_incorporated $TE_length_total $SS"
								else
									result="No enrichment"
									sequence=$(tail -n +2 ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.fasta)
									SS="No real SS/Non-canonical AG ($sequence)"
									#echo "$transcript_trinity $exon overlap_5' $TE_length_incorporated $TE_length_total $SS"
								fi 
								
								echo -e "$tissue\t$strain\t$transcript_trinity\t$transcript_stringtie\t$class\t$transcript\t$gene\t$lengthTranscript\t$totalExons\t$posExon\t$exon\t$coordExon\t$coordTE\tTE_overlap_5\t$TE_consensus\t$TE_family\t$TE_superfamily\t$TE_order\t$TE_class\t$TE_length_incorporated\t$TE_length_total\t$CP\t$geneStatus\t$expr\t$SS\t$result\t$group\t$roo_type" >> ${DIR}/clean_set_minimap2/${strain}_${tissue}/chimerics.tmp.tab
							
							else
								SS=NA
								result=NA
								echo -e "$tissue\t$strain\t$transcript_trinity\t$transcript_stringtie\t$class\t$transcript\t$gene\t$lengthTranscript\t$totalExons\t$posExon\t$exon\t$coordExon\t$coordTE\tTE_overlap_5\t$TE_consensus\t$TE_family\t$TE_superfamily\t$TE_order\t$TE_class\t$TE_length_incorporated\t$TE_length_total\t$CP\t$geneStatus\t$expr\t$SS\t$result\t$group\t$roo_type" >> ${DIR}/clean_set_minimap2/${strain}_${tissue}/chimerics.tmp.tab

							fi

						elif [ $start_TE -gt $start_exon -a $end_TE -gt $end_exon ]; then
							TE_length_total=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.TE.$i.gff | awk '{print $5 - $4 + 1}')
							TE_length_incorporated=$(( $end_exon - $start_TE + 1 ))

							if [ "$exon" != "Last exon" -a $totalExons -gt 1 ]; then
								#echo "TRINITY_DN4642_c0_g1_i1"
								SSstart=$(($end_exon-3))
								SSend=$(($end_exon+6))

								echo -e "$transcript_trinity\t$SSstart\t$SSend\t$transcript\t.\t+" > ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.bed
								sed "1s/.*/>$transcript_trinity/" ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.fa > ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.name.fa
								bedtools getfasta -fi ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.name.fa -bed ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.bed -name > ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.fasta 
								#mkdir totalOverlapTEexon/SS/fimo_${transcript}
								conda activate meme
								fimo --thresh 0.05 --verbosity 1 --oc ${DIR}/minimap2/${strain}_${tissue}/tmp/ ${DIR}/SSmotifs/GTmotif.txt ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.fasta 
								conda activate chimerics
								fimoResult=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/fimo.gff | wc -l )
								if [[ $fimoResult -gt 1 ]]; then
									result=$(tail -n +2 ${DIR}/minimap2/${strain}_${tissue}/tmp/fimo.gff | cut -f 4,6 -d';' | head -n1) 
									SS="Canonical GT"
									#echo "$transcript_trinity $exon overlap_3' $TE_length_incorporated $TE_length_total $SS"
								else
									result="No enrichment"
									sequence=$(tail -n +2 ${DIR}/minimap2/${strain}_${tissue}/tmp/${transcript_trinity}.exon_$posExon.SS.fasta)
									SS="No real SS/Non-canonical GT ($sequence)"
									#echo "$transcript_trinity $exon overlap_3' $TE_length_incorporated $TE_length_total $SS"
								fi 
									echo -e "$tissue\t$strain\t$transcript_trinity\t$transcript_stringtie\t$class\t$transcript\t$gene\t$lengthTranscript\t$totalExons\t$posExon\t$exon\t$coordExon\t$coordTE\tTE_overlap_3\t$TE_consensus\t$TE_family\t$TE_superfamily\t$TE_order\t$TE_class\t$TE_length_incorporated\t$TE_length_total\t$CP\t$geneStatus\t$expr\t$SS\t$result\t$group\t$roo_type" >> ${DIR}/clean_set_minimap2/${strain}_${tissue}/chimerics.tmp.tab
							else
								SS=NA
								result=NA
								echo -e "$tissue\t$strain\t$transcript_trinity\t$transcript_stringtie\t$class\t$transcript\t$gene\t$lengthTranscript\t$totalExons\t$posExon\t$exon\t$coordExon\t$coordTE\tTE_overlap_3\t$TE_consensus\t$TE_family\t$TE_superfamily\t$TE_order\t$TE_class\t$TE_length_incorporated\t$TE_length_total\t$CP\t$geneStatus\t$expr\t$SS\t$result\t$group\t$roo_type" >> ${DIR}/clean_set_minimap2/${strain}_${tissue}/chimerics.tmp.tab
							fi
						else
							echo "$transcript_trinity PROBLEM"
						fi
					fi

					i=$(( $i + 1 ))
				done <<< "$intersect"

		else
			echo -e "$transcript_trinity\tSSR" >> ${DIR}/clean_set_minimap2/${strain}_${tissue}/transcripts.lst
		fi

	done < <(grep analyzed ${DIR}/minimap2/${strain}_${tissue}/transcripts_status.lst)

