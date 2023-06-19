#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=minimap2
#SBATCH --output=logs/minimap2_%a.out
#SBATCH --error=logs/minimap2_%a.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate chimerics
module load foss/2021b
module load R/4.1.2

DIR="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/chimerics"
DIRDATA="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/DATA"

cd ${DIR}

strain=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/data.txt) | cut -f 1)
tissue=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/data.txt) | cut -f 2)

mkdir -p ${DIR}/minimap2/${strain}_${tissue}

# Extract chimerics that are not fussioned
awk ' $3 == 1 ' ${DIR}/check_fussion/output/${strain}/${tissue}/resultDistance.tab | cut -f1 | sort -u > ${DIR}/check_fussion/output/${strain}/${tissue}/genes_no_fussioned.lst

seqtk subseq ${DIR}/trinity/trinity_${strain}_${tissue}/Trinity.hits_bin90.fasta \
${DIR}/check_fussion/output/${strain}/${tissue}/genes_no_fussioned.lst > \
${DIR}/minimap2/${strain}_${tissue}/Trinity.hits_bin90.clean.fasta 

grep -w -f ${DIR}/check_fussion/output/${strain}/${tissue}/genes_no_fussioned.lst ${DIR}/trinity/trinity_${strain}_${tissue}/Trinity_refTrans.blastn.hist_bin90.list.names.keep.tab > ${DIR}/minimap2/${strain}_${tissue}/Trinity_refTrans.blastn.hist_bin90.chimerics.keep.tab

rm -rf ${DIR}/minimap2/${strain}_${tissue}/tmp
mkdir -p ${DIR}/minimap2/${strain}_${tissue}/tmp

> ${DIR}/minimap2/${strain}_${tissue}/transcripts_status.lst	
mkdir ${DIR}/minimap2/${strain}_${tissue}/genome_plot/

while read transcriptInfo
	do
		transcript_trinity=$(echo "$transcriptInfo" | cut -f1 )
		gene=$(echo "$transcriptInfo" | cut -f4 )
		transcript=$(echo "$transcriptInfo" | cut -f5 )
		transcript_stringtie=$(echo "$transcriptInfo" | cut -f2 )
		class=$(echo "$transcriptInfo" | cut -f3 )
		TEs=$(grep -w "$transcript_trinity" ${DIR}/RepeatMasker/${strain}_${tissue}/Trinity.hits_bin90.fasta.out.gff | cut -f 9 | cut -f2 -d' ' | tr -d '"' | cut -f2 -d':' | sort -u )

		grep -P "\tgene\t" /lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/GeneTransference/${strain}.gff | grep $gene | awk '{OFS="\t"}{$4=$4-2000; $5=$5+2000; print}' > ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.gff

		if ! [ -s ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.gff ]; then
			echo -e "$transcript_trinity\t$transcript\t$gene\t$transcript_stringtie\t$class\tgene_not_transferred" >> ${DIR}/minimap2/${strain}_${tissue}/transcripts_status.lst
			continue
		fi
		
		bedtools getfasta -s -fi ${DIRDATA}/genomeAssembly/${strain}/${strain}.fasta -bed ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.gff >  ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.fa

		samtools faidx ${DIR}/minimap2/${strain}_${tissue}/Trinity.hits_bin90.clean.fasta  $transcript_trinity  > ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.fa

		minimap2 -ax splice --secondary=no --sam-hit-only -C5 -t4 ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.fa \
		${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.fa > ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.sam

		samtools view -h -Sq 1 ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.sam > \
		${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.no_mapQ0.sam

		nAlignemnts=$(samtools view -c ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.no_mapQ0.sam)
		
		if [[ $nAlignemnts -gt 1 ]];then
		
			len=$(samtools view -h ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.sam | awk '{print length($10)'} | sort -hr | head -n1)
			samtools view -h -Sq 1 ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.sam  | awk -v len="$read_length" 'length($10)==len || NR==1 {print $0}' > ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.no_mapQ0.sam

		elif [[ $nAlignemnts -eq 1 ]];then
			samtools view -h -Sq 1 ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.sam > \
			${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.no_mapQ0.sam	
		fi
		
		bedtools bamtobed -split -i ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.no_mapQ0.sam > \
		${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.bed

		nLines=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.bed | wc -l )

		if [[ $nLines -eq 1 ]]; then
			nExons=$(grep -w $transcript /lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/GeneTransference/${strain}.gff | grep -P "\texon\t" | wc -l )
			if [[ $nExons -eq 1 ]]; then
				single_unit="ok"
			else
				single_unit="not_ok"
			fi
		elif [[ $nLines -gt 1 ]]; then
			single_unit="ok"
		fi

		if [ -s ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.bed -a $single_unit == "ok" ]; then

			python3 ~/bin/bed_to_gff3/bed_to_gff3_converter.py ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.bed \
			${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.tmp.gff 

			sed -i 's/:[0-9]*-[0-9]*(.)\t/\t/g'  ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.tmp.gff 
		
			sed "s/$transcript_trinity/CDS/" ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.tmp.gff | sed 's/bed_to_gff3_converter/chimerics/g' | sed '3,$ s/-/+/g' > ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.gff

			awk 'BEGIN{OFS="\t"} NR<=3{next} {if(NR==4){min=$4;max=$5}; if($4<min){min=$4}; if($5>max){max=$5}; print} END{print $1, "chimerics", "gene", min, max, $6, $7, $8, $9}' ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.gff > ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.TE.gff
			
			sed -i "s/$/; name=$transcript/g" ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.TE.gff

			RepeatMasker -pa 2 -e rmblast -lib ${DIRDATA}/TEs/consensuses_curated_v4.fasta \
			-dir ${DIR}/minimap2/${strain}_${tissue}/tmp/ \
			-norna -nolow -s -cutoff 250 -xsmall -no_is -gff  ${DIR}/minimap2/${strain}_${tissue}/tmp/$gene.$transcript_trinity.fa 

			#echo -e "$TEs" > ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.TEs.lst

			#grep -w -f ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.TEs.lst ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.gff > ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.TEs.gff 

			if [ -s ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.gff  ]; then

				echo -e "$TEs" > ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.TEs.lst

				grep -w -f ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.TEs.lst ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.gff > ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.TEs.gff 

				bedtools intersect -a ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.bed \
				-b ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.TEs.gff -wa -wb > ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.intersect.TE.bed

				if [ -s ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.intersect.TE.bed ];then

					i=1
					while read insertions
						do
							chr=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.gff | tail -n +4  | cut -f1 | sort -u)
							start=$(echo "$insertions" | cut -f4)
							end=$(echo "$insertions" | cut -f5)
							family=$(echo "$insertions" | cut -f 9 | cut -f2 -d' ' | tr -d '"' | cut -f2 -d':')
							echo -e "$chr\tchimeric\tgene\t$start\t$end\t.\t+\t.\tID=$family.$i; name=$family" >> ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.TE.gff
							echo -e "$chr\tchimeric\tCDS\t$start\t$end\t.\t+\t.\tID=$family.$i; name=$family" >> ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.TE.gff
							echo -e "$chr\tchimeric\trepeat_region\t$start\t$end\t.\t+\t.\tID=$family.$i; name=$family" >> ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.TE.gff
							i=$(( $i + 1 ))
						done < <(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.TEs.gff )

					Rscript ${DIR}/visualize_chimeric.R "${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.TE.gff" "$transcript_trinity" "$transcript" "${DIR}/minimap2/${strain}_${tissue}/genome_plot/"
					echo -e "$transcript_trinity\t$transcript\t$gene\t$transcript_stringtie\t$class\tanalyzed" >> ${DIR}/minimap2/${strain}_${tissue}/transcripts_status.lst

				else
					bedtools intersect -a ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.bed \
					-b ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.gff -wa -wb > ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.intersect.TE.bed

					if [ -s ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.intersect.TE.bed ];then

						orders_expected=$(grep -f ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.TEs.lst ${DIRDATA}/TEs/TE_library_family.csv | cut -f 7)
						TE_all_intersect=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.intersect.TE.bed | cut -f 15 | cut -f2 -d' ' | tr -d '"' | cut -f2 -d':' | sort -u) 
						TE_found=0
						while read TE
							do
								order=$(grep -w "$TE" ${DIRDATA}/TEs/TE_library_family.csv | cut -f 7 ) 
								if echo "$orders_expected" | grep -q $order; then
									TE_found=1
								fi
							done <<< "$TE_all_intersect"

						if [[ $TE_found -eq 1 ]];then
							
							i=1
							while read insertions
								do
									chr=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.gff | tail -n +4  | cut -f1 | sort -u)
									start=$(echo "$insertions" | cut -f4)
									end=$(echo "$insertions" | cut -f5)
									family=$(echo "$insertions" | cut -f 9 | cut -f2 -d' ' | tr -d '"' | cut -f2 -d':')
									echo -e "$chr\tchimeric\tgene\t$start\t$end\t.\t+\t.\tID=$family.$i; name=$family" >> ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.TE.gff
									echo -e "$chr\tchimeric\tCDS\t$start\t$end\t.\t+\t.\tID=$family.$i; name=$family" >> ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.TE.gff
									echo -e "$chr\tchimeric\trepeat_region\t$start\t$end\t.\t+\t.\tID=$family.$i; name=$family" >> ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.TE.gff
									i=$(( $i + 1 ))
								done < <(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.gff | tail -n +4)

							Rscript ${DIR}/visualize_chimeric.R "${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.TE.gff" "$transcript_trinity" "$transcript" "${DIR}/minimap2/${strain}_${tissue}/genome_plot/"

							echo -e "$transcript_trinity\t$transcript\t$gene\t$transcript_stringtie\t$class\tanalyzed_same_family" >> ${DIR}/minimap2/${strain}_${tissue}/transcripts_status.lst
						else
							
							i=1
							while read insertions
								do
									chr=$(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.gff | tail -n +4  | cut -f1 | sort -u)
									start=$(echo "$insertions" | cut -f4)
									end=$(echo "$insertions" | cut -f5)
									family=$(echo "$insertions" | cut -f 9 | cut -f2 -d' ' | tr -d '"' | cut -f2 -d':')
									echo -e "$chr\tchimeric\tgene\t$start\t$end\t.\t+\t.\tID=$family.$i; name=$family" >> ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.TE.gff
									echo -e "$chr\tchimeric\tCDS\t$start\t$end\t.\t+\t.\tID=$family.$i; name=$family" >> ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.TE.gff
									echo -e "$chr\tchimeric\trepeat_region\t$start\t$end\t.\t+\t.\tID=$family.$i; name=$family" >> ${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.TE.gff
									i=$(( $i + 1 ))
								done < <(cat ${DIR}/minimap2/${strain}_${tissue}/tmp/${gene}.$transcript_trinity.fa.out.gff | tail -n +4 )

							Rscript ${DIR}/visualize_chimeric.R "${DIR}/minimap2/${strain}_${tissue}/tmp/$transcript.$transcript_trinity.TE.gff" "$transcript_trinity" "$transcript" "${DIR}/minimap2/${strain}_${tissue}/genome_plot/"

							echo -e "$transcript_trinity\t$transcript\t$gene\t$transcript_stringtie\t$class\tanalyzed_other_family" >> ${DIR}/minimap2/${strain}_${tissue}/transcripts_status.lst
						fi
					else
						echo -e "$transcript_trinity\t$transcript\t$gene\t$transcript_stringtie\t$class\tno_intersect_TE" >> ${DIR}/minimap2/${strain}_${tissue}/transcripts_status.lst
					fi
				fi
			else
				echo -e "$transcript_trinity\t$transcript\t$gene\t$transcript_stringtie\t$class\tno_TE" >> ${DIR}/minimap2/${strain}_${tissue}/transcripts_status.lst	
			fi

		else
			if [[ $nLines -eq 1 ]]; then
				echo -e "$transcript_trinity\t$transcript\t$gene\t$transcript_stringtie\t$class\tsingle_transcript_unit" >> ${DIR}/minimap2/${strain}_${tissue}/transcripts_status.lst
			else
				echo -e "$transcript_trinity\t$transcript\t$gene\t$transcript_stringtie\t$class\tmap_failed" >> ${DIR}/minimap2/${strain}_${tissue}/transcripts_status.lst
			fi
		fi

	done < <(cat ${DIR}/minimap2/${strain}_${tissue}/Trinity_refTrans.blastn.hist_bin90.chimerics.keep.tab  )


