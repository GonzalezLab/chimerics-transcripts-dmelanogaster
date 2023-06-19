#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --job-name=trf
#SBATCH --output=logs/trf_%a.out
#SBATCH --error=logs/trf_%a.err


source ~/miniconda3/etc/profile.d/conda.sh
conda activate chimerics

DIR="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/chimerics"

cd ${DIR}

strain=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/data.txt) | cut -f 1)
tissue=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/data.txt) | cut -f 2)

mkdir -p ${DIR}/trf/${strain}_${tissue}/tmp

> ${DIR}/trf/${strain}_${tissue}/percentageSRR.tab

while read transcriptID
	do
		awk -v transcriptID="$transcriptID" ' $1 == transcriptID ' ${DIR}/RepeatMasker/${strain}_${tissue}/Trinity.hits_bin90.fasta.out.gff > ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.tmp.gff

		cat ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.tmp.gff | cut -f 1,4,5,9 | awk '{print $1"\t"$2"\t"$3"\t"$5}' FS="[\"]" OFS="\t" | cut -f 1,2,3,5 > ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.gff
		rm -rf ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.collapse.bed
		> ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.collapse.gff 
		cat ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.gff | cut -f 4 | sort | uniq | while read C
			do
				awk -v C=${C} '($4==C)' ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.gff | sort -t $'\t' -k1,1 -k2,2n | bedtools merge  -c 4 -o distinct >> ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.collapse.gff 
			done
		sed -i "s/Motif://g" ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.collapse.gff	
		bedtools getfasta -fi ${DIR}/RepeatMasker/${strain}_${tissue}/Trinity.hits_bin90.fasta.masked -bed ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.collapse.gff -name -fo ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.collapse.fa

		cd ${DIR}/trf/${strain}_${tissue}/tmp/

		${DIR}/trf409.linux64  ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.collapse.fa 2 3 5 80 10 20 15 -h -d

		i=0
		while read line ; do
			if echo "$line" | grep -q "Sequence"; then
				i=$(( $i + 1 ))
				echo $line > ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.trf.$i
				continue
			fi
			echo "$line" >> ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.trf.$i
			#echo $i
		done < <(cat ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.collapse.fa.2.3.5.80.10.20.15.dat)

		> ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.collapse.fa.trf		
		for seq in `seq 1 $i`
			do
				merge=$(cat ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.trf.$seq | grep "^[0-9]" | cut -f 1,2 -d' ' | sed  -e 's/^/TE /'  | tr " " "\t" | mergeBed | cut -f2,3 | tr "\t" " ")
				info=$(cat ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.trf.$seq | grep -E "^Sequence|^Parameter")
				echo -e "$info\n$merge"
			done >> ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.collapse.fa.trf	

		python ${DIR}/trf_search.py ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.collapse.fa $transcriptID > ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.percentageSRR.tab 

		cat ${DIR}/trf/${strain}_${tissue}/tmp/$transcriptID.percentageSRR.tab >> ${DIR}/trf/${strain}_${tissue}/percentageSRR.tab
		cd ${DIR}

	done < <(cut -f1 ${DIR}/RepeatMasker/${strain}_${tissue}/Trinity.hits_bin90.fasta.out.gff | grep -v "^#" | sort -u )