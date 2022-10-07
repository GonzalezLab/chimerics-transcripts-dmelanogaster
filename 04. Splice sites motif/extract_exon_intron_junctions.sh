#!/bin/bash

DIR=/homes/users/mcoronado/scratch/5GenomesProject/

cut -f9 ${DIR}/DATA/FlyBase_r6.31/dmel-chr-r6.31.gtf | cut -f6 -d' ' | tr -d '"' | tr -d ";" |sort -u > list_of_transcripts.txt


for i in `less list_of_transcripts.txt`
do
	nExon=$(grep -e "transcript_id \"$i\";" ${DIR}/DATA/FlyBase_r6.31/dmel-chr-r6.31.gtf | grep CDS |wc -l)
	if [ "$nExon" -gt 3 ]
		then grep -e "transcript_id \"$i\";" ${DIR}/DATA/FlyBase_r6.31/dmel-chr-r6.31.gtf | grep CDS |  sed '1d;$d'
fi
done  > all_exons.gtf

cat all_exons.gtf | cut -f1,2 -d' ' | sort -u > all_exons_unique.gtf

cat all_exons_unique.gtf | shuf -n 500 > dmel_500_random_exons3.gtf

echo > AGsites.bed
echo > GTsites.bed

while IFS= read -r exonLine
do

	chr=$(echo $exonLine | cut -f1 -d' ' )
	start=$(echo $exonLine | cut -f4 -d' ' )
	end=$(echo $exonLine | cut -f5 -d' ' )
	strand=$(echo $exonLine | cut -f7 -d' ' )
	info=$(echo $exonLine | cut -f10 -d' ' | tr -d '"' |tr -d ';')

	if [ "$strand" == "+" ]
	then
		AGstartPlus=$(($start-11))
		AGendPlus=$(($start+1))
		GTstartPlus=$(($end-4))
		GTendPlus=$(($end+7))
		echo -e "$chr\t$AGstartPlus\t$AGendPlus\t$info\t.\t$strand" >> AGsites.bed
		echo -e "$chr\t$GTstartPlus\t$GTendPlus\t$info\t.\t$strand" >> GTsites.bed
	fi

	if [ "$strand" == "-" ]
	then
		GTstartPlus=$(($start-8))
		GTendPlus=$(($start+3))
		AGstartPlus=$(($end-2))
		AGendPlus=$(($end+10))
		echo -e "$chr\t$AGstartPlus\t$AGendPlus\t$info\t.\t$strand" >> AGsites.bed
		echo -e "$chr\t$GTstartPlus\t$GTendPlus\t$info\t.\t$strand" >> GTsites.bed
	fi

done < dmel_500_random_exons3.gtf

bedtools getfasta -fi ${DIR}/DATA/FlyBase_r6.31/dmel-chr.fasta -bed AGsites.bed -s -name > AGsites.fasta 
bedtools getfasta -fi ${DIR}/DATA/FlyBase_r6.31/dmel-chr.fasta -bed GTsites.bed -s -name > GTsites.fasta 