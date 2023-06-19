#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name=fussion_transcripts
#SBATCH --output=logs/identity_%a.out
#SBATCH --error=logs/identity_rm_%a.err


source ~/miniconda3/etc/profile.d/conda.sh
conda activate chimerics

DIR="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/chimerics"
DIRDATA="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/DATA"

cd ${DIR}

strain=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/data.txt) | cut -f 1)
tissue=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/data.txt) | cut -f 2)

mkdir -p ${DIR}/check_fussion/input/$strain/$tissue/
mkdir -p ${DIR}/check_fussion/output/$strain/$tissue/
mkdir -p ${DIR}/check_fussion/tmp/$strain/$tissue/

grep -v "^#" ${DIR}/RepeatMasker/${strain}_${tissue}/Trinity.hits_bin90.fasta.out.gff | cut -f1 | sort -u > ${DIR}/check_fussion/input/$strain/$tissue/chimerics.lst

seqtk subseq ${DIR}/RepeatMasker/${strain}_${tissue}/Trinity.hits_bin90.fasta.masked ${DIR}/check_fussion/input/$strain/$tissue/chimerics.lst > ${DIR}/check_fussion/input/$strain/$tissue/chimerics.fasta

blastn -query ${DIR}/check_fussion/input/$strain/$tissue/chimerics.fasta \
       -db ${DIR}/DataAssemblyMergedReference/stringtie_merged.fasta \
       -out ${DIR}/check_fussion/output/$strain/$tissue/chimerics.blastn \
       -num_threads 4 \
       -evalue 1e-10 \
       -outfmt "6 qseqid qlen sseqid slen qcovs pident length mismatch gapopen qstart qend sstart send evalue bitscore"

Rscript calculate_percentage.R ${DIR}/check_fussion/output/$strain/$tissue/chimerics.blastn $tissue $strain

cat ${DIR}/check_fussion/output/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.tmp.blastn | awk -F'\t' '{ if (!seen[$1 "-" $4]) { print $0; seen[$1 "-" $4]=1 } }'  >  ${DIR}/check_fussion/output/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.blastn 

#awk ' $7 > 100  ' ${DIR}/check_fussion/output/$strain/$tissue/chimerics.blastn | tr ':' '\t' | awk -F'\t' '{ if (!seen[$1 "-" $4]) { print $0; seen[$1 "-" $4]=1 } }'  >  ${DIR}/check_fussion/output/$strain/$tissue/chimerics.filterLengthMin100.blastn


while read transcript
       do
              transcriptID=$(awk -v transcript="$transcript" ' $1 == transcript '  ${DIR}/trinity/trinity_${strain}_${tissue}/Trinity_refTrans.blastn.hist_bin90.list.names.tab | cut -f 2 | sort -u)
              nMatches=$(awk -v transcript="$transcript" ' $1 == transcript '  ${DIR}/check_fussion/output/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.blastn | cut -f4 | sort -u | wc -l )

              if [[ $nMatches -gt 1 ]]; then

                     awk -v transcript="$transcript" ' $1 == transcript '  ${DIR}/check_fussion/output/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.blastn | awk -v OFS='\t' '{print $1,$11,$12,$4}' > ${DIR}/check_fussion/tmp/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.$transcript.blastn

                     bedmap --count --echo-map-range --echo-map-id --fraction-both 0.8 --delim '\t'  ${DIR}/check_fussion/tmp/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.$transcript.blastn | cut -f2- | sort -u > ${DIR}/check_fussion/tmp/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.$transcript.merge.blastn

                     n=$(cat ${DIR}/check_fussion/tmp/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.$transcript.merge.blastn | wc -l)

                     if [[ $n -eq 1 ]];then
                            gene=$(cut -f 4 ${DIR}/check_fussion/tmp/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.$transcript.merge.blastn)
                            echo -e "$transcript\t$transcriptID\t$n\t$gene\tND\tduplicated_gene"
                     
                     elif [[ $n -lt $nMatches ]];then
                            #echo -e "$transcript\t$transcriptID\t$n\tND\tmaspeque"
                            
                            awk -v transcript="$transcript" ' $1 == transcript '  ${DIR}/check_fussion/output/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.blastn | cut -f4 | sort -u | sed 's/^/ID=/' > ${DIR}/check_fussion/output/$strain/$tissue/gene.lst
                            
                            grep -f ${DIR}/check_fussion/output/$strain/$tissue/gene.lst dmel-all-no-analysis-r6.31.gene.gff >  ${DIR}/check_fussion/output/$strain/$tissue/gene.gff
                            
                            sed -i 's/;.*//' ${DIR}/check_fussion/output/$strain/$tissue/gene.gff 

                            cat ${DIR}/check_fussion/output/$strain/$tissue/gene.gff | gff2bed | cut -f1-4 > ${DIR}/check_fussion/output/$strain/$tissue/gene.bed

                            bedtools closest -a ${DIR}/check_fussion/output/$strain/$tissue/gene.bed -b ${DIR}/check_fussion/output/$strain/$tissue/gene.bed -d -N | cut -f 4,9 | sort -u > ${DIR}/check_fussion/tmp/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.$transcript.distance.tab

                            > "${DIR}/check_fussion/tmp/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.$transcript.genes.lst"

                            for l in `seq 1 $n`
                                   do
                                          sed "${l}q;d" ${DIR}/check_fussion/tmp/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.$transcript.merge.blastn | cut -f4 | tr ";" "\n" > ${DIR}/check_fussion/tmp/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.$transcript.gene$l.lst

                                          echo "${DIR}/check_fussion/tmp/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.$transcript.gene$l.lst" >> "${DIR}/check_fussion/tmp/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.$transcript.genes.lst"
                                   done

                            python3 ${DIR}/combinations.py ${DIR}/check_fussion/tmp/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.$transcript.genes.lst ${DIR}/check_fussion/tmp/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.$transcript.genesCombination.lst

                            nCombinations=$(cat ${DIR}/check_fussion/tmp/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.$transcript.genesCombination.lst |wc -l)

                            nGenesComb=$(cat ${DIR}/check_fussion/tmp/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.$transcript.genesCombination.lst | head -n1 | tr ' ' '\n' | wc -l)

                            > ${DIR}/check_fussion/tmp/$strain/$tissue/combination.lst

                            if [[ $nCombinations -gt 200 ]];then
                                   echo -e "$transcript\t$transcriptID\t$n\tproblem\tND\ttoo_many_combinations"
                            else
                                   if [[ $nGenesComb -eq 2  ]];then
                                          while read combination
                                                 do
                                                        echo "$combination" | tr ' ' '\n'  | sed 's/^/ID=/'  > ${DIR}/check_fussion/output/$strain/$tissue/gene.lst
                                                        grep -f ${DIR}/check_fussion/output/$strain/$tissue/gene.lst dmel-all-no-analysis-r6.31.gene.gff >  ${DIR}/check_fussion/output/$strain/$tissue/gene.gff
                                                        sed -i 's/;.*//' ${DIR}/check_fussion/output/$strain/$tissue/gene.gff 
                                                        cat ${DIR}/check_fussion/output/$strain/$tissue/gene.gff | gff2bed | cut -f1-4 > ${DIR}/check_fussion/output/$strain/$tissue/gene.bed
                                                        dist=$(bedtools closest -a ${DIR}/check_fussion/output/$strain/$tissue/gene.bed -b ${DIR}/check_fussion/output/$strain/$tissue/gene.bed -d -N | cut -f 4,9 | sort -u  | cut -f2 | sort -u)
                                                        chr=$(cut -f1 ${DIR}/check_fussion/output/$strain/$tissue/gene.bed | sort -u | wc -l)

                                                        if [[ $chr -eq 1 ]];then
                                                               #echo $transcript $dist
                                                               if [[ $dist -lt 2000 ]]; then
                                                               echo -e "$transcript\t$transcriptID\t$nMatches\t$combination\t$dist\tsame_chr_pervasive" >> ${DIR}/check_fussion/tmp/$strain/$tissue/combination.lst
                                                               else
                                                                      echo -e "$transcript\t$transcriptID\t$nMatches\t$combination\t$dist\tsame_chr_away" >> ${DIR}/check_fussion/tmp/$strain/$tissue/combination.lst            
                                                               fi
                                                        else
                                                               echo -e "$transcript\t$transcriptID\t$nMatches\t$combination\t$dist\tdiff_chr" >> ${DIR}/check_fussion/tmp/$strain/$tissue/combination.lst
                                                        fi
                                                 done < ${DIR}/check_fussion/tmp/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.$transcript.genesCombination.lst
                                                 result=$(cat ${DIR}/check_fussion/tmp/$strain/$tissue/combination.lst | grep "same_chr"   | tr ' ' '-' | sort -k4rh | head -n1 | wc -l)
                                                 if [[ $result -eq 1  ]];then
                                                        cat ${DIR}/check_fussion/tmp/$strain/$tissue/combination.lst | grep "same_chr" |  tr ' ' '-' | sort -k4rh | head -n1
                                                 else
                                                        echo -e "$transcript\t$transcriptID\t$n\tND\tND\tdiff_chr"
                                                 fi

                                   else
                                          > ${DIR}/check_fussion/tmp/$strain/$tissue/combination.lst
                                          while read combination
                                          do
                                                 echo "$combination" | tr ' ' '\n' | sed 's/^/ID=/'  > ${DIR}/check_fussion/output/$strain/$tissue/gene.lst
                                                 grep -f ${DIR}/check_fussion/output/$strain/$tissue/gene.lst dmel-all-no-analysis-r6.31.gene.gff >  ${DIR}/check_fussion/output/$strain/$tissue/gene.gff
                                                 sed -i 's/;.*//' ${DIR}/check_fussion/output/$strain/$tissue/gene.gff 
                                                 cat ${DIR}/check_fussion/output/$strain/$tissue/gene.gff | gff2bed | cut -f1-4 > ${DIR}/check_fussion/output/$strain/$tissue/gene.bed
                                                 sumDist=$(bedtools closest -a ${DIR}/check_fussion/output/$strain/$tissue/gene.bed -b ${DIR}/check_fussion/output/$strain/$tissue/gene.bed -d -N | cut -f 4,9 | sort -u  | cut -f2 | awk '{ sum += $1 } END { print sum }')
                                                 dist=$(bedtools closest -a ${DIR}/check_fussion/output/$strain/$tissue/gene.bed -b ${DIR}/check_fussion/output/$strain/$tissue/gene.bed -d -N | cut -f 4,9 | sort -u  | cut -f2 | tr '\n' ',')
                                                 distMax=$(bedtools closest -a ${DIR}/check_fussion/output/$strain/$tissue/gene.bed -b ${DIR}/check_fussion/output/$strain/$tissue/gene.bed -d -N | cut -f 4,9 | sort -u  | cut -f2 | sort -rh | head -n1 )
                                                 nChr=$(cut -f1 ${DIR}/check_fussion/output/$strain/$tissue/gene.bed | sort -u | wc -l)
                                                 echo -e "$transcript\t$transcriptID\t$nMatches\t$combination\t$sumDist\t$dist\t$distMax\t$nChr" >> ${DIR}/check_fussion/tmp/$strain/$tissue/combination.lst
                                          done < ${DIR}/check_fussion/tmp/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.$transcript.genesCombination.lst
                                          nResult=$(cat ${DIR}/check_fussion/tmp/$strain/$tissue/combination.lst | tr ' ' '-'  | awk ' $8 == 1 ' | sort -k4rh | head -n1 | wc -l)
                                          if [[ $nResult -eq 1  ]];then
                                                 distMax=$(cat ${DIR}/check_fussion/tmp/$strain/$tissue/combination.lst | tr ' ' '-'  | awk ' $8 == 1 ' | sort -k4rh | head -n1| cut -f 7)
                                                 if [[ $distMax -lt 2000 ]];then
                                                 a=$(cat ${DIR}/check_fussion/tmp/$strain/$tissue/combination.lst | tr ' ' '-'  | awk ' $8 == 1 ' | sort -k4rh | head -n1 | cut -f 1,2,3,4,7 )
                                                 echo -e "$a\tsame_chr_pervasive"
                                                 else
                                                 a=$(cat ${DIR}/check_fussion/tmp/$strain/$tissue/combination.lst | tr ' ' '-'  | awk ' $8 == 1 ' | sort -k4rh | head -n1 | cut -f 1,2,3,4,7 )
                                                 echo -e "$a\tsame_chr_away"             
                                                 fi
                                          else
                                                 echo -e "$transcript\t$transcriptID\t$n\tND\tND\tdiff_chr"
                                          fi

                                   
                                   fi
                            fi  


                     else 
                            awk -v transcript="$transcript" ' $1 == transcript '  ${DIR}/check_fussion/output/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.blastn | cut -f4 | sort -u | sed 's/^/ID=/' > ${DIR}/check_fussion/output/$strain/$tissue/gene.lst

                            grep -f ${DIR}/check_fussion/output/$strain/$tissue/gene.lst dmel-all-no-analysis-r6.31.gene.gff >  ${DIR}/check_fussion/output/$strain/$tissue/gene.gff

                            sed -i 's/;.*//' ${DIR}/check_fussion/output/$strain/$tissue/gene.gff 

                            cat ${DIR}/check_fussion/output/$strain/$tissue/gene.gff | gff2bed | cut -f1-4> ${DIR}/check_fussion/output/$strain/$tissue/gene.bed

                            distance=$(bedtools closest -a ${DIR}/check_fussion/output/$strain/$tissue/gene.bed -b ${DIR}/check_fussion/output/$strain/$tissue/gene.bed -d -N | cut -f 4,9 | sort -u | grep FBgn)

                            chr=$(cut -f1 ${DIR}/check_fussion/output/$strain/$tissue/gene.bed | sort -u | wc -l)
                            dist=0
                            if [[ $chr -eq 1 ]];then
                                   while read distGene
                                          do
                                                 distanceC=$(echo "$distGene" | cut -f2 )
                                                 if [[ $distanceC -gt 2000 ]];then
                                                        dist=1    
                                                 fi
                                   done <<< "$distance"

                                          if [[ $dist -eq 1 ]];then

                                                 while read distGene
                                                        do
                                                        gene=$(echo "$distGene" | cut -f1 )
                                                        dist=$(echo "$distGene" | cut -f2 )
                                                        echo -e "$transcript\t$transcriptID\t$nMatches\t$gene\t$dist\tsame_chr_away"
                                                 done <<< "$distance"
                                          else
                                                 while read distGene
                                                        do
                                                        gene=$(echo "$distGene" | cut -f1 )
                                                        dist=$(echo "$distGene" | cut -f2 )
                                                        echo -e "$transcript\t$transcriptID\t$nMatches\t$gene\t$dist\tsame_chr_pervasive"
                                                 done <<< "$distance"

                                          fi

                                         
                            else
                                   
                                   echo -e "$transcript\t$transcriptID\t$nMatches\tND\tND\tdiff_chr"
                            fi

                     fi

              elif [[ $nMatches -eq 1 ]]; then
                     gene=$(awk -v transcript="$transcript" ' $1 == transcript '  ${DIR}/check_fussion/output/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.blastn | cut -f4 | sort -u)
                     echo -e "$transcript\t$transcriptID\t$nMatches\t$gene\tND\tok"
              else
                     echo -e "$transcript\t$transcriptID\t$nMatches\t?\t?\t?"
              fi
       
       done < <(cut -f1 ${DIR}/check_fussion/output/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.blastn | sort -u ) > ${DIR}/check_fussion/output/$strain/$tissue/resultDistance.tab

echo -e "tissue\tstrain\tnTotal\tnGenes\tnGoodMatch\tnPervasive\tnPotentialChim\tnDiffChr\tnBigDistance\tcomb" > ${DIR}/check_fussion/result_identity.tab
nGoodMatch=$(awk ' $3 == 1 && $6 == "ok" || $3 == 1 && $6 == "duplicated_gene"  ' ${DIR}/check_fussion/output/$strain/$tissue/resultDistance.tab | cut -f1 | sort -u | wc -l )
nPotentialChim=$(awk ' $3 > 1  && $5 > 2000 && $6 != "too_many_combinations"  || $6 == "diff_chr" || $6 == "same_chr_away"  '  ${DIR}/check_fussion/output/$strain/$tissue/resultDistance.tab | cut -f1 | sort -u | wc -l )
nDiffChr=$(awk ' $3 > 1  && $6 == "diff_chr" '  ${DIR}/check_fussion/output/$strain/$tissue/resultDistance.tab | cut -f1 | sort -u | wc -l )
nBigDistance=$(awk ' $3 > 1  && $5 > 2000 && $6 == "same_chr_away" '  ${DIR}/check_fussion/output/$strain/$tissue/resultDistance.tab | cut -f1 | sort -u | wc -l )

nPervasive=$(awk ' $3 > 1  && $6  == "same_chr_pervasive" '  ${DIR}/check_fussion/output/$strain/$tissue/resultDistance.tab | cut -f1 | sort -u | wc -l )
nGenes=$(cut -f1 ${DIR}/check_fussion/output/$strain/$tissue/chimerics.filterLengthMin100.filterScovMin80.blastn | sort -u | wc -l )
comb=$(awk ' $3 > 1  && $6  == "too_many_combinations" '  ${DIR}/check_fussion/output/$strain/$tissue/resultDistance.tab | cut -f1 | sort -u | wc -l )
nTotal=$(cat ${DIR}/check_fussion/input/$strain/$tissue/chimerics.lst | sort -u | wc -l )
echo -e "$tissue\t$strain\t$nTotal\t$nGenes\t$nGoodMatch\t$nPervasive\t$nPotentialChim\t$nDiffChr\t$nBigDistance\t$comb" >> ${DIR}/check_fussion/result_identity.tab


# blastn -query ${DIR}/check_fussion/input/$strain/$tissue/chimerics.fasta \
#        -db dmel-all-gene-r6.31.fasta \
#        -out ${DIR}/check_fussion/output/$strain/$tissue/chimerics.blastn \
#        -num_threads 4 \
#        -outfmt "6 qseqid qlen sseqid slen qcovs pident length mismatch gapopen qstart qend sstart send evalue bitscore"

# while read line; do
# if echo "$line" | grep -q "^>"; then
# id=$(echo $line | awk -F'ID=' '{print $2}' | awk -F';' '{print $1}')
# parent=$(echo $line | awk -F'parent=' '{print $2}' | awk -F';' '{print $1}')
# new_id="$id":"$parent"
# echo ">$new_id"
# else
# echo "$line"
# fi
# done < <(cat dmel-all-transcript-r6.31.fasta) > dmel-all-transcript-r6.31.transf.fasta