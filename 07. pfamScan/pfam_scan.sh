#!/bin/bash

# set the partition where the job will run

#SBATCH --partition=express

# set the number of nodes

# mail alert at start, end and abortion of execution

#SBATCH --mail-type=ALL

#SBATCH --cpus-per-task=2

#SBATCH --mem=12G

#SBATCH --job-name=PFAMScan

#SBATCH --output=logs/PFAMScan_%a.output

#SBATCH --error=logs/PFAMScan_%a.error

source ~/miniconda3/etc/profile.d/conda.sh
conda activate pfamDomains

# Run

DIR="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/chimerics"
DIRDATA="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/DATA"

strain=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/data.txt) | cut -f 1)
tissue=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/data.txt) | cut -f 2)


#assemblies="AKA-017"
#tissues="head"
cd $DIR/pfamDomains

#rm -r $strain/$tissue/tmp/gff
#rm -r $strain/$tissue/tmp/fasta
#rm $strain/$tissue/*output
rm -rf $strain/$tissue/tmp/gff
rm -rf $strain/$tissue/tmp/fasta
rm -rf $strain/$tissue/pfamOutput

mkdir -p $strain/$tissue/tmp/gff
mkdir -p $strain/$tissue/tmp/fasta
mkdir -p $strain/$tissue/pfamOutput

insertions=$(grep $tissue ${DIR}/clean_set_minimap2/chimerics_v2.tab | grep $strain | awk  -F $'\t'  ' $22 >= 0.39 && $23 == "mRNA" '  |  cut -f  1,2,4,5,7,9,16,22,24,23,28 | sort -u )

# insertion=$(grep $tissue ${DIR}/clean_set_minimap2/chimerics_v2.tab | grep $strain | awk ' $22 >= 0.39 && $23 == "mRNA" ' |  cut -f  1,2,4,5,6,7,9,16,22,24,23,27,13,26  | head -n1)

i=1

while IFS= read -r insertion
do
trinity=$(echo "$insertion" | cut -f3)
stringtieID=$(echo "$insertion" | cut -f5)
geneID=$(echo "$insertion" | cut -f6)
TEconsensusID=$(echo "$insertion" | cut -f7)
TEfam=$(echo "$insertion" | cut -f4)
#pos=$(echo "$insertion" | cut -f7 | tr -d ' ')
CP=$(echo "$insertion" | cut -f8)
expr=$(echo "$insertion" | cut -f10)
#SS=$(echo "$insertion" | cut -f11)
#group=$(echo "$insertion" | cut -f12)
id=$(echo "$insertion" | cut -f11)

grep -w "$trinity" $DIR/RepeatMasker/${strain}_${tissue}/Trinity.hits_bin90.fasta.out.gff | grep -w "$TEconsensusID" | mergeBed -d 25 > $strain/$tissue/tmp/gff/${trinity}_${TEfam}.${i}.gff

n=$(wc -l $strain/$tissue/tmp/gff/${trinity}_${TEfam}.${i}.gff | cut -f1 -d' ')

if [[ $n -eq 0 ]]; then

consensus=$(grep -w "$trinity" $DIR/RepeatMasker/${strain}_${tissue}/Trinity.hits_bin90.fasta.out.gff | cut -f 9 | cut -f2 -d' ' | tr -d '"' | cut -f2 -d':' | sort -u)
fam=$(grep "$consensus" ${DIRDATA}/TEs/TE_library_family.csv | cut -f 8)

if [[ $TEfam == $fam ]]; then
grep -w "$trinity" $DIR/RepeatMasker/${strain}_${tissue}/Trinity.hits_bin90.fasta.out.gff | grep -w "$consensus" | mergeBed -d 25 > $strain/$tissue/tmp/gff/${trinity}_${TEfam}.${i}.gff
n=$(wc -l $strain/$tissue/tmp/gff/${trinity}_${TEfam}.${i}.gff | cut -f1 -d' ')
echo $trinity $strain $tissue ${n}R $lengthTE >> n
lengthTE=$(awk '{sum += $3 - $2 + 1} END {print sum}' $strain/$tissue/tmp/gff/${trinity}_${TEfam}.${i}.gff )
fi

else
lengthTE=$(awk '{sum += $3 - $2 + 1} END {print sum}' $strain/$tissue/tmp/gff/${trinity}_${TEfam}.${i}.gff )

echo $trinity $strain $tissue $n $lengthTE >> n


fi


for seq in `seq 1 $n`
do
sed "${seq}q;d"  $strain/$tissue/tmp/gff/${trinity}_${TEfam}.${i}.gff  > $strain/$tissue/tmp/gff/${trinity}_${TEfam}.${seq}.${i}.gff 

bedtools getfasta -fi $DIR/RepeatMasker/${strain}_${tissue}/Trinity.hits_bin90.fasta.masked -bed $strain/$tissue/tmp/gff/${trinity}_${TEfam}.${seq}.${i}.gff | sed "s/$trinity:.*/${id}_${strain}_${tissue}_${TEfam}_${lengthTE}_${trinity}/g" > $strain/$tissue/tmp/fasta/${trinity}_${TEfam}.${seq}.${i}.fasta 
getorf -sequence $strain/$tissue/tmp/fasta/${trinity}_${TEfam}.${seq}.${i}.fasta  -outseq $strain/$tissue/tmp/fasta/${trinity}_${TEfam}.ORF.${seq}.${i}.fasta

#faidx  --transform chromsizes $strain/$tissue/tmp/fasta/${trinity}_${TEfam}.ORF.${seq}.fasta | sort -k2,2nr | head -n1 |cut -f1 > $strain/$tissue/tmp/fasta/longestIsoform_${trinity}_${TEfam}.${seq}.lst
#seqtk subseq $strain/$tissue/tmp/fasta/${trinity}_${TEfam}.ORF.${seq}.fasta $strain/$tissue/tmp/fasta/longestIsoform_${trinity}_${TEfam}.${seq}.lst > $strain/$tissue/tmp/fasta/${trinity}_${TEfam}.longORF.${seq}.fasta
#sed -i "1 s/.*/>${trinity}:${tissue}:${strain}:${stringtieID}:${geneID}:${TEfam}:g${group}/" $strain/$tissue/tmp/fasta/${trinity}_${TEfam}.longORF.${seq}.fasta
pfam_scan.pl -fasta $strain/$tissue/tmp/fasta/${trinity}_${TEfam}.ORF.${seq}.${i}.fasta  -dir pfamFiles -cpu 2 -outfile $strain/$tissue/pfamOutput/${trinity}_${TEconsensusID}_g${group}.${seq}.${i}.pfam.output
done
i=$(( i + 1 ))

done <<< "$insertions"

cat $strain/$tissue/pfamOutput/*.pfam.output |grep -v "#"  | grep -v -e '^$'  | sed -e 's/\s\+/ /g' | tr " " "\t" > globalResult_${strain}_${tissue}.tsv

# cat globalResult_*.tsv | sed 's/_/\t/'  | sed 's/_/\t/'  | sed 's/_/\t/'  | sed 's/_/\t/'  | sed 's/_TR/\tTR/' | sed 's/_/\t/5'  > globalResult.tsv


# <seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan>