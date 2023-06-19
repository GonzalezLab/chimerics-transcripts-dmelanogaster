#!/bin/bash

#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --job-name=trinity
#SBATCH --output=logs/trinity_rm_%a.out
#SBATCH --error=logs/trinity_rm_%a.err


source ~/miniconda3/etc/profile.d/conda.sh
conda activate chimerics

DIR="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/chimerics"
DIRDATA="/lustre/home/ibe/mcoronado/scratch/5GenomesProject/DATA"

cd ${DIR}

strain=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/data.txt) | cut -f 1)
tissue=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat $DIR/data.txt) | cut -f 2)

mkdir ${DIR}/input_trinity
mkdir -p trinity/trinity_${strain}_${tissue}

while read repData
do
rep=$(echo "$repData" | cut -f1)
echo -e "${tissue}_${strain}\t$rep\tDATA/fastp/${rep}_trim_1.fq.gz\tDATA/fastp/${rep}_trim_2.fq.gz"
done < <(grep "$strain" ${DIR}/samples.txt  | grep "$tissue") > ${DIR}/input_trinity/samples.${strain}.${tissue}.txt

cd /lustre/home/ibe/mcoronado/scratch/5GenomesProject/

# singularity exec --bind "/lustre/scratch-global/ibe/mcoronado/5GenomesProject/" -e trinityrnaseq.v2.15.1.simg Trinity \
# --seqType fq \
# --max_memory 78G --CPU 12 \
# --samples_file ANALYSIS/chimerics/input_trinity/samples.${strain}.${tissue}.txt \
# --jaccard_clip \
# --output ANALYSIS/chimerics/trinity/trinity_${strain}_${tissue}

# singularity exec --bind "/lustre/scratch-global/ibe/mcoronado/5GenomesProject/" -e trinityrnaseq.v2.15.1.simg /usr/local/bin/util/TrinityStats.pl ANALYSIS/chimerics/trinity/trinity_${strain}_${tissue}/Trinity.fasta > ANALYSIS/chimerics/trinity/trinity_${strain}_${tissue}/Trinity.stats

# blastn -query ${DIR}/trinity/trinity_${strain}_${tissue}/Trinity.fasta \
#  	-db ${DIR}/DataAssemblyMergedReference/stringtie_merged.fasta \
#  	-out ${DIR}/trinity/trinity_${strain}_${tissue}/Trinity_refTrans.blastn \
#     -evalue 1e-20 -dust no -task megablast -num_threads 12 \
#     -max_target_seqs 1 -outfmt 6

# singularity exec --bind "/lustre/scratch-global/ibe/mcoronado/5GenomesProject/" -e trinityrnaseq.v2.15.1.simg /usr/local/bin/util/analyze_blastPlus_topHit_coverage.pl \
#     ANALYSIS/chimerics/trinity/trinity_${strain}_${tissue}/Trinity_refTrans.blastn \
#     ANALYSIS/chimerics/trinity/trinity_${strain}_${tissue}/Trinity.fasta \
#     ANALYSIS/chimerics/DataAssemblyMergedReference/stringtie_merged.fasta

# grep -E "Bin_90|Bin_100" ${DIR}/trinity/trinity_${strain}_${tissue}/Trinity_refTrans.blastn.hist.list | cut -f2,3  > ${DIR}/trinity/trinity_${strain}_${tissue}/Trinity_refTrans.blastn.hist_bin90.list.names

# > ${DIR}/trinity/trinity_${strain}_${tissue}/Trinity_refTrans.blastn.hist_bin90.list.names.tab

# while read transcript
# do
# trinityID=$(echo "$transcript" | cut -f1)
# stringtieID=$(echo "$transcript" | cut -f2)
# transcriptID=$(awk -v stringtieID="$stringtieID" ' $5 == stringtieID ' ${DIR}/DataAssemblyMergedReference/mergedSamples.stringtie_merged.gtf.tmap | cut -f2 | sort -u )
# geneID=$(awk -v stringtieID="$stringtieID" ' $5 == stringtieID ' ${DIR}/DataAssemblyMergedReference/mergedSamples.stringtie_merged.gtf.tmap | cut -f1 | sort -u )
# type=$(awk -v stringtieID="$stringtieID" ' $5 == stringtieID ' ${DIR}/DataAssemblyMergedReference/mergedSamples.stringtie_merged.gtf.tmap | cut -f3 | sort -u )
# echo -e "$trinityID\t$stringtieID\t$type\t$geneID\t$transcriptID" >> ${DIR}/trinity/trinity_${strain}_${tissue}/Trinity_refTrans.blastn.hist_bin90.list.names.tab
# done < <(cat ${DIR}/trinity/trinity_${strain}_${tissue}/Trinity_refTrans.blastn.hist_bin90.list.names)

# keep = j k m n o x y

awk ' $3 == "=" ||  $3 == "j" ||  $3 == "k" ||  $3 == "m" ||  $3 == "n" ||  $3 == "o" ||  $3 == "x" ||  $3 == "y" ' ${DIR}/trinity/trinity_${strain}_${tissue}/Trinity_refTrans.blastn.hist_bin90.list.names.tab >  ${DIR}/trinity/trinity_${strain}_${tissue}/Trinity_refTrans.blastn.hist_bin90.list.names.keep.tab 

cut -f1 ${DIR}/trinity/trinity_${strain}_${tissue}/Trinity_refTrans.blastn.hist_bin90.list.names.keep.tab | sort -u > ${DIR}/trinity/trinity_${strain}_${tissue}/Trinity_refTrans.blastn.hist_bin90.list.names.keep 

seqtk subseq ${DIR}/trinity/trinity_${strain}_${tissue}/Trinity.fasta ${DIR}/trinity/trinity_${strain}_${tissue}/Trinity_refTrans.blastn.hist_bin90.list.names.keep > ${DIR}/trinity/trinity_${strain}_${tissue}/Trinity.hits_bin90.fasta 

mkdir -p ${DIR}/RepeatMasker/${strain}_${tissue}/

cd ${DIR}/RepeatMasker/${strain}_${tissue}/

RepeatMasker -pa 4 -e rmblast -lib ${DIRDATA}/TEs/consensuses_curated_v4.fasta -dir ${DIR}/RepeatMasker/${strain}_${tissue}/ -norna -nolow -s -cutoff 250 -xsmall -no_is -gff ${DIR}/trinity/trinity_${strain}_${tissue}/Trinity.hits_bin90.fasta 

# grep Unspecified RepeatMaskerImproved90/Trinity.hits_bin90.fasta.out | sed 's/  */ /g'  | sed "s/^ //g" | tr ' ' '\t' | sort -nr -k1  > RepeatMaskerImproved90/Trinity.hits_bin90.fasta.out.TE
