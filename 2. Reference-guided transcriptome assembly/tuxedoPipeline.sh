#!/bin/bash

# set the partition where the job will run

#SBATCH --partition=normal

# set the number of nodes

#SBATCH --cpus-per-task=12

#SBATCH --mem=75G

# mail alert at start, end and abortion of execution

#SBATCH --mail-type=ALL

#SBATCH --job-name=TuxedoPipeline

#SBATCH --output=TuxedoPipeline.output

#SBATCH --error=TuxedoPipeline.error

# send mail to this address

#SBATCH --mail-user=marta.coronado@ibe.upf-csic.es

module load Miniconda3/4.7.10

source activate /homes/users/mcoronado/.conda/envs/5GP

# Run

DIR=/homes/users/mcoronado/scratch/5GenomesProject

# ## PART 1 - ALIGNMENT
mkdir ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/annotations/


extract_splice_sites.py ${DIR}/DATA/FlyBase_r6.31/dmel-chr-r6.31.gtf > ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/annotations/reference.ss
extract_exons.py ${DIR}/DATA/FlyBase_r6.31/dmel-chr-r6.31.gtf > ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/annotations/reference.exon


mkdir ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/index/

hisat2-build -p 12 --ss ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/annotations/reference.ss --exon ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/annotations/reference.exon ${DIR}/DATA/FlyBase_r6.31/dmel-chr.fasta ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/index/reference


mkdir ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map

# Head

for assembly in `less ${DIR}/DATA/genomeAssembly/genomeAssembly.txt`
      do
          mkdir ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}
          mkdir ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}/Head
          for dataRna in `less ${DIR}/DATA/RNAseqData/files/${assembly}_RNAseq_head.txt`
              do
              	hisat2 -p 12 --dta -x ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/index/reference -1 ${DIR}/DATA/RNAseqData/Heads/${dataRna}_read1.fastq.gz -2 ${DIR}/DATA/RNAseqData/Heads/${dataRna}_read2.fastq.gz --summary-file ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}/Head/${dataRna}.log -S ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}/Head/${dataRna}.sam 
                samtools sort -@ 12 -o ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}/Head/${dataRna}.bam ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}/Head/${dataRna}.sam
                 # rm ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}/Head/${dataRna}.sam 
              done
      done

# Gut
for assembly in `less ${DIR}/DATA/genomeAssembly/genomeAssembly.txt`
    do
        mkdir ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}/Gut
        for dataRna in `less ${DIR}/DATA/RNAseqData/files/${assembly}_RNAseq_gut.txt`
            do
                hisat2 -p 12 --dta -x ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/index/reference -1 ${DIR}/DATA/RNAseqData/Gut/${dataRna}_read1.fastq.gz -2 ${DIR}/DATA/RNAseqData/Gut/${dataRna}_read2.fastq.gz --summary-file ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}/Gut/${dataRna}.log -S ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}/Gut/${dataRna}.sam
                samtools sort -@ 12 -o ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}/Gut/${dataRna}.bam ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}/Gut/${dataRna}.sam 
                # rm ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}/Gut/${dataRna}.sam 
            done
    done

# Ovary
for assembly in `less ${DIR}/DATA/genomeAssembly/genomeAssembly.txt`
    do
        mkdir ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}/Ovary
        for dataRna in `less ${DIR}/DATA/RNAseqData/files/${assembly}_RNAseq_ovary.txt`
            do
                hisat2 -p 12 --dta -x ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/index/reference -1 ${DIR}/DATA/RNAseqData/Ovary/${dataRna}_read1.fastq.gz -2 ${DIR}/DATA/RNAseqData/Ovary/${dataRna}_read2.fastq.gz --summary-file ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}/Ovary/${dataRna}.log -S ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}/Ovary/${dataRna}.sam 
                samtools sort -@ 12 -o ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}/Ovary/${dataRna}.bam ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}/Ovary/${dataRna}.sam 
                # ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}/Ovary/${dataRna}.sam 
            done
    done

# IMPROVED ASSEMBLY
mkdir ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/assemblyImproved

# Head
for assembly in `less ${DIR}/DATA/genomeAssembly/genomeAssembly.txt`
    do
        mkdir ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/assemblyImproved/${assembly}
        mkdir ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/assemblyImproved/${assembly}/Head
        for dataRna in `less ${DIR}/DATA/RNAseqData/files/${assembly}_RNAseq_head.txt`
            do
                stringtie -c 1.5 -g 51 -f 0.016 -j 2 -a 15 -M 0.95 ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}/Head/${dataRna}.bam -l ${dataRna} -p 12 -G ${DIR}/DATA/FlyBase_r6.31/dmel-chr-r6.31.gtf -o ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/assemblyImproved/${assembly}/Head/${dataRna}.gtf
            done
    done

# Gut
for assembly in `less ${DIR}/DATA/genomeAssembly/genomeAssembly.txt`
    do
        mkdir ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/assemblyImproved/${assembly}/Gut
        for dataRna in `less ${DIR}/DATA/RNAseqData/files/${assembly}_RNAseq_gut.txt`
            do
                stringtie -c 1.5 -g 51 -f 0.016 -j 2 -a 15 -M 0.95 ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}/Gut/${dataRna}.bam -l ${dataRna} -p 12 -G ${DIR}/DATA/FlyBase_r6.31/dmel-chr-r6.31.gtf -o ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/assemblyImproved/${assembly}/Gut/${dataRna}.gtf
            done
    done

# Ovary
for assembly in `less ${DIR}/DATA/genomeAssembly/genomeAssembly.txt`
    do
        mkdir ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/assemblyImproved/${assembly}/Ovary
        for dataRna in `less ${DIR}/DATA/RNAseqData/files/${assembly}_RNAseq_ovary.txt`
            do
                stringtie -c 1.5 -g 51 -f 0.016 -j 2 -a 15 -M 0.95 ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/map/${assembly}/Ovary/${dataRna}.bam -l ${dataRna} -p 12 -G ${DIR}/DATA/FlyBase_r6.31/dmel-chr-r6.31.gtf -o ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/assemblyImproved/${assembly}/Ovary/${dataRna}.gtf
            done
    done




stringtie --merge -F 0 -T 10 -g 0 -c 1.5 -f 0.016 -p 12 -G ${DIR}/DATA/FlyBase_r6.31/dmel-chr-r6.31.gtf -o ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/assemblyImproved/stringtie_merged.gtf ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/mergelistImproved.txt
mkdir ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/assemblyImproved/mergedSamples
cd ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/assemblyImproved/mergedSamples
gffcompare -r ${DIR}/DATA/FlyBase_r6.31/dmel-chr-r6.31.gtf -G -o mergedSamples ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/assemblyImproved/stringtie_merged.gtf

gffread -w ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/assemblyImproved/stringtie_merged.fasta -g ${DIR}/DATA/FlyBase_r6.31/dmel-chr.fasta ${DIR}/ANALYSES/deNovoTranscriptAnnotationReference/assemblyImproved/stringtie_merged.gtf -F
