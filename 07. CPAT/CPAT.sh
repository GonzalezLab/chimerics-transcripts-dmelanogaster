#!/bin/bash

# set the partition where the job will run

#SBATCH --partition=normal

# set the number of nodes

# mail alert at start, end and abortion of execution

#SBATCH --mail-type=ALL

#SBATCH --cpus-per-task=2

#SBATCH --mem=8G

#SBATCH --job-name=PFAMScan

#SBATCH --output=PFAMScan.output

#SBATCH --error=PFAMScan.error

# send mail to this address

#SBATCH --mail-user=marta.coronado@ibe.upf-csic.es

module load Miniconda3/4.7.10

source activate /homes/users/mcoronado/.conda/envs/trinity

# Run
# To install CPAT first (version 3.0.3)
# conda install cpat
# pip3 install CPAT --upgrade

# PREBUILT MODELS FOR DROSOPHILA
# https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/fly_Hexamer.tsv/download
# https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/Fly_logitModel.RData/download
DATA=/homes/users/mcoronado/scratch/5GenomesProject/DATA/FlyBase_r6.31
DIR=/homes/users/mcoronado/scratch/5GenomesProject/ANALYSES/trinity-repeatMasker
assemblies="AKA-017 JUT-011 MUN-016 SLA-001 TOM-007"
tissues="head gut ovary"
#assemblies="AKA-017"
#tissues="head"

# min ORF default=75

# for assembly in $assemblies
# 	do
# 		for tissue in $tissues
# 			do
# 				mkdir -p CPAT/${assembly}/${tissue}
# 				cpat.py -x flyPrebuild/fly_Hexamer.tsv -d flyPrebuild/Fly_logitModel.RData \
# 				--top-orf=100 \
# 				--antisense \
# 				-g ${DIR}/Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/Trinity.hits_bin90.fasta \
# 				-o CPAT_${assembly}_${tissue}_output

# 				mv CPAT_run_info.log 
# 				mv CPAT_${assembly}_${tissue}_output.no_ORF.txt CPAT/${assembly}/${tissue}
# 				mv CPAT_${assembly}_${tissue}_output.ORF_prob.best.tsv CPAT/${assembly}/${tissue}
# 				mv CPAT_${assembly}_${tissue}_output.ORF_prob.tsv CPAT/${assembly}/${tissue}
# 				mv CPAT_${assembly}_${tissue}_output.ORF_seqs.fa CPAT/${assembly}/${tissue} 
# 				mv CPAT_${assembly}_${tissue}_output.r CPAT/${assembly}/${tissue}
# 				mv CPAT_run_info.log CPAT/${assembly}/${tissue}
# 			done
# 	done

# transcripts with middle exon
cat ../../TE_chimeric_global_REVISED_v2.tab |  grep  "Middle exon" |  cut -f4 |sort -u

# only chimerics
for assembly in $assemblies
do
	for tissue in $tissues
	do
		mkdir -p CPAT_TE/${assembly}/${tissue}
		grep ${assembly} ../../TE_chimeric_global_REVISED_v2.tab | grep ${tissue} | grep  "Middle exon" |  cut -f4 |sort -u > CPAT_TE/${assembly}/${tissue}/genesTE_${assembly}_${tissue}.lst
		seqtk subseq ${DIR}/Trinity_${assembly}_${tissue}/RepeatMasker_Blast_Merged_Reference/Trinity.hits_bin90.fasta CPAT_TE/${assembly}/${tissue}/genesTE_${assembly}_${tissue}.lst > CPAT_TE/${assembly}/${tissue}/Trinity.hits_bin90_TE_${assembly}_${tissue}.fasta
		cpat.py -x flyPrebuild/fly_Hexamer.tsv -d flyPrebuild/Fly_logitModel.RData \
		--top-orf=100 \
		--antisense \
		-g CPAT_TE/${assembly}/${tissue}/Trinity.hits_bin90_TE_${assembly}_${tissue}.fasta \
		-o CPAT_${assembly}_${tissue}_output

		mv CPAT_${assembly}_${tissue}_output.no_ORF.txt CPAT_TE/${assembly}/${tissue}
		mv CPAT_${assembly}_${tissue}_output.ORF_prob.best.tsv CPAT_TE/${assembly}/${tissue}
		mv CPAT_${assembly}_${tissue}_output.ORF_prob.tsv CPAT_TE/${assembly}/${tissue}
		mv CPAT_${assembly}_${tissue}_output.ORF_seqs.fa CPAT_TE/${assembly}/${tissue} 
		mv CPAT_${assembly}_${tissue}_output.r CPAT_TE/${assembly}/${tissue}
		mv CPAT_run_info.log CPAT_TE/${assembly}/${tissue}
	done
done

for assembly in $assemblies
do
	for tissue in $tissues
	do
		cat CPAT_TE/${assembly}/${tissue}/CPAT_${assembly}_${tissue}_output.no_ORF.txt
	done
done | wc -l
#We have info for all transcripts analyzed (no WARNING, No ORFs found for XXX)
> CPAT_TE.ORF_prob.best.INFO.tsv
for assembly in $assemblies
do
	for tissue in $tissues
	do
		header=$(head -n1 CPAT_TE/${assembly}/${tissue}/CPAT_${assembly}_${tissue}_output.ORF_prob.best.tsv)
		echo -e "$header\tstatus\tstringtieID\tFlyBaseGeneRef\ttissue\tassembly" > CPAT_TE/${assembly}/${tissue}/CPAT_${assembly}_${tissue}_output.ORF_prob.best.INFO.tsv

		echo $tissue $assembly

		while IFS= read -r ORF
		do
			transcript=$(echo "$ORF" | cut -f1)
			FlyBaseID=$(grep ${assembly}  ../../TE_chimeric_global_REVISED_v2.tab | grep ${tissue} | grep -wi ${transcript} | cut -f 7 | sort -u)
			FlyBaseIDgene=$(grep ${assembly}  ../../TE_chimeric_global_REVISED_v2.tab | grep ${tissue} |  grep -wi ${transcript} | cut -f 8 | sort -u)
			stringtieID=$(grep ${assembly}  ../../TE_chimeric_global_REVISED_v2.tab | grep ${tissue} | grep -wi ${transcript} | cut -f 6 | sort -u)
			status=$(grep "$FlyBaseID" ${DATA}/dmel-chr-r6.31.gtf | head -n1 | cut -f3)
			echo -e "$ORF\t$status\t$stringtieID\t$FlyBaseIDgene\t$tissue\t$assembly"
		done < <(tail -n+2 CPAT_TE/${assembly}/${tissue}/CPAT_${assembly}_${tissue}_output.ORF_prob.best.tsv) >> CPAT_TE/${assembly}/${tissue}/CPAT_${assembly}_${tissue}_output.ORF_prob.best.INFO.tsv

		tail -n +2 CPAT_TE/${assembly}/${tissue}/CPAT_${assembly}_${tissue}_output.ORF_prob.best.INFO.tsv >> CPAT_TE.ORF_prob.best.INFO.tsv
	done
done 

header=$(head -n1 CPAT_TE/AKA-017/head/CPAT_AKA-017_head_output.ORF_prob.best.INFO.tsv)	
sed -i "1i $header\tCP" CPAT_TE.ORF_prob.best.INFO.tsv

mv CPAT_TE.ORF_prob.best.INFO.tsv CPAT_TE.ORF_prob.best.INFO.tmp.tsv
awk -F"\t" 'NR>1 {$0 = $0 FS (($11 >= 0.39) ? "codingPotential" : "nonCodingPotential")} 1' CPAT_TE.ORF_prob.best.INFO.tmp.tsv > CPAT_TE.ORF_prob.best.INFO.tsv 

sed -i "1 s/.*/$header\tCP/" CPAT_TE.ORF_prob.best.INFO.tsv 

for assembly in $assemblies
do
	for tissue in $tissues
	do
		cut -f12 CPAT_TE/${assembly}/${tissue}/CPAT_${assembly}_${tissue}_output.ORF_prob.best.INFO.tsv | sort -u
	done
done | sort -u

# mrna
# pseudogene, ncRNA

echo -e "assembly\ttissue\t# transcripts TE\ttranscripts coding potential (>= 0.39)\ttranscripts non-coding potential (<0.39)\ttranscripts coding potential (mRNA)\ttranscripts coding potential (ncRNA, pre_miRNA, snoRNA, tRNA, pseudogene)\ttranscripts non-coding potential (mRNA)\tntranscripts non-coding potential  (ncRNA, pre_miRNA, snoRNA, tRNA, pseudogene)\tpercentage coding\tp-value (enrichment of non-coding in <0.39)" > CPAT_results.tsv

for assembly in $assemblies
do
	for tissue in $tissues
	do
	#echo $tissue $assembly
	nTranscripts=$(tail -n +2 CPAT_TE/${assembly}/${tissue}/CPAT_${assembly}_${tissue}_output.ORF_prob.best.INFO.tsv  | wc -l )
	nTranscriptsCodingPot=$(tail -n +2 CPAT_TE/${assembly}/${tissue}/CPAT_${assembly}_${tissue}_output.ORF_prob.best.INFO.tsv | awk '$11>=0.39' | wc -l)
	nTranscriptsCodingPotProtein=$(tail -n +2 CPAT_TE/${assembly}/${tissue}/CPAT_${assembly}_${tissue}_output.ORF_prob.best.INFO.tsv | awk '$11>=0.39' | grep "mRNA" | wc -l )
	nTranscriptsCodingPotNotProtein=$(tail -n +2 CPAT_TE/${assembly}/${tissue}/CPAT_${assembly}_${tissue}_output.ORF_prob.best.INFO.tsv | awk '$11>=0.39' | grep -v "mRNA" | wc -l)
	nTranscriptsNonCodingPot=$(tail -n +2 CPAT_TE/${assembly}/${tissue}/CPAT_${assembly}_${tissue}_output.ORF_prob.best.INFO.tsv  | awk '$11<0.39' | wc -l)
	nTranscriptsNonCodingPotProtein=$(tail -n +2 CPAT_TE/${assembly}/${tissue}/CPAT_${assembly}_${tissue}_output.ORF_prob.best.INFO.tsv  | awk '$11<0.39' | grep "mRNA" | wc -l)
	nTranscriptsNonCodingPotNotProtein=$(tail -n +2 CPAT_TE/${assembly}/${tissue}/CPAT_${assembly}_${tissue}_output.ORF_prob.best.INFO.tsv  | awk '$11<0.39' | grep -v "mRNA" | wc -l)
	RES=$(Rscript proportions.R $nTranscripts $nTranscriptsCodingPot $nTranscriptsCodingPotProtein $nTranscriptsCodingPotNotProtein $nTranscriptsNonCodingPot $nTranscriptsNonCodingPotProtein $nTranscriptsNonCodingPotNotProtein)
	per=$(echo "$RES" | cut -f2 -d' ' | tr -d '"')
	pval=$(echo "$RES" | cut -f3 -d' ' | tr -d '"')

	echo -e "$assembly\t$tissue\t$nTranscripts\t$nTranscriptsCodingPot\t$nTranscriptsNonCodingPot\t$nTranscriptsCodingPotProtein\t$nTranscriptsCodingPotNotProtein\t$nTranscriptsNonCodingPotProtein\t$nTranscriptsNonCodingPotNotProtein\t$per\t$pval"
done
done >> CPAT_results.tsv




