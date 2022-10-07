#!/bin/bash

# set the partition where the job will run

#SBATCH --partition=normal

# set the number of nodes

# mail alert at start, end and abortion of execution

#SBATCH --mail-type=ALL

#SBATCH --cpus-per-task=1

#SBATCH --mem=4G

#SBATCH --job-name=parseGO

#SBATCH --output=parseGO.output

#SBATCH --error=parseGO.error

# send mail to this address

#SBATCH --mail-user=marta.coronado@ibe.upf-csic.es

module load Miniconda3/4.7.10

source activate /homes/users/mcoronado/.conda/envs/trinity

# Run

DIR=/homes/users/mcoronado/scratch/5GenomesProject

tissues="head gut ovary"
#tissue="gut"

#nClusters=$(grep "Annotation Cluster" ${tissue}ResultGO.tsv | wc -l )
nClusters=5
echo -e "cluster\ttissue\tnGenes\tenrichmentScore\tshortName\tGOnames" > resultParse.tsv
for tissue in $tissues
	do
		for n in $(seq 1 $nClusters); do
			nextCluster=$((n + 1))
			#awk "/Annotation Cluster $n/{flag=3;next}/Annotation Cluster $nCluster/{flag=0}flag" ${tissue}ResultGO.tsv
			nGenes=$(grep -ozP "(?s)Annotation Cluster ${n}\t.*(Annotation Cluster ${nextCluster}\t)" ${tissue}ResultGO.tsv | grep "FBGN" | cut -f 6 | sed "s/, /\n/g" | sort -u | wc -l)
			nameGO=$(grep -ozP "(?s)Annotation Cluster ${n}\t.*(Annotation Cluster ${nextCluster}\t)" ${tissue}ResultGO.tsv | grep GOTERM_BP_ALL | cut -f2 | cut -f2 -d'~' | tr "\n" ";")
			enrichment=$(grep -P "Annotation Cluster ${n}\t" ${tissue}ResultGO.tsv | cut -f2 | cut -f2 -d":" | cut -c2-)
			echo -e "$n\t$tissue\t$nGenes\t$enrichment\t\t$nameGO" >> resultParse.tsv
			echo "#######"
			echo cluster $n $tissue
			#grep -ozP "(?s)Annotation Cluster ${n}\t.*(Annotation Cluster ${nextCluster}\t)" ${tissue}ResultGO.tsv | grep GOTERM_BP_ALL | cut -f2 | cut -f2 -d'~' 
			words=$(grep -ozP "(?s)Annotation Cluster ${n}\t.*(Annotation Cluster ${nextCluster}\t)" ${tissue}ResultGO.tsv | grep GOTERM_BP_ALL | cut -f2 | cut -f2 -d'~' | tr " " "\n" | awk ' { tot[$0]++ } END { for (i in tot) print tot[i],i } '  | sort)
			#echo $words
			echo ""
		done
	done 

n=1
for tissue in $tissues
	do
		echo $tissue
		> ${tissue}_revigo.input
	    #nClusters=$(grep "Enrichment Score" ${tissue}ResultGO.tsv  | awk '$6>1.3' | cut -f 1 | cut -f3 -d' ' | tail -n1)
	    #nClusters=$(($nClusters + 1))
		nClusters=6
		while read line
			do
				GOID=$(echo "$line" | cut -f2  | cut -f1 -d'~')
				pvalue=$(echo "$line" | cut -f5)
				echo $GOID $pvalue >> ${tissue}_revigo.input
			done < <(grep -ozP "(?s)Annotation Cluster ${n}\t.*(Annotation Cluster ${nClusters}\t)" ${tissue}ResultGO.tsv | grep GOTERM_BP_ALL)
	done