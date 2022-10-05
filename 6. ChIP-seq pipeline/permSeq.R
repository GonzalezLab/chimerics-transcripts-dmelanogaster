#To execute this script the following modules should be loaded:

# module load R/3.5.1-foss-2018b
# module load csem/2.3-Perl-5.26.1-foss-2018a
# module load SRA-Toolkit/2.8.2-1-centos_linux64
# module load FastQC/0.11.7-Java-1.8.0_74
# module load BBMap/38.00
# module load SAMtools/1.6-foss-2016b
# module load Bowtie/1.2.2-foss-2016b
# module load Trimmomatic/0.36-Java-1.8.0_92

suppressMessages(library(parallel))
suppressMessages(library(permseq))
suppressMessages(library(ShortRead))

args <- commandArgs(trailingOnly=TRUE)
assembly <- args[1]
tissue <- args[2]
histone <- args[3]


#27 was selected as the minimum size given that the shortest reads are this long
bowtieIndex <- paste0("/homes/users/mcoronado/scratch/5GenomesProject/DATA/genomeAssembly/",assembly,"/bowtieIndex_",assembly,"/",assembly)
bowtiebin <- "/homes/aplic/noarch/software/Bowtie/1.2.2-foss-2016b/bin/"

wDir = "/homes/users/mcoronado/scratch/5GenomesProject"

sampleData <- read.table(paste0("files/",tissue,"_",assembly,"_",histone,".txt"), stringsAsFactors=F)$V1

print(paste0(tissue,"-",assembly,"-",histone))
print(sampleData)
print(paste0("bowtieIndex: ", bowtieIndex))


permseqMapping <- mclapply(sampleData, function(x){
  readAllocate <- readAllocate(object=NULL,
  	chipFile=paste0(wDir, "/DATA/ChIPSeq_5genomes/",tissue,"FastP/",x,".clean.fastq"),
  	outfileLoc=paste0(wDir, "/ANALYSES/encodeChIPseqPipeline/mapping/", assembly,"/",tissue,"/",histone,"/",x,"_mapping"),
  	outputFormat="BED",
  	chipThres = 500,
  	bowtieIndex=bowtieIndex,
  	csemDir="/aplic/noarch/software/csem/2.3-Perl-5.26.1-foss-2018a/",
  	bowtieDir=bowtiebin,
  	pBowtie = 30)
}, mc.cores = 25)

