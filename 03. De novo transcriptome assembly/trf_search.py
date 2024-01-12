import sys
import os
from Bio import SeqIO
import subprocess

new_ref_tes = sys.argv[1]
transcript = sys.argv[2]

trf_path = "/lustre/home/ibe/mcoronado/scratch/5GenomesProject/ANALYSIS/chimerics/trf409.linux64"

# Search for single repeats in TEs
dicc_sr = {}
with open(os.path.basename(new_ref_tes)+".trf", 'r') as inFile:
    nbInSeq = 0
    for line in inFile:
        row = line.split(" ")
        if len(row) > 1 and "Sequence:" in row[0]:
            nbInSeq += 1
            seqName = row[1][:-1]
        if len(row) == 2 and not "Sequence:" in row[0]:
            start = row[0]
            end = row[1]
            if seqName in dicc_sr.keys():
                dicc_sr[seqName] += int(end) - int(start) + 1
            else:
                dicc_sr[seqName] = int(end) - int(start) + 1


for te in SeqIO.parse(new_ref_tes, "fasta"):
    if te.id in dicc_sr.keys():
        te.name = te.id.split("::")[0]
        te_len = len(str(te.seq))
        lenSR = dicc_sr[te.id]
        print(transcript + "\t" + te.name + "\t" + str((lenSR * 100) / te_len))
    else:
        te.name = te.id.split("::")[0]
        print(transcript + "\t" + te.name + "\t" + "0")

