<a name="readme-top"></a>
# Chimeric gene-TE transcripts in _D. melanogaster_ natural strains

Scripts and data for the paper "Transposons contribute to the diversification of the head, gut, and ovary transcriptomes across Drosophila natural strains".

## Directories
- **Chimeric gene-TE transcripts**: folder that contains a tsv file with the 2,169 chimeric transcripts reported
- **Reference-guided transcriptome assembly**: contains the script to perform the reference-guided transcriptome assembly following Pertea et al. 2016
- **Detection of chimeric genes: *De novo* transcriptome assembly**, **RepeatMasker** and **Minimap2**: contains the script to perform the de novo transcriptome assembly using Trinity, and run RepeatMasker on the *de novo* transcripts and align the transcript to the genome with Minimap2
- **Splice sites motifs**: scripts using the MEME suite to find the splice sites motif in D. melanogaster exon-intron junctions
- **Expression data**: script to generate expression matrix using trinity and salmon
- **ChIP-seq encode pipeline**: a sample script to run the encode pipeline (mapping and peak calling)
- ***Roo* analyses**: scripts to perform the blast analysis with the roo region
- **GO clustering**: script to parse GO clustering David tool output
- **CPAT**: script to run CPAT
- **pfamScan**: script to run PFAM
- **Scripts and data for figures**: contains RData and scripts to create the figures shown in the paper

## Citation
Coronado-Zamora, M. and González, J. Transposons contribute to the diversification of the head, gut, and ovary transcriptomes across Drosophila natural strains. BioRxiv.

## Contact

Project: [https://github.com/GonzalezLab/chimerics-transcripts-dmelanogaster](https://github.com/GonzalezLab/chimerics-transcripts-dmelanogaster)

<p align="right">(<a href="#readme-top">back to top</a>)</p>
