#### Chimeric gene-TE file description

Column description:

- **Tissue**: head, gut, ovary - body part in which the chimeric transcript was detected
- **Assembly**: AKA-017, JUT-011, MUN-016, SLA-001, TOM-007 - strain in which the chimeric transcript was detected
- **Category**: exonOverlapTE, partialOverlap, TEoverlapExon: if the TE was detected within an exon, overlapping with the exon, or the exon was detected within a TE (not used for downstream analyses)
- **Transcript**: ID according to Trinity
- **Length transcript**: length of the assembled transcript in bp
- **stringtieID**: id of the transcript according to the reference-guided assembly. Most genes have a FBtr ID, the ones starting with MSTRG are transcripts not found in the reference annotation
- **FlyBaseRef**: id of the transcript according to FlyBase. 
- **FlyBaseGeneRef**: id of the gene in FBgn format
- **nExon**: number of exons detected in the transcript after performing the alignment to the genome with minimap2
- **coordTE**: coordinates where the TE was detected
- **lengthExon**: coordinates of the exon that has a matching TE
- **exonPos**: position in which the TE was detected
- **TEconsensusID**: name of the TE according to the MCTE library
- **TEfamily**: family of the TE
- **classCode**: transcript type code according to the reference-guided assembly
- **definition**: transcript type description according to the reference-guided assembly
- **exonPosDescription**: position in which the TE is detected: first, last or middle exon (corresponds to the 3'/5' UTR, internal exon classification) or (TE inside) when is detected within
-**distanceIntron**: for some transcripts, distances in bp between the first or the last exon with the next one (not used in downstream analyses)
-**FIMOenrichment**: motif and p-value if a motif was found with FIMO
-**SS**: a code specifing the meaning of FIMOenrichment column
-**lengthTE**: length of the TE insertion in bp
-**score**: score of the TE according to RepeatMasker (not used for downstream analyses)
-**AS**: logical column according if the AS sites where detected or not according to `FIMOenrichment` and `SS` columns
-**group**: 1: overlap and AS insertions group, 2: internal insertions group
