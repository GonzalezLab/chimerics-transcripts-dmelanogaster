#### Chimeric gene-TE file description

Column description:

- **Tissue**: head, gut, ovary - body part in which the chimeric transcript was detected
- **Strain**: AKA-017, JUT-011, MUN-016, SLA-001, TOM-007 - strain in which the chimeric transcript was detected
- **Transcript**: id of the transcript according to FlyBase. 
- **Transcript Trinity**: ID according to Trinity
- **TEfamily**: family of the TE
- **TE coordinate**: coordinates where the TE was detected in the assembled region
- **Transcript Stringtie**: id of the transcript according to the reference-guided assembly. Most genes have a FBtr ID, the ones starting with MSTRG are transcripts not found in the reference annotation
- **Class**: transcript type code according to the reference-guided assembly
- **Gene**: id of the gene in FBgn format
- **Length transcript**: length of the assembled transcript in bp
- **Total number of exons**: number of exons detected in the transcript after performing the alignment to the genome with minimap2
- **Expon position**: position in which the TE was detected (number)
- **Exon position description**: position in which the TE is detected: first, last or middle exon (corresponds to the 3'/5' UTR, internal exon classification) or (TE inside) when is detected within
- **Type**: exon_within_TE, exon_within_TE_3, exon_within_TE_5, TE_overlap_3, TE_overlap_5, TE_within_exon
- **TE consensus**: name of the TE according to the MCTE library
- **TE superfamily**: superfamily of the TE
- **TE order**: order of the TE
- **TE class**: class of the TE
- **length TE incorporated**: length of the TE fragment detected in the chimeric transcipts in bp
- **Length TE total**: length of the TE insertion in bp
- **CP**: coding potential according to CPAT
- **status**: type of transcript according to FlyBase annotations (RNA,ncRNA,pseudogene,pre_miRNA,tRNA,snoRNA)
- **avgExpr**: average level of expression (in TMM, averaging the expression level of three replicates)
- **SS**: a code specifing the meaning of FIMOenrichment column
- **FIMOenrichment**: motif and p-value if a motif was found with FIMO
- **group**: 1: overlap and AS insertions group, 2: internal insertions group
- **score**: score of the TE according to RepeatMasker (not used for downstream analyses)
- **position of roo**: position in the roo in which there's a match (low complexity region)
- **roo type**: match of the roo fragment in the consensus (low complexity region [repeat], or LTR)