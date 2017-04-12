This directory stores pre-processed data from the study.

Specifically the following datasets are included:
RNA-seq counts - Libraries were prepared from total RNA using Illumina's Stranded mRNA kit and sequenced on NextSeq 500 as paired 75 bp reads. The reads were aligned to both human and mouse genomes using subjunc. Aligned reads were then sorted by looking at the total number of matched base pairs for each pair of reads. In order for the read pair to be assigned as human, we required 20 more matched pairs in human alignment compared to the mouse alignment. Human reads were then summarised with featureCounts as the gene level to in-built RefSeq annotation. 
rna.seq.counts.unfiltered.rds - this file contains RNA-seq data from this study in the form of edgeR's DGEList object and contains count matrix, gene annotation and sample annotation.
combined.TCGA.normal.and.PDX.counts.rds - this file additionally contains summarised counts for RNA-seq data from TCGA. These were re-aligned using subread and obtained from the follwing GEO repository: GSE62944.

Copy number calls.

Somatic mutation calls.
