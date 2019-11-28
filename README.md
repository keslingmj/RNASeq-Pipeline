# RNASeq Pipeline

This pipeline has 2 files.  The [pipeline.sh](pipeline.sh) is a bash script that
1. Maps the reads in the samples' fastq files onto the genome using hisat2
2. Using samtools sorts the bam files and re-saves them as bam files
3. Indexes the bam files using samtools for viewing reads in the IGV viewer
4. Creates QA files using samtools flagstat
5. Assembles transcript structures and quantitates these levels using stringtie but focuses on known genes
6. Stringtie expands transcript set to those having experimental evidence but not in the RefSeq set.  It is this set that is used for Salmon.
7. All transcripts, including new ones, are quantitated.  This is only for comparing to Salmon.
8. The fasta file of all known and newly established transcript structures is created with gffread
9. Salmon indexes these transcript sequences
10. salmon quant is run to estimate transcript-level TPM levels for each transcript-sample pair
11. [post-process-salmon.R](post-process-salmon.R) is run to gather TPM levels across all salmon quant.sf files and sum across transcripts of the same gene.

This pipeline is a modified version of one written by the Perteas in the Salzberg lab.
In order to run the pipeline, it is assumed that there exists the following directory structure containing the following files.