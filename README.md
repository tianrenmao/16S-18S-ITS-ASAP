# Amplicon Sequence Analysis Pipeline (ASAP 1.4, May 2018)

This is an automatic pipeline for analysis of amplicon sequence data including 16S, 18S and ITS. It wraps QIIME commands and complements them with additional analysis where QIIME is not good at, such as combine multiple sequencing runs, OTU clustering and chimeric removal with UPARSE, alignment filtering with Gblock, removing Chloroplast sequences.

Below is the work flow:
- Search datasets (It accepts multiple sequencing runs)
- Read quality check - FastQC
- Merge pair-end sequence files - PEAR
- Demultiplexing (split library) - QIIME (1.9.1)
- Trim forward and reverse primer and any extra bases
- Merge sequence of different datasets
- Merge map files
- Dereplicate reads and calculate abundance - VSEARCH
- OTU clustering, chimera and singleton removal - UPARSE
- Make OTU table - VSEARCH
- Taxonomic classification of OTU â€“ RDP Classifier
- Remove Chloroplast and Mitochondria from OTU table
- Alignment and filtering â€“ MAFFT & Gblock
- Phylogenetic tree construction â€“ QIIME
- Convert and resample OTU table - QIIME
- Read summary
- Make taxonomic summary chart â€“ QIIME
- Make heat map (phylum to family) â€“ R script
- Alpha diversity - QIIME
- Beta diversity and plots - QIIME
- Summary

Running ASAP is quite easy and speedy. Users just need to put their sequence files and corresponding map files in different directories, copy the wrapper script and change parameters, and run the script. Wait for 1-10 hours and harvest the results. See the folder test_data for examples.

1. Put your fastq files (R1, R2 and I1, I1 is optional, see the settings in the wrapper script) and corresponding map file in dataset_1, dataset_2, ...
2. Copy the wrapper script and change parameters and configuration
3. nohup Perl ASAP_wrapper.pl &

The results includes OTUs sequences, OTU table, phylogentic tree, alpha diversity indexes (- OTU, shannon, simpson, dominance, goods coverage, PD whole), rarefaction curve of all indexes, PCoA charts (bray curtis, euclidean, unifrac), heat maps, bar charts of community composition at all levels etc.

After running, ASAP summarize statistics such as total read number, chimeric sequence number, singltons, average sequence number, OTUs number, Chloroplast OTUs, which can be put in manuscript easily. Method description is attached.

ASAP chooses the most commonly used and cited softwares. It has been testified using real data of our lab and give highly similar results to other pipelines.

Dependency:
- QIIME	1.9.1
- RDP Classifier	>=2.12
- FastQC	>=0.11.5
- MAFFT	>=3.8.3
- Perl	>=5.16

Installing QIIME virtual environment:
- conda create --yes -n qiime191 -c bioconda python=2.7 numpy=1.10 matplotlib=1.4.3 mock nose qiime
- source activate qiime191 (redo this if you deactivated the environment)

Contact Tian Renmao (rtian@ou.edu, University of Oklahoma) for questions.


