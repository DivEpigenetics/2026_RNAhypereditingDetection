# 2026_RNAhypereditingDetection

## RNA Hyperediting Detection Pipeline README file
This repository contains the implementation of the RNA hyperediting detection workflow used in our study of A-to-I RNA hyperediting in the marbled crayfish (Procambarus virginalis). The workflow is based on the RNA hyperediting detection strategy developed by Porath et al. (2014) and has been adapted for marbled crayfish. The purpose of this repository is to facilitate reproduction of our analyses and application of the workflow to additional datasets.

### Pipeline overview
The pipeline consists of six sequential steps:

    1. Genome alignment and extraction of unmapped reads
        ◦ Paired-end RNA-seq reads are aligned to the reference genome using BWA.
        ◦ Read pairs for which both mates remain unmapped are extracted for downstream analysis.
    2. Filtering of unmapped reads
        ◦ Low-quality reads and reads containing excessive ambiguous bases or low-complexity sequences are removed.
    3. Generation of nucleotide-converted reference genomes
        ◦ Twelve mismatch-masked reference genomes corresponding to all possible nucleotide substitution classes are generated.
        ◦ The substitution classes are defined in config/classes.tsv.
    4. Alignment to nucleotide-converted references
        ◦ Filtered reads are converted according to each substitution class and aligned independently against the corresponding converted reference genome.
    5. Detection of hyperedited reads, editing sites and editing clusters
        ◦ Alignments are evaluated following the criteria described by Porath et al. (2014)
        ◦ Hyperedited reads, individual editing sites and editing clusters are reported as BED files.
    6. Normalization
        ◦ The number of detected editing sites is normalized as:
            ▪ SPMM (Sites Per Million primary mapped reads)
            ▪ SPMB (Sites Per Million mapped bases)

### Repository structure
#### scripts/
    01_initial_map_and_unmapped.sh
    02_filter_unmapped.sh
    03_build_masked_refs.sh
    04_masked_align.sh
    05_call_and_cluster.sh
    05_calling_per_sites.py
    06_sites_per_million.sh

#### config/
    classes.tsv

### Requirements
    • Python 3
    • BWA (v0.7.18)
    • SAMtools (v1.20)
    • GNU awk
    • gzip

### Input
    The pipeline requires:
The pipeline requires:
        ◦ paired-end RNA-seq FASTQ files
        ◦ a reference genome in FASTA format

### Output
The pipeline produces:
        ◦ filtered unmapped reads
        ◦ nucleotide-converted reference genomes
        ◦ alignments against each substitution class
        ◦ BED files containing:
            ◦ hyperedited reads
            ◦ editing sites
            ◦ editing clusters
        ◦ a summary table containing the number of editing sites together with SPMM and SPMB values for each substitution class.

### Citation
If you use this repository, please cite both:
    1. Porath, H. T., Carmi, S. & Levanon, E. Y. A genome-wide map of hyper-edited RNA reveals numerous new sites. Nat. Commun. 5, 4726 (2014).
who originally developed the RNA hyperediting detection methodology.
    2. Diaz et al., which describes the implementation and application of this workflow.
