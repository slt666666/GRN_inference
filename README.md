# GRN_inference

0. (DEG analysis)

1. QC and alignment of ATAC-seq sequence reads (snakemake part) 

    * Script    : ATAC-seq_read_QC_and_alignment.sh
    * Input     : fastq files of ATAC-seq
    * Output    : bam files
 
2. Call differentially accessible regions (DARs) (Ding et al., 2021)

    * Script    : call_DARs.sh
    * Input     : bam files generated by 1.
    * Output    : bed files of DARs (ex. wt_a4_vs_un_peaks.narrowPeak.bed)

3. make position weight matrices (PWMs) for motif mapping (Kulkarni et al., 2019)

    * Script    : meme2jaspar.py
    * Input     : meme format files from motif database (from https://meme-suite.org/meme/db/motifs)
    * Output    : meme format file & jaspar format file for each motif

4. PWM Mapping using promoter regions (Kulkarni et al., 2019)

    * Script    : PWM_mapping.sh
    * Input     : reference fasta file & motif PWM files generated by 3.
    * Output    : bed files of ensembled mapping result

5. extract TAIR10 promoter regions for hypergeometric test (Kulkarni et al., 2019)

    * Script    : get_promoter_regions.sh
    * Input     : gff & fasta files of reference genome
    * Output    : bed file of promoter region

6. hypergeometric test to determine the overlap between TFBSs and THSs (Kulkarni et al., 2019)

    * Script    : hypergeometric_test_overlap_TFBS_THS.sh
    * Input     : reference data, bed files generated by 4. & 5.
    * Output    : TF list with q_values

7. extract differentially expressed TFs and their target genes
    
    * Script: extract_DE_TFs_and_targets.py
    * Input : TF list with q_values & DEG results generated by 0. & gff file of reference genome
    * Output: differentially expressed TF list with their differentially expressed target genes 

    -> GRN is infered based on connections between TFs and target genes identified by 7. 
    