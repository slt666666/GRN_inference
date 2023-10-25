# GRN_inference

0. (DEG analysis)

1. QC and alignment of ATAC-seq sequence reads (snakemake part) 

    Script:

    * ATAC-seq_read_QC_and_alignment.sh
 
2. Call differentially accessible regions (DARs) (Ding et al., 2021)

    Script:

    * call_DARs.sh

3. make position weight matrices (PWMs) for motif mapping (Kulkarni et al., 2019)

    Script:

    * meme2jaspar.py

4. PWM Mapping using promoter regions (Kulkarni et al., 2019)

    Script:

    * PWM_mapping.sh

5. extract TAIR10 promoter regions for hypergeometric test (Kulkarni et al., 2019)

    Script:

    * get_promoter_regions.sh

6. hypergeometric test to determine the overlap between TFBSs and THSs (Kulkarni et al., 2019)

    Script:

    * hypergeometric_test_overlap_TFBS_THS.sh

7. extract differentially expressed TFs and their target genes
    Script:
    * extract_DE_TFs_and_targets.py
    