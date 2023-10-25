# input data  : aligment bam files for treatment
# output data : macs3 outputs
# ex) wt_a4_rep1.sorted.bam, wt_a4_rep2.sorted.bam vs wt_un_rep1.sorted.bam, wt_un_rep2.sorted.bam

# call accessible chromatin regions that more enriched in treatA compared to treatB
macs3 callpeak -t treatA_rep1.bam treatA_rep2.bam -c treatB_rep1.bam treatB_rep2.bam -f BAM \
               -n treatA_vs_treatB -g 135000000 -B -q 0.05
