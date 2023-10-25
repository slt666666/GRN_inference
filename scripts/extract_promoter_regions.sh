# extract TAIR10 promoter regions (Kulkarni et al., 2019)
# input  : reference gff & fasta
# output : promoters bed & fasta

python3 get_TAIR10_promoter_regions.py TAIR10_GFF3_genes.gff.gz
bedtools merge -i TAIR10_promoters.bed > TAIR10_promoters.merged.bed
bedtools getfasta -fi TAIR10_chr_all.fas -bed TAIR10_promoters.merged.bed > TAIR10_promoters.fasta