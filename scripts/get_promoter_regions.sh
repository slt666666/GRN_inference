# extract TAIR10 promoter regions for hypergeometric test (Kulkarni et al., 2019)
# 5000bp upstream & 1000bp downstream.
# input  : TAIR10 reference gff file
# output : promoter bed file

# extract promoter regions as bed format
python3 get_TAIR10_promoter_regions.py TAIR10_GFF3_genes.gff.gz
# sort & merge bed file
bedtools sort -i TAIR10_promoters.bed > TAIR10_promoters.sorted.bed
bedtools merge -i TAIR10_promoters.sorted.bed > TAIR10_promoters.sorted.merged.bed