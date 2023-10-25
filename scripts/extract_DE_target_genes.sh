# extract differentially expressed TFs and their target genes

# extract overlap regions between TFBS & THS for identification of target genes of TFs
extract_overlap_regions.sh

# extract differentially expressed TFs & target genes based on RNA-seq result
python3 extract_DE_TFs_and_targets.py