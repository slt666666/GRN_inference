import pandas as pd
import numpy as np
import sys

def extract_forward(chr_gff, i, distance, feature):
    start = chr_gff.iloc[i, 3]
    end = chr_gff.iloc[i, 4]
    strand = chr_gff.iloc[i, 6]
    gene_id = chr_gff.iloc[i, 8][3:12]
    if i == 0:
        if start > distance:
            region = [each_chr, start-distance, start-1, gene_id, feature]
        else:
            region = [each_chr, 1, start-1, gene_id, feature]
    else:
        if chr_gff[(chr_gff.iloc[:, 3] < start) & (chr_gff.iloc[:, 4] > start)].shape[0] > 0:
            region = ["NA", "NA", "NA"]
        else:
            end_prev = chr_gff[chr_gff.iloc[:, 4] < start].iloc[:, 4].max()
            if ((start-1) - (end_prev+1)) >= distance:
                region = [each_chr, start-distance, start-1, gene_id, feature]
            else:
                region = [each_chr, end_prev+1, start-1, gene_id, feature]
    if region[1] >= region[2]:
        region = "NA"
    return region

def extract_backward(chr_gff, chr_end, i, distance, feature):
    start = chr_gff.iloc[i, 3]
    end = chr_gff.iloc[i, 4]
    strand = chr_gff.iloc[i, 6]
    gene_id = chr_gff.iloc[i, 8][3:12]
    if i == chr_gff.shape[0]-1:
        if chr_end - end > distance:
            region = [each_chr, end+1, end+distance, gene_id, feature]
        else:
            region = [each_chr, end+1, chr_end, gene_id, feature]
    else:
        if chr_gff[(chr_gff.iloc[:, 3] < end) & (chr_gff.iloc[:, 4] > end)].shape[0] > 0:
            region = ["NA", "NA", "NA"]
        else:
            start_next = chr_gff[chr_gff.iloc[:, 3] > end].iloc[:, 3].min()
            if ((start_next-1) - (end+1)) >= distance:
                 region = [each_chr, end+1, end+distance, gene_id, feature]
            else:
                region = [each_chr, end+1, start_next-1, gene_id, feature]
    if region[1] >= region[2]:
        region = "NA"
    return region


gff_file = sys.argv[0]
gff = pd.read_csv(gff_file, sep="\t", comment="#", header=None)
chr_len = gff[gff.iloc[:, 2] == "chromosome"]
gff = gff[gff.iloc[:, 2] == "gene"]

promoter = []
for each_chr in ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]:
    chr_gff = gff[gff.iloc[:, 0] == each_chr]
    chr_gff = chr_gff.sort_values(4)
    chr_gff = chr_gff.reset_index(drop=True)
    chr_end = chr_len[chr_len.iloc[:, 0] == each_chr].iloc[:, 4].values[0]
    for i in range(chr_gff.shape[0]):
        strand = chr_gff.iloc[i, 6]
        if strand == "+":
            region = extract_forward(chr_gff, i, 5000, "upstream")
            if region != "NA": promoter.append(region)
            region = extract_backward(chr_gff, chr_end, i, 1000, "downstream")
            if region != "NA": promoter.append(region)
        else:
            region = extract_backward(chr_gff, chr_end, i, 5000, "upstream")
            if region != "NA": promoter.append(region)
            region = extract_forward(chr_gff, i, 1000, "downstream")
            if region != "NA": promoter.append(region)
promoter = pd.DataFrame(promoter)
# display(promoter)
promoter.to_csv("TAIR10_promoters.bed", sep="\t", index=None, header=None)