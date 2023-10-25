import pandas as pd
import numpy as np

# extract differentially expressed genes using RNA-seq data
def get_DEGs(file, logFC_threshold, FDR_threshold):
    DEGs = pd.read_csv(file, index_col=0)
    DEGs = DEGs[(DEGs["logFC"] > logFC_threshold) & (DEGs["FDR"] < FDR_threshold)]
    return DEGs.index.values

# exract target genes of TFs based on distance
def get_closest_genes(Motif_id, treat, DEGs, gene_gff):
    TFBS = pd.read_csv("../data/PWM_mapping/Overlap_TFBS/{}_{}_overlap.bed".format(Motif_id, treat), sep="\t", header=None)
    gene_gff["AtID"] = gene_gff.iloc[:, 8].str[3:12]
    target_genes = []
    for each_chr in ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]:
        chr_gff = gene_gff[gene_gff.iloc[:, 0] == each_chr]
        chr_TFBS = TFBS[TFBS.iloc[:, 0] == each_chr]
        for i in range(chr_TFBS.shape[0]):
            start = chr_TFBS.iloc[i, 1]
            end = chr_TFBS.iloc[i, 2]
            closest_forward = chr_gff[chr_gff.iloc[:, 4] < start].iloc[:, 4].max()
            closest_backward = chr_gff[chr_gff.iloc[:, 3] > end].iloc[:, 3].min()
            if start - closest_forward < closest_backward - end:
                target_genes.append(chr_gff[chr_gff.iloc[:, 4] == closest_forward].AtID.values[0])
            elif start - closest_forward > closest_backward - end:
                target_genes.append(chr_gff[chr_gff.iloc[:, 3] == closest_backward].AtID.values[0])
            elif start - closest_forward == closest_backward - end:
                target_genes.append(chr_gff[chr_gff.iloc[:, 4] == closest_forward].AtID.values[0])
                target_genes.append(chr_gff[chr_gff.iloc[:, 3] == closest_backward].AtID.values[0])
    target_genes = np.unique(list(set(target_genes) & set(DEGs)))
    return target_genes

# Gene annotation information
gene_gff = pd.read_csv("../data/TAIR10_GFF3_genes.gff", sep="\t", comment="#", header=None)
gene_gff = gene_gff[gene_gff.iloc[:, 2] == "gene"]

# get differentially expressed genes
DEGs = get_DEGs(files[condition][0], 1, 0.01)
# get TFs that showed TFBS enrichment
TFBS = pd.read_csv("~~~~q-value.csv")
TFBS = TFBS[TFBS["q"] < 0.01].reset_index(drop=True)
# extract differentially expressed TFs
DE_TFBS = TFBS[TFBS["AtID"].isin(DEGs)].reset_index(drop=True)
print("Differentially expressed TFs that showed TFBS enrichment in {}:".format(condition), len(DE_TFBS.AtID.unique()))
# annotate target genes of TFs based on genetic distance
DE_target_genes = [get_closest_genes(Motif_id, files[condition][2], DEGs, gene_gff) for Motif_id in DE_TFBS.Motif_ID]
DE_TFBS["DE_target_genes"] = [",".join(DE_target_gene) for DE_target_gene in DE_target_genes]
# save
DE_TFBS.to_csv("../data/{}_TF_and_target.csv".format(condition))