import pandas as pd
import numpy as np
import sys


# calculate q values from p values
def calculate_q(p_seq):
    p_arr = np.asarray(p_seq)
    N = len(p_arr)
    idx_arr = np.argsort(p_arr)
    q_arr = p_arr[idx_arr] * N / (np.arange(N) + 1)
    q_arr = np.minimum.accumulate(q_arr[::-1])[::-1]
    q_arr[idx_arr] = q_arr.copy()
    return q_arr

# extract TFs that showed significant enrichment
def get_enrichment_TFBS(file, Motif_ID_list, q_threshold):
    TFBS = pd.read_csv(file, sep="\t", header=None)
    TFBS.columns = ["Motif_ID", "treat", "realNum", "nShuffle", "biggerNum", "p"]
    TFBS["q"] = calculate_q(TFBS["p"])
    TFBS = TFBS[[i in Motif_ID_list["Motif_ID"].values for i in TFBS["Motif_ID"]]]
    TFBS["AtID"] = [Motif_ID_list[Motif_ID_list["Motif_ID"] == i].AtID.values[0] for i in TFBS["Motif_ID"]]
    TFBS["Symbol"] = [Motif_ID_list[Motif_ID_list["Motif_ID"] == i].Motif_Symbol.values[0] for i in TFBS["Motif_ID"]]    
    TFBS["Source"] = [Motif_ID_list[Motif_ID_list["Motif_ID"] == i].Source.values[0] for i in TFBS["Motif_ID"]]    
    return TFBS[TFBS["q"] < q_threshold].reset_index(drop=True)


# Motif ID vs AtID information
Motif_ID_list = pd.read_csv("../data/motif_database/motif_IDs.csv", index_col=0)
TFBS = get_enrichment_TFBS(sys.argv[0], Motif_ID_list, 0.01)
TFBS.to_csv(sys.argv[1])