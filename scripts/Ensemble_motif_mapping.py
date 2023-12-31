import pandas as pd
from pandas.errors import EmptyDataError
import numpy as np
import glob

def get_cb_top7000(file_path):
    cb_scores = []
    with open(file_path, "r") as f:
        lines = f.readlines()
        for i in lines:
            if i[:22] == "Cluster-Buster version":
                break
            else:
                if i[0] == ">":
                    region = i[1:-1]
                    chrom = i[1:5]
                    start = int(i[6:i.find(" ")].split("-")[0])
                elif (i[0] == "#") | (i[0] == "\n"):
                    pass
                else:
                    scores = i[:-1].split("\t")
                    cb_scores.append([chrom, start+int(scores[1])-1, start+int(scores[2])-1, float(scores[3])])

    cb_scores = pd.DataFrame(cb_scores)
    cb_scores.columns = ["chr", "start", "end", "score"]
    cb_scores = cb_scores.sort_values(by="score", ascending=False).iloc[:7000, :]
    cb_scores["tool"] = "CB"
    return cb_scores

def get_fimo_list(file_path):
    try:
        fimo = pd.read_csv(file_path, sep="\t", comment="#")
        fimo_scores = []
        for i in fimo.itertuples():
            seq_name = i[3]
            chrom = seq_name[:seq_name.find(":")]
            start = int(seq_name[seq_name.find(":")+1:seq_name.find("-")])
            fimo_scores.append([chrom, start+int(i[4])-1, start+int(i[5])-1, float(i[7])])
        fimo_scores = pd.DataFrame(fimo_scores)
        fimo_scores.columns = ["chr", "start", "end", "score"]
        fimo_scores["tool"] = "FIMO"
    except EmptyDataError:
        fimo_scores = pd.DataFrame()
        
    return fimo_scores

motifs = [i[i.rfind("/")+1:] for i in glob.glob("../data/PWM_mapping/FIMO/*")]
for motif in motifs:
    fimo_list = get_fimo_list("../data/PWM_mapping/FIMO/{}/fimo.tsv".format(motif))
    cb_top7000 = get_cb_top7000("../data/PWM_mapping/CB/{}.txt".format(motif))
    if fimo_list.shape[0] != 0:
        ensemble = pd.concat([fimo_list, cb_top7000]).reset_index(drop=True)
    else:
        ensemble = cb_top7000
    ensemble.to_csv("../data/PWM_mapping/Ensemble/{}.bed".format(motif), index=None, header=None, sep="\t")