# make PWM files for FIMO and Cluster-Buster

import sys
import re
import glob
import copy
import numpy as np

def read_meme_file(filepath):
	pwms={}
	data = open(filepath).read()

	pos=0
	while(1):
		rec_loc=data.find('\nMOTIF', pos)
		if rec_loc<0:
			break
		nl=data.find('\n', rec_loc+1)
		motif=data[rec_loc+6:nl].strip().split(':')[0]
		mat_header_start=data.find('letter-probability', rec_loc)
		mat_header_end=data.find('\n', mat_header_start+1)

		i = data.find(':', mat_header_start)
		kp = data[i+1:mat_header_end]
		r = re.compile(r"(\w+)=\s?([\w.\+]+)\s?")
		attribs=dict(r.findall(kp))

		ist=mat_header_end+1
		iend=-1

		j=0
		z=[]
		while(j<int(attribs["w"])):
			iend = data.find('\n', ist)
			z.append(list(map(float, data[ist:iend].strip().split())))
			ist = iend+1
			j+=1

		pwms[motif]=np.array(z)
		pos=ist
	return pwms

meme_files = glob.glob("../data/motif_database/source_data/*.meme")
for meme_file in meme_files:
    pwms=read_meme_file(meme_file)
    for k, p in pwms.items():
        output_file = "../data/motif_database/jaspar_format/"+k[:k.find(" ")]+".jaspar"
        with open(output_file, "a") as f:
            f.write(">"+k[:k.find(" ")]+"\n")
            np.savetxt(f, p.T, delimiter='\t', fmt='%0.5f')

meme_files = glob.glob("../data/motif_database/source_data/*.meme")
for meme_file in meme_files:
    with open(meme_file, "r") as f:
        lines = f.readlines()
        base_ind = []
        for i, k in enumerate(lines):
            if k[:5] == "MOTIF":
                break
            else:
                base_ind.append(i)
                
        Motif_ind = []
        for i, k in enumerate(lines):
            if k[:5] == "MOTIF":
                Motif_ind.append(i)
        Motif_ind.append(len(lines)-1)
        
        for ind in range(len(Motif_ind)-1):
            # range = Motif_ind[ind], Motif_ind[ind+1]-1
            motif_name = lines[Motif_ind[ind]].split(" ")[1]
            output_file = "../data/motif_database/meme_format/"+motif_name+".meme"
            with open(output_file, "a") as f:
                each_motif_ind = copy.copy(base_ind)
                each_motif_ind.extend(list(range(Motif_ind[ind], Motif_ind[ind+1])))
                for i in each_motif_ind:
                    f.write(lines[i])