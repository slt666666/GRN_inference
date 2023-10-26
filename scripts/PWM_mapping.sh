# PWM mapping by Ensemble approach (Kulkarni et al., 2019)

# make PWM file of each TF for FIMO and Cluster-Buster
python3 meme2jaspar.py

# PWM maping to reference genome to identify TFBS by FIMO & Cluster-Buster
for each_motif in all_motifs; do
    fimo -o $each_motif ${each_motif}.meme TAIR10_promoters.sorted.merged.fasta
    cbust ${each_motif}.jaspar TAIR10_promoters.sorted.merged.fasta -c 0 -f 1 > ${each_motif}.txt
done

# ensemble motif mapping result of FIMO & CB
# (merge all match in FIMO & top 7000 match in CB)
python3 Ensemble_motif_mapping.py
for each_motif in all_motifs; do
    bedtools merge -i each_motif_ensembled.bed > each_motif_ensembled.merged.bed
done