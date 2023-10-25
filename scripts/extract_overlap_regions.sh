# extract overlap regions between TFBS & THS

for treat in "wt_a4_vs_un_peaks" "wt_kv_vs_un_peaks" "setiwt_e2_vs_mk_peaks";do
    regionBed="DARs/${treat}.narrowPeak.bed"
    for motiffile in ensembled_motif_mapping_files;do
        motif=${motiffile%.merged.bed}
        bedtools intersect -a $motiffile -b $regionBed -u -f 0.5 > ../Overlap_TFBS/${motif}_${treat}_overlap.bed
    done
done