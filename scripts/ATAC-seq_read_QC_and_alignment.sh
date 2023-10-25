# ATAC-seq read QC & alignment
# input  : fastq files of ATAC-seq
# output : sorted bam files

## Quality check
fastqc --extract -f fastq -o sample_X_rep_X sample_X_rep_X.fastq

## Possibly Adapter removal like : https://informatics.fas.harvard.edu/atac-seq-guidelines.html

## Alignment by bowtie2
bowtie-build -f TAIR10_refrence.fasta TAIR10
bowtie2 --threads 8 -x TAIR10 --very-sensitive --end-to-end --maxins 20000 --no-discordant --no-mixed --fr --time --no-unal --qc-filter -1 sample_X_rep_X_R1.fastq -2 sample_X_rep_X_R2.fastq |
samtools view -b --reference TAIR10_refrence.fasta --threads 8 |
samtools sort -o sample_X_rep_X.sorted.bam && samtools index sample_X_rep_X.sorted.bam

## Put flag of duplications to bam file
picard MarkDuplicates I=sample_X_rep_X.sorted.bam O=sample_X_rep_X.sorted.mark_dup.bam M=sample_X_rep_X_markdup.txt REMOVE_SEQUENCING_DUPLICATES=false TAGGING_POLICY=OpticalOnly TAG_DUPLICATE_SET_MEMBERS=true
