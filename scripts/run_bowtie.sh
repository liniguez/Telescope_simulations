#!/bin/sh

cd /home/ubuntu/efs/simulations/telescope_june_2018
Rscript get_fastq_sim.R

for i in $(find -type d -name 'HML2_SIM*')
do 
 perl fasta_to_fastq.pl ${i}/sample_01_1.fasta > ${i}/sample_01_1.fastq &
 perl fasta_to_fastq.pl ${i}/sample_01_2.fasta > ${i}/sample_01_2.fastq &
 wait
 bowtie2 --no-unal --score-min L,0,1.6 -p $(nproc) -k 100 --very-sensitive-local -x /home/ubuntu/references/bowtie2/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ${i}/sample_01_1.fastq -2 ${i}/sample_01_2.fastq > ${i}/mapping_k100.sam
 telescope assign --tempdir ${i} --quiet --outdir ${i} --updated_sam ${i}/mapping_k100.sam transcripts_MB_hg38.gtf
done

