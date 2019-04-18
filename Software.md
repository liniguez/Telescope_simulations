Runing Software
================
Luis P Iniguez
4/12/2019

This is how every counting method algorithm was run per simulation:

### Basics mapping

``` bash
SAMP=HML2_SIM_1 #needs to change for each simulation
conda activate RpETEtr
bowtie2 -f -p $(nproc) --no-unal -x bowtie2index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ${SAMP}/sample_01_1.fasta -2 ${SAMP}/sample_01_2.fasta -S ${SAMP}/default.sam
samtools view -bS ${SAMP}/default.sam > ${SAMP}/default.bam
samtools sort ${SAMP}/default.bam > ${SAMP}/default_sorted.bam

bowtie2 -f -p $(nproc) -k 100 --no-unal --very-sensitive-local --score-min L,0,1.6 -x bowtie2index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ${SAMP}/sample_01_1.fasta -2 ${SAMP}/sample_01_2.fasta -S ${SAMP}/k100.sam
samtools view -b ${SAMP}/k100.sam > ${SAMP}/k100.bam
samtools sort ${SAMP}/k100.bam > ${SAMP}/k100_sorted.bam
```

### RepEnrich

``` bash
python RepEnrich2/RepEnrich2_subset.py ${SAMP}/default.bam 30 ${SAMP}/RepEnrich_results --pairedend TRUE
python RepEnrich2/RepEnrich2.py --is_bed TRUE --cpus $(nproc) --pairedend TRUE transcripts_RepEnrich.bed  ${SAMP}/ RepEnrich_results repenrich_HML2 ${SAMP}/RepEnrich_results_multimap_R1.fastq --fastqfile2 ${SAMP}/RepEnrich_results_multimap_R2.fastq ${SAMP}/RepEnrich_results_unique.bam 
```

### TEtranscript

``` bash
TEtranscripts -t ${SAMP}/k100.bam -c ${SAMP}/k100.bam --GTF gencode.v30.basic.annotation.gtf --TE transcripts_TEtrasncripts.gtf --stranded no --project ${SAMP}/TEtranscript --verbose 1 -i 1000
grep -P '.+\:.+\:' ${SAMP}/TEtranscript.cntTable > ${SAMP}/TEtranscript_count.txt
rm ${SAMP}/TEtranscript.cntTable ${SAMP}/TEtranscriptDESeq.R

conda deactivate
```

### RSEM

``` bash
conda activate RSEMSLM
rsem-calculate-expression --bowtie2 --bowtie2-sensitivity-level very_sensitive --no-qualities --paired-end ${SAMP}/sample_01_1.fasta ${SAMP}/sample_01_2.fasta rsem_DB/reference ${SAMP}/rsem
```

### SalmonTE

``` bash
perl -e '{ open FILE, $ARGV[0]; while(<FILE>) {
 chomp $_;
 if ($_ =~ /^>(.+)/) {
  if($header ne "") {print "\@".$header."\n".$sequence."\n"."+"."\n".$sequence_quality."\n";
  }$header = $1;$sequence = "";$sequence_length = "";$sequence_quality = "";
 }else{$sequence .= $_;$sequence_length = length($_);
  for(my $i=0; $i<$sequence_length; $i++) {$sequence_quality .= "I"}}}close FILE;
  print "\@".$header."\n".$sequence."\n"."+"."\n".$sequence_quality."\n";}' ${SAMP}/sample_01_1.fasta > ${SAMP}/sample_01_1.fastq
perl -e '{ open FILE, $ARGV[0]; while(<FILE>) {
 chomp $_;
 if ($_ =~ /^>(.+)/) {
  if($header ne "") {print "\@".$header."\n".$sequence."\n"."+"."\n".$sequence_quality."\n";
  }$header = $1;$sequence = "";$sequence_length = "";$sequence_quality = "";
 }else{$sequence .= $_;$sequence_length = length($_);
  for(my $i=0; $i<$sequence_length; $i++) {$sequence_quality .= "I"}}}close FILE;
  print "\@".$header."\n".$sequence."\n"."+"."\n".$sequence_quality."\n";}' ${SAMP}/sample_01_2.fasta > ${SAMP}/sample_01_2.fastq

  
SalmonTE/SalmonTE.py quant --reference=salmonTE_DB ${SAMP}/sample_01_1.fastq ${SAMP}/sample_01_2.fastq --outpath=${SAMP}/salmonTE_all --exprtype=count
conda deactivate
```

### Telescope

``` bash
conda activate TELESCOPE
telescope assign --outdir ${SAMP}/ --reassign_mode exclude ${SAMP}/k100.bam transcripts.gtf
conda deactivate
```

###### Note:

All samples could be run at once with a for:

``` bash
for i in $(ls -d HML2_SIM *)
do
SAMP=$i

...
...
...

done
```
