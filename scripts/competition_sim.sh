#!/bin/sh

COND=$1


SA1=/home/ubuntu/efs/simulations/telescope_june_2018/${COND}/sample_01_1.fasta
SA2=/home/ubuntu/efs/simulations/telescope_june_2018/${COND}/sample_01_2.fasta
OUTBA=/home/ubuntu/efs/simulations/telescope_june_2018/${COND}/mapping_k100
BT1=/home/ubuntu/efs/simulations/annotation/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set
BT2=/home/ubuntu/references/hg38/bowtie2/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index
GTF=/home/ubuntu/efs/simulations/annotation/transcripts_MB_hg38.gtf
GTFGE=/home/ubuntu/efs/simulations/annotation/gencode.v28.basic.annotation.gtf

#Best Alignment#
bowtie2 -p $(nproc) --very-sensitive-local -S ${OUTBA}_best.sam --score-min L,0,1.6 -x 1200 -f -x ${BT2} -1 ${SA1} -2 ${SA2} 
samtools view -@ $(nproc) -b ${OUTBA}_best.sam > ${OUTBA}_best.bam
samtools sort -@ $(nproc) ${OUTBA}_best.bam > ${OUTBA}_best_sorted.bam
rm ${OUTBA}_best.bam ${OUTBA}_best.sam
#K100 Alignment#
#bowtie2 -p $(nproc) -k 100 -S ${OUTBA}_k100.sam --very-sensitive-local --score-min L,0,1.6 -x 1200 -f -x ${BT2} -1 ${SA1} -2 ${SA2} 
samtools view -@ $(nproc) -b ${OUTBA}.sam > ${OUTBA}.bam 
#rm ${OUTBA}_k100.sam

#Counting#
#Best#
htseq-count --format=bam --order=pos --stranded=no --idattr=locus --minaqual=0 ${OUTBA}_best_sorted.bam $GTF > ${OUTBA}_best_count.txt
#All#
htseq-count --format=bam --order=name --stranded=no --idattr=locus --minaqual=0 ${OUTBA}.bam $GTF > ${OUTBA}_all_count.txt
#Primary#
samtools view -@ $(nproc) -F 0x100 ${OUTBA}.bam | htseq-count --format=sam --order=pos --stranded=no --idattr=locus --minaqual=0 - $GTF > ${OUTBA}_primary_count.txt
#Unique#
samtools view -@ $(nproc) -f 0x100 ${OUTBA}.bam | cut -f1 | sort -u > ${OUTBA}_temp_reads.txt
samtools view -@ $(nproc) ${OUTBA}.bam > ${OUTBA}_temp.sam
perl -e '{ open (IN,"$ARGV[0]"); while ($l=<IN>){chomp $l; $h{$l}=1;}close(IN);
 open(SAM,"$ARGV[1]"); while(<SAM>){@vec=split("\t",$_); if(!$h{$vec[0]}){print $_;}}}' ${OUTBA}_temp_reads.txt ${OUTBA}_temp.sam | htseq-count --format=sam --order=pos --stranded=no --idattr=locus --minaqual=0 - ${GTF} > ${OUTBA}_unique_count.txt
rm ${OUTBA}_temp_reads.txt ${OUTBA}_temp.sam

#TEtranscript

cp ${GTF} ${GTF}_4tetranscripts.gtf
sed -i 's/repName/family_id/' ${GTF}_4tetranscripts.gtf 
sed -i 's/intModel/family_id/' ${GTF}_4tetranscripts.gtf
sed -i 's/category/class_id/' ${GTF}_4tetranscripts.gtf
sed -i 's/repClass/class_id/' ${GTF}_4tetranscripts.gtf
grep -v '#' ${GTF}_4tetranscripts.gtf > ${GTF}_4tetranscripts2.gtf

TEtranscripts -t ${OUTBA}.bam -c ${OUTBA}.bam --GTF ${GTFGE} --TE ${GTF}_4tetranscripts2.gtf --stranded no --project ${OUTBA}_TEtranscript --verbose 1 -i 1000
grep -P '.+\:.+\:' ${OUTBA}_TEtranscript.cntTable > ${OUTBA}_TEtranscript_count.txt

rm ${GTF}_4tetranscripts.gtf ${GTF}_4tetranscripts2.gtf ${OUTBA}_TEtranscript.cntTable ${OUTBA}_TEtranscript_DESeq.R




#RepEnrich#
bowtie --threads $(nproc) --time -f -m 1 --sam -X 600 --max ${OUTBA}_multimap.txt --chunkmbs 512 ${BT1} -1 ${SA1} -2 ${SA2} ${OUTBA}_uniq_REPENRICH.sam
samtools view -u ${OUTBA}_uniq_REPENRICH.sam | samtools sort > ${OUTBA}_uniq_REPENRICH_sorted.bam
samtools index ${OUTBA}_uniq_REPENRICH_sorted.bam

perl -e '{ open (FA, $ARGV[0]); while($h=<FA>){$seq=<FA>; chomp $h; chomp $seq; @vech=split("",$h);$vech[0]="@";$h=join("",@vech); @vecs=split("",$seq);$size=scalar(@vecs); @qual=("I")x$size; $qual=join("",@qual);print "$h\n$seq\n+\n$qual\n";}}' ${OUTBA}_multimap_1.txt > ${OUTBA}_multimap_1.fastq
perl -e '{ open (FA, $ARGV[0]); while($h=<FA>){$seq=<FA>; chomp $h; chomp $seq; @vech=split("",$h);$vech[0]="@";$h=join("",@vech); @vecs=split("",$seq);$size=scalar(@vecs); @qual=("I")x$size; $qual=join("",@qual);print "$h\n$seq\n+\n$qual\n";}}' ${OUTBA}_multimap_2.txt > ${OUTBA}_multimap_2.fastq

RREF=/home/ubuntu/efs/repenrich/hg38_repeatmasker_clean.txt
WHA=hg38_repeatmasker/
WHA2=hg38_repeatmasker
python /home/ubuntu/efs/repenrich/RepEnrich.py --cpus $(nproc) --pairedend TRUE ${RREF} ${OUTBA}_repenrich_${WHA} ${WHA2} /home/ubuntu/efs/repenrich/RepEnrich_setup_hg38 ${OUTBA}_multimap_1.fastq --fastqfile2 ${OUTBA}_multimap_2.fastq ${OUTBA}_uniq_REPENRICH_sorted.bam &


RREF=/home/ubuntu/efs/repenrich/HML2_mod.repenrich.txt
WHA=HML2/
WHA2=HML2_hg38
python /home/ubuntu/efs/repenrich/RepEnrich.py --is_bed TRUE --cpus $(nproc) --pairedend TRUE ${RREF} ${OUTBA}_repenrich_${WHA} ${WHA2} /home/ubuntu/efs/repenrich/HML2.hg38_2 ${OUTBA}_multimap_1.fastq --fastqfile2 ${OUTBA}_multimap_2.fastq ${OUTBA}_uniq_REPENRICH_sorted.bam &

#RREF=/home/ubuntu/efs/repenrich/transcripts_MB_hg38.repenrich.txt
#WHA=HERV/
#WHA2=HERV_hg38
#python /home/ubuntu/efs/repenrich/RepEnrich.py --is_bed TRUE --cpus $(nproc) --pairedend TRUE ${RREF} ${OUTBA}_repenrich_${WHA} ${WHA2} /home/ubuntu/efs/repenrich/HERV_hg38 ${OUTBA}_multimap_1.fastq --fastqfile2 ${OUTBA}_multimap_2.fastq ${OUTBA}_uniq_REPENRICH_sorted.bam &


wait

rm ${OUTBA}_uniq_REPENRICH.sam ${OUTBA}_uniq_REPENRICH_sorted.bam ${OUTBA}_uniq_REPENRICH_sorted.bam.bai ${OUTBA}_multimap_2.txt ${OUTBA}_multimap_2.fastq ${OUTBA}_multimap_1.txt ${OUTBA}_multimap_1.fastq
