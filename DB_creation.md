Data Bases
================
Luis P Iniguez
4/10/2019

All the annotations used for creating all databases can be found [here](https://github.com/mlbendall/telescope_annotation_db/tree/master/builds/HERV_rmsk.hg38.v2) under transcripts.gtf version hg38.v2

Installing and getting all the algorithms for Transposable Elements quantification:

``` bash
git clone https://github.com/nerettilab/RepEnrich2.git
git clone https://github.com/hyunhwaj/SalmonTE
conda create -n TeleSim python=2.7.12 bowtie2 bedtools=2.25.0 samtools=1.3.1 biopython=1.66 pysam=0.14.1 tetoolkit rsem docopt snakemake pandas r-tidyverse r-scales r-writexls r-biocmanager
conda create -n TELESCOPE python=3.6 numpy=1.13 scipy=0.19.0 pysam=0.12 cython intervaltree samtools=1.5 openssl=1.0.2o intervaltree=2.1.0 pyyaml=3.12
```

Downloading necessary files for data bases:

``` bash
conda activate TeleSim
#Human Genome
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa

#Human Gene Annotations
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.basic.annotation.gtf.gz
gunzip gencode.v30.basic.annotation.gtf.gz
#Human TE Annotations
wget https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/HERV_rmsk.hg38.v2/transcripts.gtf
```

RepEnrich Database:

``` bash
grep -v '#' transcripts.gtf | grep -P '\texon\t' | awk '{split($10,a,"\"");print $1,$4,$5,a[2],"LTR","HERV";}' OFS='\t' > transcripts_RepEnrich.bed
python RepEnrich2/RepEnrich2_setup.py transcripts_RepEnrich.bed hg38.fa repenrich_DB/ --is_bed TRUE --threads $(nproc)
```

TEtranscripts gtf:

``` bash
cp transcripts.gtf transcripts_TEtrasncripts.gtf
sed -i 's/repName/family_id/' transcripts_TEtrasncripts.gtf 
sed -i 's/intModel/family_id/' transcripts_TEtrasncripts.gtf
sed -i 's/category/class_id/' transcripts_TEtrasncripts.gtf
sed -i 's/repClass/class_id/' transcripts_TEtrasncripts.gtf
```

RSEM Database:

``` bash
mkdir rsem_DB
rsem-prepare-reference -gtf transcripts.gtf --bowtie2 hg38.fa rsem_DB/reference
```

SalmonTE Database:

``` bash
mkdir salmonTE_DB
awk '{if($0~/^>/){print $1,"HERV","Homo sapiens"}else print $0;}' OFS="\t" rsem_DB/reference.transcripts.fa > salmonTE_DB/reference_transcripts.fa
echo "HERV,HERV,Endogenous Retrovirus,Transposable Element" >> SalmonTE/reference/clades_extended.csv
SalmonTE/SalmonTE.py index --input_fasta=salmonTE_DB/reference_transcripts.fa --ref_name=salmonTE_DB
```

``` bash
conda deactivate
```
