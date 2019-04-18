Data Bases
================
Luis P Iniguez
4/10/2019

Installing and getting all the algorithms for Transposable Elements quantification:

``` bash
git clone https://github.com/nerettilab/RepEnrich2.git
git clone https://github.com/hyunhwaj/SalmonTE
conda create -n RpETEtr python=2.7.12 bowtie2 bedtools=2.25.0 samtools=1.3.1 biopython=1.66 pysam=0.14.1 tetoolkit
conda create -n RSEMSLM rsem docopt snakemake pandas r-tidyverse r-scales r-writexls r-biocmanager 
conda create -n TELESCOPE python=3.6 numpy=1.13 scipy=0.19.0 pysam=0.12 cython intervaltree samtools=1.5 openssl=1.0.2o intervaltree=2.1.0 pyyaml=3.12
conda activate TELESCOPE
pip install git+git://github.com/mlbendall/telescope.git
conda deactivate
```

Downloading necessary files for data bases. All the annotations used for creating databases can be found [here](https://github.com/mlbendall/telescope_annotation_db/tree/master/builds/HERV_rmsk.hg38.v2) under transcripts.gtf version hg38.v2

``` bash
#Human Genome
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa

#Human Gene Annotations
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.basic.annotation.gtf.gz
gunzip gencode.v30.basic.annotation.gtf.gz

#Bowtie2 index
mkdir bowtie2index
cd bowtie2index
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
tar -xvzf GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
rm GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
cd ..

#Human TE Annotations and HML2 only in principal chromosomes
wget https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/HERV_rmsk.hg38.v2/transcripts.gtf
grep -v '#' transcripts.gtf | awk '{ if ($1 !~ "_")print$0;}' OFS='\t' |grep 'HML2' > HML2.gtf
```

RepEnrich Database: RepEnrich is ment to track the expression of whole families but not single locus resulution. Therfore a "hacked" version of the RepEnrich database was build here. In this version each HML2 is take into account as a single family.

``` bash
conda activate RpETEtr
grep -v '#' HML2.gtf | grep -P '\texon\t' | awk '{split($10,a,"\"");print $1,$4,$5,a[2],"LTR","HERV";}' OFS='\t' > transcripts_RepEnrich.bed
python RepEnrich2/RepEnrich2_setup.py transcripts_RepEnrich.bed hg38.fa repenrich_HML2/ --is_bed TRUE --threads $(nproc)
conda deactivate
```

TEtranscripts gtf:

``` bash
grep -v '#' transcripts.gtf > transcripts_TEtrasncripts.gtf
sed -i 's/repName/family_id/' transcripts_TEtrasncripts.gtf 
sed -i 's/intModel/family_id/' transcripts_TEtrasncripts.gtf
sed -i 's/category/class_id/' transcripts_TEtrasncripts.gtf
sed -i 's/repClass/class_id/' transcripts_TEtrasncripts.gtf
```

RSEM Database:

``` bash
conda activate RSEMSLM
mkdir rsem_DB
rsem-prepare-reference -gtf transcripts.gtf --bowtie2 hg38.fa rsem_DB/reference
```

SalmonTE Database:

``` bash
mkdir salmonTE_DB
awk '{if($0~/^>/){print $1,"HERV","Homo sapiens"}else print $0;}' OFS="\t" rsem_DB/reference.transcripts.fa > salmonTE_DB/reference_transcripts.fa
echo "HERV,HERV,Endogenous Retrovirus,Transposable Element" >> SalmonTE/reference/clades_extended.csv
SalmonTE/SalmonTE.py index --input_fasta=salmonTE_DB/reference_transcripts.fa --ref_name=salmonTE_DB
conda deactivate
```

Telescope does not need a pre-build database and works with the gtf provided.
