# Telescope_simulations

In order to prove Telescope, we simulate the expression of ten randomly choosen HML2. The level of expression of those HERV's started at 30 counts and increased by 30 read counts up to 300, the expression level was also randomly asignated. Paired-end reads were generated with POLYESTER bioconductor package. (script for generating reads: get_fastq_sim.R)

After reads were simulated bowtie2 was used for mapping all reads (`bowtie2 --score-min L,0,1.6 -k 100 --very-sensitive-local` script for bowtie2 run_bowtie.sh). Telescope was run inmediatly after using annotation file as the whole hg38 HERV annotation (transcripts_MB_hg38.gtf)

The results of Telescope were compared to other count methods: Best count, Unique count, RepEnrich (https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-583), TEtranscripts (https://academic.oup.com/bioinformatics/article/31/22/3593/240793). For best count bowtie2 was runed again with default parameters, but very-sensitive-local and score-min as previously(`bowtie2 --very-sensitive-local --score-min L,0,1.6`). For Unique count and TEtranscript the same mapping file as used for Telescoped was used. 
RepEnrich runs differently, first it needs a pseudo genome for repeated elements which is explained below, and it also needs two rounds of mapping. The first one identifies unique mapping reads and the second round looks for multiple mapping reads (as described in the RepEnrich methodology). The script used for the comparison of transcript count is competition_sim.sh

Two pseudo genomes were used for RepEnrich, first the default hg38 (https://github.com/nskvir/RepEnrich), and second a "hacked" version. These last pseudogenome each HML2 locus was considered a repeated element. 

Results:

![Unique counts](results/Unique_counts.tiff?raw=true "Title")

![Best counts](results/Best_counts.tiff?raw=true "Title")

![TEtranscripts counts](results/TEtranscripts_counts.tiff?raw=true "Title")

![Telescope counts](results/Telescope_counts.tiff?raw=true "Title")

![RepEnric counts](results/RepEnrich_counts.tiff?raw=true "Title")

![Presicion and Recall](results/Presicion_Recall.tiff?raw=true "Title")

![F1 score density](results/F1_score_density.tiff?raw=true "Title")
