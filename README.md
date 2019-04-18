# Statistical Performance of TE Quantification Methods

All scripts needed to examine the sensitivity and biases of computational approaches for quantifying TE expression: 1) unique counts, 2) best counts, 3) [RepEnrich](https://github.com/nerettilab/RepEnrich2), 4) [TEtranscripts](http://hammelllab.labsites.cshl.edu/software/), 5) [RSEM](http://deweylab.github.io/RSEM/), 6) [SalmonTE](https://github.com/LiuzLab/SalmonTE), and 7) [Telescope](https://github.com/mlbendall/telescope)

* [Data Bases](DB_creation.md)
* [Simulations](Simulations.md)
* [Software](Software.md)
* [Counting](Count_Plots.md)
* [Plots](Images.md)

### Telescope simulations

We simulated 25 independent data sets, each consisting of 10 randomly chosen HML-2 proviruses. In addition, we included 3 non-HML-2 loci in each data set to represent possible shortcomings with the annotation; in empirical data we often observe fragments that map outside the annotation. The non-HML-2 loci were selected randomly from a list of genomic regions not annotated as HML-2, but share high similarity with members of the HML-2 subfamily. Each HML-2 locus was expressed at a different level, ranging from 30 to 300 fragments per locus, while non-HML-2 loci were expressed at 150 fragments each. Using this expression pattern, we simulated sequencing fragments with the Bioconductor package for RNA-seq simulation, [Polyester](https://bioconductor.org/packages/release/bioc/html/polyester.html) . All simulations used the parameters of read length: 75 bp; average fragment size: 250; fragment size standard deviation: 25; and an Illumina error model with an error rate of 5e-3.

After reads were simulated bowtie2 was used for mapping all reads as well as all other mapping algorithms needed for other counting methods. [Telescope HERV annotations](https://github.com/mlbendall/telescope_annotation_db/tree/master/builds) on hg38 were used as reference for all counting methods

The results of Telescope were compared to other count methods: Best count, Unique count, RepEnrich, TEtranscripts, RSEM and SalmonTE. 
