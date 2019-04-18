Counting HML2
================
Luis P Iniguez
4/15/2019

After running all TE counting algortihms the HML2 count was pulled as following:

Best and unique count were calculated with htseq-count. For best count, default parameters from bowtie2 were used. In contrast, for unique reads an additional step of cleaning multiple mapped reads out of the bam from -k 100 was needed and then the cleaned sam file was used as the input for htseq-count.

``` bash
for i in $(ls -d HML2_*)
do
SAMP=$i

#Best#
htseq-count --format=bam --order=pos --stranded=no --idattr=locus --minaqual=0 ${SAMP}/default_sorted.bam transcripts.gtf > ${SAMP}/best_count.txt

#Unique#
samtools view -@ $(nproc) -f 0x100 ${SAMP}/k100_sorted.bam | cut -f1 | sort -u > ${SAMP}/temp_reads.txt #get not primary reads, which means multiple mapped reads. 
perl -e '{ open (IN,"$ARGV[0]");
  while ($l=<IN>){chomp $l; $h{$l}=1;}close(IN);
  open(SAM,"$ARGV[1]"); while(<SAM>){
    @vec=split("\t",$_);
    if(!$h{$vec[0]}){ # if the read is not one of the multiple mapped reads then print it a new sam file
    print $_;}}}' ${SAMP}/temp_reads.txt ${SAMP}/k100.sam | htseq-count --format=sam --order=name --stranded=no --idattr=locus --minaqual=0 - transcripts.gtf > ${SAMP}/unique_count.txt
rm ${SAMP}/temp_reads.txt 

done
```

The rest of the softwares had their unique output and in order to put everything together this was performed:

``` r
for (i in grep(list.dirs(recursive = F), pattern="HML2_", value=T)){

#Simulated count matrix aka real
tempf<-file.path(i,"countmat.txt") 
countmat <- read.table(tempf)
countmat["__no_feature",] <- sum(countmat[grep("HML2", rownames(countmat),invert=T),]) #Summing all non-TE regions read counts, __no_feature means either not mapped or mapped outside the annotated regions

countmat$locus<-rownames(countmat)
countmat<-countmat[grep("HML2|__no_feature", rownames(countmat),),] #remove non-TE regions counts
colnames(countmat)<-c("real","locus")

#Best Count
tempf<-file.path(i,"best_count.txt")
bestc <- read.table(tempf, col.names=c('locus','best_counts'),sep='\t', stringsAsFactors=F)
bestc <- bestc[bestc[,2]>0,] #all features with 0 count are removed
bestc[bestc$locus == '__no_feature',2] <- bestc[bestc$locus == '__no_feature',2]+sum(countmat$real)-sum(bestc[,2]) #calculate __no_feature counts

#Unique Count
tempf<-file.path(i,"unique_count.txt")
unic <-read.table(tempf, col.names=c('locus','uni_counts'),sep='\t', stringsAsFactors=F)
unic<- unic[unic[,2]>0,]
unic[unic$locus == '__no_feature',2] <- unic[unic$locus == '__no_feature',2]+sum(countmat$real)-sum(unic[,2])

#RepEnrich
tempf<-file.path(i,"RepEnrich_results_fraction_counts.txt")
rec1 <-read.table(tempf, col.names=c('locus', 'class', 'family', 'repenrich_counts'), sep='\t', stringsAsFactors=F)
rec1 <- rec1[,c('locus','repenrich_counts')]
rec1[,2]<- round(rec1[,2]/2,0) #RepEnrich conts paired-end reads double as single end reads, therfore  division is needed
rec1 <- rec1[rec1[,2]>0,]
rec1<-rbind(rec1,c("__no_feature",(sum(countmat$real)-sum(rec1[,2]))))

#TEtrasncript
tempf<-file.path(i,"TEtranscript_count.txt")
tetra<-read.table(tempf, col.names=c('names','tetr_counts','tetr_counts2'), stringsAsFactors=F)[,c(1,2)]
tetra$locus<-as.factor(sapply(strsplit(tetra$names, split=":"), `[[`, 1))
tetra[,2]<-as.numeric(gsub(" ", "", tetra[,2], fixed = TRUE))
tetra<-aggregate(tetr_counts ~ locus, tetra, sum)
tetra$locus<-as.character(tetra$locus)
tetra<- tetra[tetra[,2]>0,]
tetra<-rbind(tetra,c("__no_feature",(sum(countmat$real)-sum(tetra[,2]))))

#RSEM
tempf<-file.path(i,"rsem.genes.results")
rstable<-read.table(tempf,sep='\t', header=T, stringsAsFactors=F)
rstable<-rstable[,c(1,5)]
colnames(rstable) <- c("locus","RSEM")
rstable$RSEM <- round(rstable$RSEM,0)
rstable<- rstable[rstable[,2]>0,]
rstable<-rbind(rstable,c("__no_feature",(sum(countmat$real)-sum(rstable[,2]))))

#SalmonTE
tempf<-file.path(i,"salmonTE_all/EXPR.csv")
stable<-read.table(tempf,sep=',', header=T, stringsAsFactors=F)
colnames(stable) <- c("locus","SalmonTE")
stable$SalmonTE<- round(stable$SalmonTE,0)
stable<- stable[stable[,2]>0,]
stable<-rbind(stable,c("__no_feature",(sum(countmat$real)-sum(stable[,2]))))

#Telescope
tempf<-file.path(i,"telescope-telescope_report.tsv")
ttable<-read.table(tempf,sep='\t', header=T, stringsAsFactors=F, skip=1)
ttable<-ttable[,c(1,3)]
colnames(ttable) <- c("locus","telescope")
ttable<- ttable[ttable[,2]>0,]
ttable[ttable$locus == '__no_feature',2] <- ttable[ttable$locus == '__no_feature',2]+sum(countmat$real)-sum(ttable[,2])

merged <- bestc
merged <- merge(merged, unic, by='locus', all=T)
merged <- merge(merged, rec1, by='locus', all=T)
merged <- merge(merged, ttable, by='locus', all=T)
merged <- merge(merged, tetra, by='locus', all=T)
merged <- merge(merged, rstable, by='locus', all=T)
merged <- merge(merged, stable, by='locus', all=T)
merged <- merge(merged,countmat,by='locus',all=T)
merged[is.na(merged)] <- 0

merged<-matrix(as.numeric(unlist(merged[,-1])), ncol=8, dimnames=list(merged$locus,colnames(merged[,-1])))

others<-colSums(merged[merged[,8]==0,]) # this are other HML2 that presented expression rather the ones simulated

merged <- merged[!(merged[,8]==0) & !(rownames(merged) %in% c('__ambiguous','__not_aligned')),]
merged <- rbind(merged,'Others'=as.vector(others))
merged <- merged[order(as.numeric(merged[,8]), decreasing=T),]

tempf<-file.path(i,"table.txt")
write.table(merged, tempf, sep="\t", quote=F, row.names=T)
}
```
