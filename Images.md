Images
================
Luis P Iniguez
4/15/2019

First all counts need to be together in a file:

``` bash
head -n 1 HML2_SIM_1/table.txt > allSim_numbers.txt
awk 'FNR>1' HML2_SIM_*/table.txt >> allSim_numbers.txt
```

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ggplot2", quietly = TRUE,warn.conflicts=F)) BiocManager::install("ggplot2", version = "3.8")
if (!requireNamespace("cowplot", quietly = TRUE,warn.conflicts=F)) BiocManager::install("cowplot", version = "3.8")
if (!requireNamespace("RColorBrewer", quietly = TRUE,warn.conflicts=F)) BiocManager::install("RColorBrewer", version = "3.8")
```

``` r
library(ggplot2)
library(cowplot)
library("RColorBrewer")


data<-read.table("allSim_numbers.txt", skip=1, col.names = c("locus","best_counts","uni_counts","repenrich_counts","telescope","tetr_counts","RSEM","SalmonTE","real"))
data<-cbind(data,cathegorie=data$real)
data$cathegorie[data$locus=="__no_feature"]<-"Not HERV/\nNot Mapped"
data$cathegorie[data$locus=="Others"]<-"Other HERV"
data$cathegorie<-factor(data$cathegorie,levels=c("Not HERV/\nNot Mapped","Other HERV",as.character(seq(from=30,to=300,by=30))))

colorp<-brewer.pal(n = 7, name = "Dark2")

tem<- theme( text = element_text(size=10),
              legend.text=element_text(size=7),
              panel.grid.major = element_line(colour = "gray80"),
              panel.border = element_blank(),
              axis.text.x = element_text(angle=45,hjust = 1,lineheight = 0.7),
              #strip.text.x = element_text(angle = 0),
              #axis.ticks = element_blank(),
              #axis.title.x=element_blank(),
              #axis.text.x=element_blank(),
              panel.grid.minor.x=element_blank(),
              panel.grid.major.x=element_blank(),
              strip.text.x = element_text(size=7),
              strip.text.y = element_text(size=7),
              plot.margin=margin(l=0,unit="cm", b=0,r=0),
              strip.background = element_rect(colour="white",fill="gray80"),
              legend.position="bottom")

plot_box_sim<-function(data2plot,x2plot,y2plot,title,color2plot,tem2plot=tem,limit){
  res<-ggplot(data=data2plot,mapping= aes_string(y=y2plot,x=x2plot))+
    geom_boxplot(outlier.shape = NA, fill=color2plot, alpha=0.5)+ggtitle(title)+
    geom_jitter(shape=16, position=position_jitter(0.25), colour= "black", alpha=0.75)+
    geom_segment(aes(x=2,y=0,xend=12.5,yend=315),colour=rgb(1,0,0),linetype="dashed",size=1.1)+
    geom_segment(aes(x=0.5,y=limit,xend=1.5,yend=limit),colour=rgb(1,0,0),linetype="dashed",size=1.1)+
    ylab("Read Count")+ xlab("")+theme_bw()+tem2plot
  return(res)
}
```

``` r
pbest<-plot_box_sim(data2plot=data,x2plot="cathegorie",y2plot="best_counts", title="Best Random",color2plot=colorp[1],limit=data$real[data$locus=="__no_feature"][1])
puniq<-plot_box_sim(data2plot=data,x2plot="cathegorie",y2plot="uni_counts", title="Unique",color2plot=colorp[2],limit=data$real[data$locus=="__no_feature"][1])
pRepEn<-plot_box_sim(data2plot=data,x2plot="cathegorie",y2plot="repenrich_counts", title="RepEnrich",color2plot=colorp[3],limit=data$real[data$locus=="__no_feature"][1])
pTele<-plot_box_sim(data2plot=data,x2plot="cathegorie",y2plot="telescope", title="Telescope",color2plot=colorp[4],limit=data$real[data$locus=="__no_feature"][1])
pTEt<-plot_box_sim(data2plot=data,x2plot="cathegorie",y2plot="tetr_counts", title="TEtranscripts",color2plot=colorp[5],limit=data$real[data$locus=="__no_feature"][1])
pREM<-plot_box_sim(data2plot=data,x2plot="cathegorie",y2plot="RSEM", title="RSEM",color2plot=colorp[6],limit=data$real[data$locus=="__no_feature"][1])
pSalm<-plot_box_sim(data2plot=data,x2plot="cathegorie",y2plot="SalmonTE", title="SalmonTE",color2plot=colorp[7],limit=data$real[data$locus=="__no_feature"][1])

plot_grid(puniq,pbest,pRepEn,pTEt,pREM,pSalm,pTele,  nrow = 2)
```

![](Images_files/figure-markdown_github/fig1-1.png)

All reads can be divided in two categories, depending from where they come from, HML2 reads or non-TE reads. Those reads can then be correctly or incorrectly mapped,depending of the outcome of the counting method, leading to 4 different categories: a) reads map to HML2 correctly (True Positive) b) reads map to HML2 incorrectly (False Positive) c) reads not mapped correctly (True Negative) d) reads not mapped incorrectly (False Negative)

It is worth to mention that all the assumptions were made based on the counted reads and not on the proper mapping, since the majority of the algorithms tested provide only the count and not the mapping file.

``` r
data2est<-data[-grep("__no_feature",data$locus),] #all information to mapped to any HERV
dataNOMA<-data[grep("__no_feature",data$locus),colnames(data)[2:9]] #all information of non-mapped reads

mcor<-matrix(ncol=8, nrow=length(data2est$best_counts)) #map correctly
mincor<-matrix(ncol=8, nrow=length(data2est$best_counts)) #map incorrectly
logunma<-matrix(ncol = 8,nrow = 25) #logical more unmapped than simulated

x<-colnames(data)[-1]

for (i in 1:8){
 logi<-(data2est[,x[i]] <= data2est$real)   #has more or less than real read count
 mcor[logi,i]<-data2est[logi,x[i]]            #those with less, all reads are mapped correctly to the HML2
 mcor[!logi,i]<-data2est$real[!logi]          #those with more, only the real reads are mapped correctly
 mincor[!logi,i]<-data2est[!logi,x[i]]-data2est$real[!logi] #those with more, who many were incorrectly mapped
 mincor[logi,i]<-0                            #those with less, non reads are mapped incorrectly
 logunma[,i]<-as.logical(dataNOMA[,x[i]]>dataNOMA$real) #has more unmapped reads than real?
}

tempseq<-vector()
tempVN<-numeric()

mcorHML2<-matrix(ncol = 8,nrow = 25)    #reads map to HML2 correcty
mapincor<-matrix(ncol = 8,nrow = 25)    #reads map to HML2 incorrecty
unmapnoHML2<-matrix(ncol = 8,nrow = 25) #reads not mapped correctly
unmapHML2<-matrix(ncol = 8,nrow = 25)   #reads not mapped incorrectly 

for (i in 1:25){                              #calculation for each simulation
  tempseq<- seq(from=(i-1)*11+1,to = i*11,by=1)
  mcorHML2[i,]<-colSums(mcor[tempseq,1:8])    #sum of map correctly
  mapincor[i,]<-colSums(mincor[tempseq,1:8])  #sum of map incorrectly
  
  unmapnoHML2[i,logunma[i,]] <- dataNOMA$real[i] #if there are more than real unmapped reads then all real reads from non-TE are correctly assinged = not mapped
  unmapnoHML2[i,!logunma[i,]] <- as.numeric(dataNOMA[i,!logunma[i,]])  #if there are less than real unmapped reads only the data from the counting method is counted 
  unmapHML2[i,logunma[i,]] <- as.numeric(dataNOMA[i,logunma[i,]])-dataNOMA$real[i] #how many are incorrectly mapped outside the annotations
  unmapHML2[i,!logunma[i,]] <- 0 #if unmapped reads are less than real unmapped reads there are no HML2 reads incorrectly assinged to a non-TE or unmapped
} 
```

![](Images_files/figure-markdown_github/fig2a-1.png)

``` r
presicion_samp<-mcorHML2[,1:7]/(mcorHML2[,1:7]+mapincor[,1:7])
recall_samp<-mcorHML2[,1:7]/(mcorHML2[,1:7]+unmapHML2[,1:7])
F1_score<- stack(as.data.frame(2*presicion_samp*recall_samp/(presicion_samp+recall_samp)))
presicion_samp<- stack(as.data.frame(presicion_samp))
recall_samp<- stack(as.data.frame(recall_samp))
colnames(presicion_samp)<-c("Presicion","Counting")
colnames(recall_samp)<-c("Recall","Counting")
colnames(F1_score)<-c("F1","Counting")

presicion<-stack(as.data.frame(t(colSums(mcorHML2[,1:7])/(colSums(mcorHML2[,1:7])+colSums(mapincor[,1:7])))))
recall<-stack(as.data.frame(t(colSums(mcorHML2[,1:7])/(colSums(mcorHML2[,1:7])+colSums(unmapHML2[,1:7])))))

colnames(presicion)<-c("Presicion","Counting")
colnames(recall)<-c("Recall","Counting")

forplt<-data.frame(Presicion=presicion_samp$Presicion,
                   Recall=recall_samp$Recall,
                   Counting=factor(presicion_samp$Counting))
forplt2<-data.frame(Presicion=presicion$Presicion,
                    Recall=recall$Recall,
                    Counting=factor(presicion$Counting))

pres<-ggplot(data=forplt,aes(x=Recall,y=Presicion, colour=Counting))+
  geom_point(alpha=0.75)+
  scale_color_manual(values=colorp)+
  geom_point(data=forplt2,aes(x=Recall,y=Presicion),shape=8, size=8)+
  guides(col=guide_legend(nrow=2, byrow = T))+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.margin=margin(l=0,unit="cm", b=0,r=0),
        strip.text.x = element_blank(),
        legend.spacing.x = unit(0.005,"native"),
        legend.margin = margin(t=0,r=0,b=0,l=0))
pres
```

![](Images_files/figure-markdown_github/fig3-1.png)

``` r
levels(F1_score$Counting)<-c("Best", "Unique", "RepEnrich", "Telescope", "TEtranscripts", "RSEM","SalmonTE")
f1<-ggplot(F1_score, aes(x=F1, color=Counting, fill=Counting))+
  geom_density()+
  scale_color_manual(values=colorp)+
  scale_fill_manual(values=alpha(colorp,0.5))+
  xlab("F1 score")+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        #legend=guide_legend(nrow = 2),
        strip.text.x = element_blank(),
        legend.spacing.x = unit(0.005,"npc"),
        plot.margin=margin(l=0,unit="cm", b=0,r=0),
        axis.text.x = element_text(angle=45,hjust = 1))+
  guides(col=guide_legend(nrow=2, byrow = T, keywidth = 1), size=2)
f1
```

![](Images_files/figure-markdown_github/fig4-1.png)

The paper plot:

``` r
plot_grid(puniq,pbest,pRepEn,pTEt,pREM,pSalm,pTele,pres,f1, labels="AUTO", nrow = 3)
```

![](Images_files/figure-markdown_github/paper-1.png)
