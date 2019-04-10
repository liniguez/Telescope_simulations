library(ggplot2)
library(cowplot)
library(viridis)


totreads=2100
notHERV=3
notHERVreads=450

data<-read.delim("allSim_numbers_7.txt")
data$real[data$real==15]<-0
data$real[data$locus == "NOFEAT"]<-notHERVreads
#x<-c(seq(from=300,to=30,by=-30),rep("delete",notHERV),"Other HERV","Not HERV")
x<-data$real
data<-cbind(data,cathegorie=x)
data<-data[!(data$locus=="nonHML2"),]
data$cathegorie[data$locus=="NOFEAT"]<-"Not HERV/\nNot Mapped"
data$cathegorie[data$locus=="Others"]<-"Other HERV"

#data<-data[!(data$cathegorie=="delete"),]

data$cathegorie<-factor(data$cathegorie,levels=c("Not HERV/\nNot Mapped","Other HERV",as.character(seq(from=30,to=300,by=30))))
for(i in seq(from=1,to=300,by=12)){
  data$SalmonTE[i+11]<-(totreads-sum(data$SalmonTE[i:(i+10)]))
  data$tetr_counts[i+11]<-(totreads-sum(data$tetr_counts[i:(i+10)]))
  
  data$uni_counts[i+11]<-(totreads-sum(data$uni_counts[i:(i+10)]))
  data$repenrich_counts[i+11]<-max(c(0,(totreads-sum(data$repenrich_counts[i:(i+10)]))))
  data$telescope[i+11]<-(totreads-sum(data$telescope[i:(i+10)]))
  
  data$RSEM[i+11]<-(totreads-sum(data$RSEM[i:(i+10)]))
}

colores<-plasma(7)
colores2<-plasma(7,alpha=0.5)
#colores<-c(rgb(0,0,1),rgb(0,1,0),rgb(1,0.4,0),rgb(1,0,1),rgb(0,1,1),rgb(1,0.07,0.5),rgb(0.5,0.16,0.88))
#colores2<-c(rgb(0,0,1,alpha=0.5),rgb(0,1,0,alpha=0.5),rgb(1,0.4,0,alpha=0.5),rgb(1,0,1,alpha=0.5),rgb(0,1,1,alpha=0.5),rgb(1,0.07,0.5,alpha=0.5),rgb(0.5,0.16,0.88,alpha=0.5))
tema<- theme( text = element_text(size=24),
              legend.text=element_text(size=10),
              panel.grid.major = element_line(colour = "gray80"),
              panel.border = element_blank(),
              axis.text.x = element_text(angle=45,hjust = 1,lineheight = 0.7),
              #strip.text.x = element_text(angle = 0),
              #axis.ticks = element_blank(),
              #axis.title.x=element_blank(),
              #axis.text.x=element_blank(),
              panel.grid.minor.x=element_blank(),
              panel.grid.major.x=element_blank(),
              strip.text.x = element_text(size=15),
              strip.text.y = element_text(size=15),
              strip.background = element_rect(colour="white",fill="gray80"),
              legend.position="bottom")



plot_box_sim<-function(data2plot,x2plot,y2plot,title,color2plot,tema2plot=tema){
  res<-ggplot(data=data2plot,mapping= aes_string(y=y2plot,x=x2plot))+ 
    geom_boxplot(outlier.shape = NA)+ggtitle(title)+
    geom_jitter(shape=16, position=position_jitter(0.25), colour= color2plot,alpha=0.75, fill=color2plot)+
    geom_segment(aes(x=2,y=0,xend=12.5,yend=315),colour=rgb(1,0,0),linetype="dashed",size=1.1)+
    geom_segment(aes(x=0.5,y=notHERVreads,xend=1.5,yend=notHERVreads),colour=rgb(1,0,0),linetype="dashed",size=1.1)+
    ylab("Read Count")+ xlab("")+theme_bw()+tema2plot
  return(res)
}
  

pbest<-plot_box_sim(data2plot=data,x2plot="cathegorie",y2plot="best_counts", title="Best Random",color2plot=colores[1])
ggsave("Best_Random.tiff",plot=pbest)
pbest<-plot_box_sim(data2plot=data,x2plot="cathegorie",y2plot="uni_counts", title="Unique",color2plot=colores[2])
ggsave("Unique_counts.tiff",plot=pbest)
pbest<-plot_box_sim(data2plot=data,x2plot="cathegorie",y2plot="repenrich_counts", title="RepEnrich",color2plot=colores[3])
ggsave("RepEnrich_counts.tiff",plot=pbest)
pbest<-plot_box_sim(data2plot=data,x2plot="cathegorie",y2plot="telescope", title="Telescope",color2plot=colores[4])
ggsave("Telescope_counts.tiff",plot=pbest)
pbest<-plot_box_sim(data2plot=data,x2plot="cathegorie",y2plot="tetr_counts", title="TEtranscripts",color2plot=colores[5])
ggsave("TEtranscripts_counts.tiff",plot=pbest)
pbest<-plot_box_sim(data2plot=data,x2plot="cathegorie",y2plot="RSEM", title="RSEM",color2plot=colores[6])
ggsave("RSEM_counts.tiff",plot=pbest)
pbest<-plot_box_sim(data2plot=data,x2plot="cathegorie",y2plot="SalmonTE", title="SalmonTE",color2plot=colores[7])
ggsave("SalmonTE_counts.tiff",plot=pbest)





data2est<-data[-grep("NOFEAT",data$locus),]
dataNOMA<-data[grep("NOFEAT",data$locus),c(2,5,6,7,8,9,10,11)]

mcor<-matrix(ncol=8, nrow=length(data2est$best_counts))
mincor<-matrix(ncol=8, nrow=length(data2est$best_counts))
logunma<-matrix(ncol = 7,nrow = 25)

logi<-(data2est$best_counts - data2est$real)<=0
mcor[logi,1]<-data2est$best_counts[logi]
mcor[!logi,1]<-data2est$real[!logi]
mincor[!logi,1]<-abs(data2est$best_counts[!logi] - data2est$real[!logi])
mincor[logi,1]<-0
logunma[,1]<-dataNOMA$best_counts>dataNOMA$real


logi<-(data2est$uni_counts - data2est$real)<=0
mcor[logi,2]<-data2est$uni_counts[logi]
mcor[!logi,2]<-data2est$real[!logi]
mincor[!logi,2]<-abs(data2est$uni_counts[!logi] - data2est$real[!logi])
mincor[logi,2]<-0
logunma[,2]<-dataNOMA$uni_counts>dataNOMA$real

logi<-(data2est$repenrich_counts - data2est$real)<=0
mcor[logi,3]<-data2est$repenrich_counts[logi]
mcor[!logi,3]<-data2est$real[!logi]
mincor[!logi,3]<-abs(data2est$repenrich_counts[!logi] - data2est$real[!logi])
mincor[logi,3]<-0
logunma[,3]<-dataNOMA$repenrich_counts>dataNOMA$real

logi<-(data2est$telescope - data2est$real)<=0
mcor[logi,4]<-data2est$telescope[logi]
mcor[!logi,4]<-data2est$real[!logi]
mincor[!logi,4]<-abs(data2est$telescope[!logi] - data2est$real[!logi])
mincor[logi,4]<-0
logunma[,4]<-dataNOMA$telescope>dataNOMA$real

logi<-(data2est$tetr_counts - data2est$real)<=0
mcor[logi,5]<-data2est$tetr_counts[logi]
mcor[!logi,5]<-data2est$real[!logi]
mincor[!logi,5]<-abs(data2est$tetr_counts[!logi] - data2est$real[!logi])
mincor[logi,5]<-0
logunma[,5]<-dataNOMA$tetr_counts>dataNOMA$real

logi<-(data2est$RSEM - data2est$real)<=0
mcor[logi,6]<-data2est$RSEM[logi]
mcor[!logi,6]<-data2est$real[!logi]
mincor[!logi,6]<-abs(data2est$RSEM[!logi] - data2est$real[!logi])
mincor[logi,6]<-0
logunma[,6]<-dataNOMA$RSEM>dataNOMA$real

logi<-(data2est$SalmonTE - data2est$real)<=0
mcor[logi,7]<-data2est$SalmonTE[logi]
mcor[!logi,7]<-data2est$real[!logi]
mincor[!logi,7]<-abs(data2est$SalmonTE[!logi] - data2est$real[!logi])
mincor[logi,7]<-0
logunma[,7]<-dataNOMA$SalmonTE>dataNOMA$real

logi<-(data2est$real - data2est$real)<=0
mcor[logi,8]<-data2est$real[logi]
mcor[!logi,8]<-data2est$real[!logi]
mincor[!logi,8]<-abs(data2est$real[!logi] - data2est$real[!logi])
mincor[logi,8]<-0

tempseq<-vector()
tempVN<-numeric()

mcorHML2<-matrix(ncol = 7,nrow = 25)
mapincor<-matrix(ncol = 7,nrow = 25)
unmapHML2<-matrix(ncol = 7,nrow = 25)
unmapnoHML2<-matrix(ncol = 7,nrow = 25)



for (i in 1:25){
  tempseq<- seq(from=(i-1)*11+1,to = i*11,by=1)
  mcorHML2[i,]<-colSums(mcor[tempseq,1:7])
  mapincor[i,]<-colSums(mincor[tempseq,1:7])
  
  unmapHML2[i,as.numeric(which(logunma[i,]))] <- as.numeric(dataNOMA[i,which(logunma[i,])]-dataNOMA$real[i])
  unmapnoHML2[i,as.numeric(which(logunma[i,]))] <- dataNOMA$real[i]
  unmapHML2[i,as.numeric(which(!logunma[i,]))] <- 0
  unmapnoHML2[i,as.numeric(which(!logunma[i,]))] <- as.numeric(dataNOMA[i,as.numeric(which(!logunma[i,]))])
  
} 

colnames(mcorHML2)<-c("Best Random","Unique","RepEnrich","Telescope","TEtranscripts","RSEM","SalmonTE")
mcorHML2_st<-stack(as.data.frame(mcorHML2))
mcorHML2_st$ind<-factor(mcorHML2_st$ind,levels = c("Best Random","Unique","RepEnrich","Telescope","TEtranscripts","RSEM","SalmonTE"))
ggplot(data= mcorHML2_st)+
  geom_violin(aes(x=ind,y=values,fill=ind),trim = T,scale = "width")+
  geom_hline(yintercept=1650,colour=rgb(1,0,0),linetype="dashed",size=1.1)+
  scale_fill_manual(values=colores)+
  ylab("Correct mapped reads")+
  #scale_y_continuous(trans='log10')+
  theme(legend.position="none",axis.text.x = element_text(angle=45,hjust = 1,lineheight = 0.7),axis.title.x=element_blank())
ggsave("Correctly_mapped_HML2read.tiff")


colnames(mapincor)<-c("Best Random","Unique","RepEnrich","Telescope","TEtranscripts","RSEM","SalmonTE")
mapincor_st<-stack(as.data.frame(mapincor))
mapincor_st$ind<-factor(mapincor_st$ind,levels = c("Best Random","Unique","RepEnrich","Telescope","TEtranscripts","RSEM","SalmonTE"))
ggplot(data= mapincor_st)+
  geom_violin(aes(x=ind,y=values,fill=ind),trim = T,scale = "width")+
  geom_hline(yintercept=0,colour=rgb(1,0,0),linetype="dashed",size=1.1)+
  scale_fill_manual(values=colores)+
  ylab("Incorrect mapped reads")+
  #scale_y_continuous(trans='log10')+
  theme(legend.position="none",axis.text.x = element_text(angle=45,hjust = 1,lineheight = 0.7),axis.title.x=element_blank())
ggsave("Incorrectly_mapped_HML2read.tiff")

colnames(unmapHML2)<-c("Best Random","Unique","RepEnrich","Telescope","TEtranscripts","RSEM","SalmonTE")
unmapHML2_st<-stack(as.data.frame(unmapHML2))
unmapHML2_st$ind<-factor(unmapHML2_st$ind,levels = c("Best Random","Unique","RepEnrich","Telescope","TEtranscripts","RSEM","SalmonTE"))
ggplot(data= unmapHML2_st)+
  geom_violin(aes(x=ind,y=values,fill=ind),trim = T,scale = "width")+
  geom_hline(yintercept=0,colour=rgb(1,0,0),linetype="dashed",size=1.1)+
  scale_fill_manual(values=colores)+
  ylab("Incorrectly not mapped reads")+
  #scale_y_continuous(trans='log10')+
  theme(legend.position="none",axis.text.x = element_text(angle=45,hjust = 1,lineheight = 0.7),axis.title.x=element_blank())
ggsave("Incorrectly_mapped_NON-HML2read.tiff")

colnames(unmapnoHML2)<-c("Best Random","Unique","RepEnrich","Telescope","TEtranscripts","RSEM","SalmonTE")
unmapnoHML2_st<-stack(as.data.frame(unmapnoHML2))
unmapnoHML2_st$ind<-factor(unmapnoHML2_st$ind,levels = c("Best Random","Unique","RepEnrich","Telescope","TEtranscripts","RSEM","SalmonTE"))
ggplot(data= unmapnoHML2_st)+
  geom_violin(aes(x=ind,y=values,fill=ind),trim = TRUE,scale = "width")+
  geom_hline(yintercept=450,colour=rgb(1,0,0),linetype="dashed",size=1.1)+
  scale_fill_manual(values=colores)+
  ylab("Correctly not mapped reads")+
  #scale_y_continuous(trans='log10')+
  theme(legend.position="none",axis.text.x = element_text(angle=45,hjust = 1,lineheight = 0.7),axis.title.x=element_blank())
ggsave("Correctly_mapped_NON-HML2read.tiff")

all<-rbind(cbind(mcorHML2_st,color="Correctly mapped reads"),
           cbind(mapincor_st,color="Incorrectly mapped reads"),
           cbind(unmapHML2_st,color="Incorrectly unmapped reads"),
           cbind(unmapnoHML2_st,color="Correctly unmapped reads"))
levels(all$color)<-c("Correctly mapped reads","Incorrectly mapped reads","Incorrectly unmapped reads","Correctly unmapped reads")
ggplot(data=all)+
  geom_violin(aes(x=ind,y=values,fill=color),trim = TRUE,scale = "width")+
  facet_grid(cols=vars(ind), scales = "free_x",space = "free")+
  ylab("Reads")+
  theme(legend.position="bottom",
        legend.title = element_blank(),
        #legend=guide_legend(nrow = 2),
        strip.text.x = element_blank(),
        axis.text.x = element_text(angle=45,hjust = 1),
        axis.title.x=element_blank(), panel.spacing.x = unit(1,"line"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
ggsave("All_correctly.tiff")


presicion_samp<-mcorHML2/(mcorHML2+mapincor)
recall_samp<-mcorHML2/(mcorHML2+unmapHML2)
F1_score<- 2*presicion_samp*recall_samp/(presicion_samp+recall_samp)
specificity_samp<-unmapnoHML2/(unmapnoHML2+unmapHML2)

presicion<-colSums(mcorHML2)/(colSums(mcorHML2)+colSums(mapincor))
recall<-colSums(mcorHML2)/(colSums(mcorHML2)+colSums(unmapHML2))
specificity<-colSums(unmapnoHML2)/(colSums(unmapnoHML2)+colSums(unmapHML2))
tiff("Presicion_Recall.tiff", res=300, width=960, height=960, pointsize=7)
plot(recall_samp, presicion_samp,xlim=c(0,1),ylim=c(0,1), ylab="Presicion",xlab="Recall", pch=18,cex=1,
     col=c(rep(colores[1],25),rep(colores[2],25),rep(colores[3],25),rep(colores[4],25),
           rep(colores[5],25),rep(colores[6],25),rep(colores[7],25)))
points( recall,presicion, col=colores, pch=8,cex=3)
legend("bottomleft", legend = c("Best", "Unique", "RepEnrich", "Telescope", "TEtranscripts", "RSEM","SalmonTE"),
       col = colores, pch=18,bty = "n",pt.cex=1.5)
dev.off()

F1_score_st<-stack(as.data.frame(F1_score))
levels(F1_score_st$ind)<-c("Best", "Unique", "RepEnrich", "Telescope", "TEtranscripts", "RSEM","SalmonTE")
ggplot(F1_score_st, aes(x=values, color=ind, fill=ind))+
  geom_density()+
  scale_color_manual(values=colores)+
  scale_fill_manual(values=colores2)+
  xlab("F1 score")+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        #legend=guide_legend(nrow = 2),
        strip.text.x = element_blank(),
        legend.spacing.x = unit(0.005,"npc"),
        axis.text.x = element_text(angle=45,hjust = 1))+
  guides(col=guide_legend(nrow=1, byrow = T, keywidth = 1))
ggsave("F1_score_density.tiff")
