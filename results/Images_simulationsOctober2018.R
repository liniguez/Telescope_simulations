
data<-read.delim("allSim_numbers_2.txt")

tiff("Best_counts.tiff", res=300, width=960, height=960, pointsize=7)
boxplot(data$best_counts~data$real,  main="Best Random", xlab = "Simulated Expression", ylab= "Method Count")
stripchart(best_counts~real, vertical = TRUE, data =data, jitter=0.3,method = "jitter", add = TRUE, pch = 20,cex=0.75, col = rgb(0,0,1,alpha=0.50))
abline(a=-30,b=30, col=rgb(1,0,0,alpha=0.75), lty=4, lwd=3)
dev.off()

tiff("Unique_counts.tiff", res=300, width=960, height=960, pointsize=7)
boxplot(data$uni_counts~data$real,  main="Unique", xlab = "Simulated Expression", ylab= "Method Count",ylim=c(0,300))
stripchart(uni_counts~real, vertical = TRUE, data =data, jitter=0.3,method = "jitter", add = TRUE, pch = 20,cex=0.75, col = rgb(0,1,0,alpha=0.50))
abline(a=-30,b=30, col=rgb(1,0,0,alpha=0.75), lty=4, lwd=3)
dev.off()

tiff("RepEnrich_counts.tiff", res=300, width=960, height=960, pointsize=7)
boxplot(data$repenrich_counts~data$real,  main="RepEnrich", xlab = "Simulated Expression", ylab= "Method Count")
stripchart(repenrich_counts~real, vertical = TRUE, data =data, jitter=0.3,method = "jitter", add = TRUE, pch = 20,cex=0.75, col = rgb(1,0.4,0,alpha=0.50))
abline(a=-30,b=30, col=rgb(1,0,0,alpha=0.75), lty=4, lwd=3)
dev.off()

tiff("Telescope_counts.tiff", res=300, width=960, height=960, pointsize=7)
boxplot(data$telescope~data$real,  main="Telescope", xlab = "Simulated Expression", ylab= "Method Count")
stripchart(telescope~real, vertical = TRUE, data =data, jitter=0.3,method = "jitter", add = TRUE, pch = 20,cex=0.75, col = rgb(1,0,1,alpha=0.50))
abline(a=-30,b=30, col=rgb(1,0,0,alpha=0.75), lty=4, lwd=3)
dev.off()

tiff("TEtranscripts_counts.tiff", res=300, width=960, height=960, pointsize=7)
boxplot(data$tetr_counts~data$real,  main="TEtranscripts", xlab = "Simulated Expression", ylab= "Method Count")
stripchart(tetr_counts~real, vertical = TRUE, data =data, jitter=0.3,method = "jitter", add = TRUE, pch = 20,cex=0.75, col = rgb(0,1,1,alpha=0.50))
abline(a=-30,b=30, col=rgb(1,0,0,alpha=0.75), lty=4, lwd=3)
dev.off()


tiff("SalmonTE_counts.tiff", res=300, width=960, height=960, pointsize=7)
boxplot(data$SalmonTE~data$real,  main="TEtranscripts", xlab = "Simulated Expression", ylab= "Method Count")
stripchart(SalmonTE~real, vertical = TRUE, data =data, jitter=0.3,method = "jitter", add = TRUE, pch = 20,cex=0.75, col = rgb(1,0.07,0.5,alpha=0.50))
abline(a=-30,b=30, col=rgb(1,0,0,alpha=0.75), lty=4, lwd=3)
dev.off()

tiff("RSEM_counts.tiff", res=300, width=960, height=960, pointsize=7)
boxplot(data$RSEM~data$real,  main="TEtranscripts", xlab = "Simulated Expression", ylab= "Method Count")
stripchart(RSEM~real, vertical = TRUE, data =data, jitter=0.3,method = "jitter", add = TRUE, pch = 20,cex=0.75, col = rgb(0.5,0.16,0.88,alpha=0.50))
abline(a=-30,b=30, col=rgb(1,0,0,alpha=0.75), lty=4, lwd=3)
dev.off()



vp<-matrix(ncol=7, nrow=length(data$best_counts))
fn<-matrix(ncol=7, nrow=length(data$best_counts))
fp<-matrix(ncol=7, nrow=length(data$best_counts))

logi<-(data$best_counts - data$real)<=0
vp[logi,1]<-data$best_counts[logi]
vp[!logi,1]<-data$real[!logi]
fn[logi,1]<-abs(data$best_counts[logi] - data$real[logi])
fn[!logi,1]<-0
fp[logi,1]<-0
fp[!logi,1]<-data$best_counts[!logi] - data$real[!logi]

logi<-(data$uni_counts - data$real)<=0
vp[logi,2]<-data$uni_counts[logi]
vp[!logi,2]<-data$real[!logi]
fn[logi,2]<-abs(data$uni_counts[logi] - data$real[logi])
fn[!logi,2]<-0
fp[logi,2]<-0
fp[!logi,2]<-data$uni_counts[!logi] - data$real[!logi]

logi<-(data$repenrich_counts - data$real)<=0
vp[logi,3]<-data$repenrich_counts[logi]
vp[!logi,3]<-data$real[!logi]
fn[logi,3]<-abs(data$repenrich_counts[logi] - data$real[logi])
fn[!logi,3]<-0
fp[logi,3]<-0
fp[!logi,3]<-data$repenrich_counts[!logi] - data$real[!logi]

logi<-(data$telescope - data$real)<=0
vp[logi,4]<-data$telescope[logi]
vp[!logi,4]<-data$real[!logi]
fn[logi,4]<-abs(data$telescope[logi] - data$real[logi])
fn[!logi,4]<-0
fp[logi,4]<-0
fp[!logi,4]<-data$telescope[!logi] - data$real[!logi]

logi<-(data$tetr_counts - data$real)<=0
vp[logi,5]<-data$tetr_counts[logi]
vp[!logi,5]<-data$real[!logi]
fn[logi,5]<-abs(data$tetr_counts[logi] - data$real[logi])
fn[!logi,5]<-0
fp[logi,5]<-0
fp[!logi,5]<-data$tetr_counts[!logi] - data$real[!logi]

logi<-(data$RSEM - data$real)<=0
vp[logi,6]<-data$RSEM[logi]
vp[!logi,6]<-data$real[!logi]
fn[logi,6]<-abs(data$RSEM[logi] - data$real[logi])
fn[!logi,6]<-0
fp[logi,6]<-0
fp[!logi,6]<-data$RSEM[!logi] - data$real[!logi]

logi<-(data$SalmonTE - data$real)<=0
vp[logi,7]<-data$SalmonTE[logi]
vp[!logi,7]<-data$real[!logi]
fn[logi,7]<-abs(data$SalmonTE[logi] - data$real[logi])
fn[!logi,7]<-0
fp[logi,7]<-0
fp[!logi,7]<-data$SalmonTE[!logi] - data$real[!logi]


tempseq<-vector()
presicion_samp<-matrix(ncol = 7,nrow = 25)
recall_samp<-matrix(ncol = 7,nrow = 25)
F1_score <-matrix(ncol = 7,nrow = 25)
for (i in 1:25){
 tempseq<- seq(from=(i-1)*11+1,to = i*11,by=1)
 presicion_samp[i,]<-colSums(vp[tempseq,])/(colSums(vp[tempseq,])+colSums(fp[tempseq,]))
 recall_samp[i,]<-colSums(vp[tempseq,])/(colSums(vp[tempseq,])+colSums(fn[tempseq,]))
 F1_score[i,] <- 2* presicion_samp[i,]*recall_samp[i,]/(presicion_samp[i,]+recall_samp[i,])
} 

presicion<-colSums(vp)/(colSums(vp)+colSums(fp))
recall<-colSums(vp)/(colSums(vp)+colSums(fn))

tiff("Presicion_Recall.tiff", res=300, width=960, height=960, pointsize=7)
plot(recall_samp, presicion_samp,xlim=c(0,1),ylim=c(0,1), ylab="Presicion",xlab="Recall", pch=18,cex=1,
     col=c(rep(rgb(0,0,1,alpha=0.50),25),
           rep(rgb(0,1,0,alpha=0.50),25),
           rep(rgb(1,0.4,0,alpha=0.50),25),
           rep(rgb(1,0,1,alpha=0.50),25),
           rep(rgb(0,1,1,alpha=0.50),25),
           rep(rgb(0.5,0.16,0.88,alpha=0.50),25),
           rep(rgb(1,0.07,0.5,alpha=0.50),25)
           ))
points( recall,presicion, col=c(rgb(0,0,1), rgb(0,1,0),rgb(1,0.4,0),rgb(1,0,1),rgb(0,1,1),rgb(0.5,0.16,0.88),rgb(1,0.07,0.5)), pch=8,cex=3)
legend("bottomleft", legend = c("Best", "Unique", "RepEnrich", "Telescope", "TEtranscripts","RSEM","SalmonTE"),
  col = c(rgb(0,0,1), rgb(0,1,0),rgb(1,0.4,0),rgb(1,0,1),rgb(0,1,1),rgb(0.5,0.16,0.88),rgb(1,0.07,0.5)), pch=18,bty = "n",pt.cex=1.5)
dev.off()



F1_score2<-cbind(as.vector(F1_score),c(rep(rgb(0,0,1),25),rep(rgb(0,1,0),25),rep(rgb(1,0.4,0),25),rep(rgb(1,0,1),25),rep(rgb(0,1,1),25),rep(rgb(0.5,0.16,0.88),25),rep(rgb(1,0.07,0.5),25)))
F1_score3<-F1_score2[order(F1_score2[,1]),]

tiff("F1_score.tiff", res=300, width=960, height=960, pointsize=7)
barplot(as.numeric(F1_score3[,1]), col=F1_score3[,2],border = NA, ylab="F1 score")
legend("bottom",legend = c("Best", "Unique", "RepEnrich", "Telescope", "TEtranscripts","RSEM","SalmonTE"),
       col = c(rgb(0,0,1), rgb(0,1,0),rgb(1,0.4,0),rgb(1,0,1),rgb(0,1,1),rgb(0.5,0.16,0.88),rgb(1,0.07,0.5)),pch=15,pt.cex=1.5,bty = "n", xpd = TRUE,ncol=2, inset = c(0,-0.25),cex=1)
dev.off()



