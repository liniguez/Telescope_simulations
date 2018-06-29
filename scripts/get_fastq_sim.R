##########################
#
#files needed:
#	params_Sim.txt (file with parameters for Simulations (Matthew's script)
#	HML2.gtf (gtf of all HML2 of Matthew's annotation hg38)
#	/home/ubuntu/efs/simulations/telescope_june_2018/fastas (folder with a single fasta file for each chromosome hg38)
#
#######################

library(polyester)
library("Biostrings")

param_file<-"params_Sim.txt"
fastas_folder<- "/home/ubuntu/efs/simulations/telescope_june_2018/fastas"
gtf_file<- "HML2.gtf"
# Load the parameters
p <- read.table(param_file, sep='\t', stringsAsFactors=F)
params <- lapply(1:nrow(p), function(r) {
  v <- p[r, 2]
  if(grepl(',', v)) {
    vl <- as.integer(unlist(strsplit(v, ',')))
    return(vl)
  }
  suppressWarnings(
    if (is.na(as.numeric(v))) {
      if (is.na(as.logical(v))) {
        ret <- v 
      } else {
        ret <- as.logical(v)
      }
    } else {
      ret <- as.numeric(v)
    }
  )
  ret
})
names(params) <- p$V1


sequen<-seq_gtf(gtf=gtf_file, seqs=fastas_folder)
for (i in 1:25){
  set.seed(i+params$seed)
  transc<-sample(names(sequen))[1:10]
  set.seed(i+params$seed)
  expr<-sample(seq(30,300,by=30))
  countmat<-data.frame(row.names=transc,expression=expr)
  writeXStringSet(x=sequen[rownames(countmat)], filepath="temp_seq4sim.fa")
  outdir_temp<-paste0(params$outdir,"_",i)
  file_count<-paste0(outdir_temp,"/countmat.txt")
  simulate_experiment_countmat(readmat = as.matrix(countmat),fasta="temp_seq4sim.fa", outdir=outdir_temp, fraglen = params$fraglen, fragsd = params$fragsd, readlen = params$readlen, error_model = params$error_model,error_rate = params$error_rate,seed = (i+params$seed))
  write.table(countmat,file=file_count,quote=F,sep="\t")
}

