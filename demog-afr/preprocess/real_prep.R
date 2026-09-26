### preparation: get the best 5000 regions for real data
'%ni%' <- Negate('%in%')
options("scipen"=100)

### take the 5000 best windows
basdir="/lisc/data/scratch/admixlab/mk_data/"
library(GenomicRanges)
qra<-IRanges(start=seq(1,250000001,10000),end=seq(50000,250000000+50000,10000))
alovr<-list();ovrval<-list()
for (chrom in (1:22)) {
  print(chrom)
  chv<-read.table(paste(basdir,"/san/",chrom,".ver.txt.gz",sep=""),sep="\t",header=F)
  rr<-IRanges(start=chv[,1],end=chv[,1])
  ovr<-countOverlaps(qra,rr)
  alovr[[chrom]]<-qra[which(ovr>45000),]
  ovrval[[chrom]]<-ovr
}

ovrval<-unlist(ovrval)

allreg<-list()
for (chrom in (1:22)) { allreg[[chrom]]<-GRanges(seqnames=paste("chr",chrom,sep=""),ranges=qra) }
allreg<- unlist(GRangesList(allreg))

ovrid<-order(ovrval,decreasing=T)
ovrst<-ovrval;rereg<-allreg;vals<-c()
bestreg<-list()
i=0
repeat {
  i=i+1
  cn<-rereg[ovrid[1],]; vals=c(vals,ovrst[ovrid[1]])
  torm<-c(1,which(as.vector(poverlaps(cn,rereg))==T))
  ovrst<-ovrst[-torm]; ovrid<-order(ovrst,decreasing=T);rereg<-rereg[-c(torm),]
  bestreg[[i]]<-cn
  if (i==5000) { break }
}


bereg<-reduce( unlist(GRangesList(bestreg)))
save(bereg,file=paste(basdir,"/san/bestregions.Robject",sep=""))
bebed<-as.data.frame(bereg)[,c(1:3)]
bebed<-bebed[order(bebed[,2]),]
bebed<-bebed[order(as.numeric(do.call(rbind,strsplit(as.character(bebed[,1]),split="r"))[,2])),]
bebed[,2]<-bebed[,2]-1

write.table(bebed,file=paste(basdir,"/san/bestregions.bed",sep=""),sep="\t",row.names=F,col.names=F,quote=F)

