library(demografr)
library(slendr)
init_env(uv=T)
.libPaths( c("/lisc/data/scratch/admixlab/mk_data/rlib/", .libPaths()) )  
#library(tibble, quietly = TRUE)
library(e1071, quietly = TRUE)
library(GenomicRanges, quietly = TRUE)
library(reticulate, quietly = TRUE)
library(dplyr)
library(scales)
py_require("scipy")

'%ni%' <- Negate('%in%')
options("scipen"=100)

## load simulations (should be 25k in total)
allda<-list()
it<-1; mis=c()
for (iter in c(1:3125)) {
  tt<-try(load(file=paste("/lisc/data/scratch/admixlab/mk_data/san/simul/simul_",iter,".Robject",sep="")))
  if (inherits(tt,"try-error")) { mis=c(mis,iter);next }
  allda[[it]]<-data
  it=it+1
}


allda<-combine_data(allda)

## load real data
load(file="/lisc/data/scratch/admixlab/mk_data/san/real_data_stats.Robject")
allda$observed<-observed
sallda<-allda

# select stats that make sense
patrn<-c("diversity_.*_med","divergence_.*_med","f2","fst_.*_med","td_.*_med","segs_.*_med","afss_med_|afss_mad_|afss_cov_|afss_fixed_","private.*_mean","sstar_.*-mean")
selc<-list();slc<-list()
for (j in (1:length(observed)))
{ selc[[j]]<-colnames(allda$simulated)[grep(patrn[j],colnames(allda$simulated))];slc[[j]]<-grep(patrn[j],colnames(allda$simulated)) }

# subset the observed and simulated data accordingly
nobserved<-list()
alcol<-colnames(allda$simulated)
j=0
while(j < length(observed)) { 
  j=j+1
  lnt<-nrow(observed[[j]])
  vls<-alcol[1:lnt]
  observed[[j]][,1]<-vls
  alcol<-alcol[-c(1:lnt)]
  nobserved[[j]]<-observed[[j]][which(observed[[j]][,1]%in%selc[[j]]),]
  }

obsv<-nobserved[[1]][,2];osn<-nobserved[[1]][,1]
for (j in (2:length(nobserved))) { obsv<-c(obsv,nobserved[[j]][,2]);osn=c(osn,nobserved[[j]][,1])}
names(obsv)<-osn
suballda<-sallda
suballda$observed<-nobserved
names(suballda$observed)<-names(sallda$observed)
suballda$simulated<-allda$simulated[,unlist(selc)]

## normalize values to values between 0 and 1 (across simulation and observation)
normvals <- as.data.frame(rbind(allda$simulated[,unlist(selc)],obsv)) %>%
  mutate(across(where(is.numeric), ~ rescale(.x, to = c(0, 1), na.rm = TRUE)))

simnorm<-normvals[-nrow(normvals),]
obsv<-unlist(normvals[nrow(normvals),])

fosn<-list(data.frame(set=osn[1:4],val=obsv[1:4]),data.frame(set=osn[5:10],val=obsv[5:10]), data.frame(pop=osn[11:16],f2=obsv[11:16]),
           data.frame(set=osn[17:22],val=obsv[17:22]),data.frame(set=osn[23:26],val=obsv[23:26]),
           data.frame(set=osn[27:30],val=obsv[27:30]),data.frame(names=osn[31:48],value=obsv[31:48]),
           data.frame(set=osn[49:52],val=obsv[49:52]),data.frame(set=osn[53:56],val=obsv[53:56]))
names(fosn)<-names(sallda$observed)

allda$observed<-fosn
allda$simulated<-simnorm


## run ABC & save values
myabc <- run_abc(allda, engine = "abc", tol = 0.05, method = "neuralnet",numnet=100) ## subset of stats with normalization
myabc2 <- run_abc(suballda, engine = "abc", tol = 0.05, method = "neuralnet",numnet=100) ## subset of stats

ee<-extract_summary(myabc)
save(myabc,file="/lisc/data/scratch/admixlab/mk_data/san/abc_inference.Robject")
write.table(ee,file="/lisc/data/scratch/admixlab/mk_data/san/abc_values.tsv",sep="\t",col.names=T,row.names=T,quote=F)


### diagnostic plots

#pdf("~/demog-afr/plots/post.pdf",6,6)
#plot_posterior(myabc, param = "gf") + ggplot2::coord_cartesian()
#plot_posterior(myabc, param = "Ne") + ggplot2::coord_cartesian()
#plot_posterior(myabc, param = "T") + ggplot2::coord_cartesian()
#dev.off()


load(file="/lisc/data/scratch/admixlab/mk_data/san/abc_inference.Robject")
cv_abc <- cross_validate(myabc, nval = 100, tols = c(0.01, 0.025,0.05,0.075,0.1,0.25), method = "neuralnet")
save(cv_abc,file="/lisc/data/scratch/admixlab/mk_data/san/abc_cross_val.Robject")



## confusion matrices
cv<-list(); posts<-list()
load(file="/lisc/data/scratch/admixlab/mk_data/san/simul_red/red_selection.Robject")
cv[["reduced"]]<-cv_sel
posts[["reduced"]]<-post_sel
load(file="/lisc/data/scratch/admixlab/mk_data/san/simul_wss/wss_selection.Robject")
cv[["wss"]]<-cv_sel
posts[["wss"]]<-post_sel
load(file="/lisc/data/scratch/admixlab/mk_data/san/simul_nog/nog_selection.Robject")
cv[["noghost"]]<-cv_sel
posts[["noghost"]]<-post_sel

library(grid)


