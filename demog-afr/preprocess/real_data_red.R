library(demografr)
library(slendr)
init_env(uv=T)
.libPaths( c("/lisc/data/scratch/admixlab/mk_data/rlib/", .libPaths()) )  
#library(tibble, quietly = TRUE)
library(e1071, quietly = TRUE)
library(GenomicRanges, quietly = TRUE)
library(reticulate, quietly = TRUE)
py_require("scipy")

'%ni%' <- Negate('%in%')
options("scipen"=100)

########## here for REAL data! --- for Nama & MSL only
### metadata & functions slightly adjusted for processing the ts file from real data

## names of combinations
idx_pairs <- combn(seq_along(1:2), 2, simplify = FALSE)
idx_nam<-paste(do.call(rbind,idx_pairs)[,1],do.call(rbind,idx_pairs)[,2],sep="-")
uninam=c("Nama","MSL")
combnam<-c("MSL-Nama","Nama-MSL")
names<-read.table("~/demog-afr/preprocess/Subset_WSS.txt",sep="\t",header=T)[1:8,c(4)]
samples<-list()
lens<-c()
popv<-c()
for (j in (uninam)) { samples[[j]]<-names[grep(j,names)];lens[j]<-length(samples[[j]]);popv<-c(popv,rep(j,lens[j])) }
hapsamp<-rep(names,each=2)


## the statistics
# nucleotide diversity in each population in windows (median & mad)
compute_diversity <- function(ts) {
  samples <- ts_names(ts, split = "pop")
  slen<-ts$sequence_length
  winsize=50000
  tm<-ts_diversity(ts, sample_sets = samples,windows=seq(winsize,slen-winsize,winsize))
  tm$med<-unlist(lapply(tm$diversity,median))*10000
  tm$mad<-unlist(lapply(tm$diversity,mad))*10000
  tm<-data.frame(set=paste(rep(tm$set,2),c(rep("med",2),rep("mad",2)),sep="_"),val=c(tm$med,tm$mad))
  return(tm)
}

# pairwise divergence between populations (dxy) in windows  (median & mad)
compute_divergence <- function(ts) {
  samples <- ts_names(ts, split = "pop")
  slen<-ts$sequence_length
  winsize=50000
  tm<-ts_divergence(ts, sample_sets = samples,windows=seq(winsize,slen-winsize,winsize))
  tm$med<-unlist(lapply(tm$divergence,median))*10000
  tm$mad<-unlist(lapply(tm$divergence,mad))*10000
  tm<-data.frame(set=paste(rep(paste(tm$x,tm$y,sep="-"),2),c(rep("med",1),rep("mad",1)),sep="_"),val=c(tm$med,tm$mad))
  return(tm)
}

# pairwise f2-statistics genome-wide
compute_f2 <- function(ts) {
  calcf2<-function(pops,ts) { ts_f2(ts,A=pops[[1]], B=pops[[2]]) }
  samples <- ts_names(ts, split = "pop")
  combnam<-c("Nama-MSL")
  idx_pairs <- combn(seq_along(1:2), 2, simplify = FALSE)
  sampair <- lapply(idx_pairs, function(ix) samples[ix])
  tm<-do.call(rbind,lapply(sampair,calcf2,ts=ts))
  tm$f2<-tm$f2*10000
  op<-data.frame(pop=combnam,f2=tm$f2)
  return(op)
}

# pairwise FST in windows  (median & mad)
compute_fst <- function(ts) {
  samples <- ts_names(ts, split = "pop")
  slen<-ts$sequence_length
  winsize=50000
  tm<-ts_fst(ts, sample_sets = samples,windows=seq(winsize,slen-winsize,winsize))
  tm$med<-unlist(lapply(tm$Fst,median,na.rm=T))*100
  tm$mad<-unlist(lapply(tm$Fst,mad,na.rm=T))*100
  tm<-data.frame(set=paste(rep(paste(tm$x,tm$y,sep="-"),2),c(rep("med",1),rep("mad",1)),sep="_"),val=c(tm$med,tm$mad))
  return(tm)
}

# population-wise Tajima's D in windows (median & mad)
compute_td <- function(ts) {
  samples <- ts_names(ts, split = "pop")
  slen<-ts$sequence_length
  winsize=50000
  tm<-ts_tajima(ts, sample_sets = samples,windows=seq(winsize,slen-winsize,winsize))
  tm$med<-unlist(lapply(tm$D,median,na.rm=T))
  tm$mad<-unlist(lapply(tm$D,mad,na.rm=T))
  tm<-data.frame(set=paste(rep(tm$set,2),c(rep("med",2),rep("mad",2)),sep="_"),val=c(tm$med,tm$mad))
  return(tm)
}

# number of segregating sites per population in windows (median & mad)
compute_segs <- function(ts) {
  samples <- ts_names(ts, split = "pop")
  slen<-ts$sequence_length
  winsize=50000
  tm<-ts_segregating(ts, sample_sets = samples,windows=seq(winsize,slen-winsize,winsize))
  tm$med<-unlist(lapply(tm$segsites,median))/100
  tm$mad<-unlist(lapply(tm$segsites,sd))/100
  tm<-data.frame(set=paste(rep(tm$set,2),c(rep("med",2),rep("mad",2)),sep="_"),val=c(tm$med,tm$mad))
  return(tm)
}

# AFS-based statistics per population, in subsets of 3: median & mad of AFS, pairwise correlation between populations, alleles at frequency 1 per population, fixed alleles per population
compute_afss <- function(ts) {
  library(e1071, quietly = TRUE)
  pstat<-function(input) { c(median(input)*1000,mad(input)*1000, kurtosis(input)) }
  bstat<-function(pr,afs) { cov(afs[[pr[1]]],afs[[pr[2]]],method="spearman") }
  samples <- ts_names(ts, split = "pop")
  smples <- lapply(samples, head, 3)
  names(smples)<-names(samples)
  afsr<-list(); afsl<-list()
  for (i in (1:length(smples))) {afsr[[i]]<-ts_afs(ts,sample_sets=list(smples[[i]])); afsl[[i]]<-afsr[[i]]/sum(afsr[[i]]) }
  names(afsl)<-names(samples)
  perp<-do.call(rbind,lapply(afsl,pstat))  
  colnames(perp)<-c("med","mad","kurt")
  idx_pairs <- combn(seq_along(1:2), 2, simplify = FALSE)
  idx_nam<-paste(do.call(rbind,idx_pairs)[,1],do.call(rbind,idx_pairs)[,2],sep="-")
  pcov<-unlist(lapply(idx_pairs,bstat,afs=afsl))
  names(pcov)<-idx_nam
  afspe<-do.call(cbind,afsr)[c(2,length(afsr[[1]])),]/c(1000,100)
  rownames(afspe)<-c("singleton","fixed")
  colnames(afspe)<-names(samples)
  retva<-unlist(list(perp,pcov,afspe))
  retva<-data.frame(names=c(paste("med_",names(samples),sep=""),paste("mad_",names(samples),sep=""),paste("kurt_",2:3,sep=""),paste("cov_",4,sep=""),paste("singl_",names(samples),sep=""),paste("fixed_",names(samples),sep="")),value=retva)
  return(retva)
}

# mean and sd private SNPs 
compute_private <- function(ts) {
  library(GenomicRanges, quietly = TRUE)
  compute_priv <- function(ts,sets,winsize) {
    privcomp<-function(region,posvec) { 
      return(rowSums(as.data.frame(mcols(subsetByOverlaps(posvec,IRanges(start=region[1],end=region[2]),)))>0))
    }
    samples <- ts_names(ts, split = "pop")
    slen<-ts$sequence_length
    winsize=50000
    ref=unlist(samples[sets[[2]]])
    tgt=unlist(samples[sets[[1]]])
    gtm<-ts_genotypes(ts,quiet=T)
    posvec<-IRanges(start=as.numeric(unlist(gtm[,1])),end=as.numeric(unlist(gtm[,1])))
    winvec<- split(cbind(seq(0,slen-winsize,winsize),seq(winsize,slen,winsize)),seq(length(seq(winsize,slen,winsize))))
    gtm<-gtm[,-1]
    odd<-seq_len(ncol(gtm))%%2
    gtfin<-gtm[,which(odd==1)]+gtm[,which(odd==0)]
    colnames(gtfin)<-gsub("_chr1|_chr2","",colnames(gtfin))
    selec<-which(rowSums(gtfin[,ref])>0)
    gtfin<-gtfin[selec,tgt]
    posvec<-posvec[selec,]
    mcols(posvec)$gtfin<-gtfin
    allvals<-unlist(lapply(winvec,privcomp,posvec=posvec))
    tm<-data.frame(set=paste(rep(names(samples)[sets[[1]]],2),c("mean","sd"),sep="_"),val=c(mean(allvals),sd(allvals)))
    return(tm)
  }
  samples <- ts_names(ts, split = "pop")
  combinset<-list(list(1,2),list(2,1))
  priv<-do.call(rbind,lapply(combinset,compute_priv,ts=ts,winsize=winsize))
  return(priv)
}    

## sstar mean and sd across individuals across windows; for different combinations
compute_sstar <- function(ts) {
  reticulate::py_run_string(r"(

import numpy as np
import pandas as pd
from typing import List, Tuple
from scipy.spatial import distance_matrix

def calc_sstar_summary(ts, window_size=50_000,
                      match_bonus: int = 5000, max_mismatch: int = 5,
                      mismatch_penalty: int = -10000) -> pd.DataFrame:
    """
    Calculate mean and standard deviation of sstar statistics for all specified target-reference combinations.

    Args:
        ts: tskit TreeSequence object
        window_size: Size of windows to calculate sstar in
        match_bonus: Bonus for matching genotypes
        max_mismatch: Maximum number of mismatches allowed
        mismatch_penalty: Penalty for mismatches

    Returns:
        DataFrame with columns: set, val
        where set is in format "KS-RHGn-mean" or "KS-RHGn-sd"
        and val is the corresponding value
    """
    # Get sample information
    samples = ts.nsamples

    # Define all target-reference combinations to analyze
    combinations = [
        ("Nama", "MSL"),
        ("MSL", "Nama")
    ]

    # Define the sstar calculation functions
    def _create_matrices(gt: np.ndarray, positions: np.ndarray, idx: int):
        hap = gt[:, idx]
        mask = hap != 0
        pos_sub = positions[mask]
        geno_matrix = gt[mask]
        keep = geno_matrix.sum(axis=1) != 1
        geno_matrix = geno_matrix[keep]
        pos_sub = pos_sub[keep]
        gd_mat = distance_matrix(geno_matrix, geno_matrix, p=1)
        p = pos_sub.astype(np.float64)
        pos_mat = np.tile(p, (p.size, 1))
        pd_mat = np.abs(pos_mat.T - pos_mat)
        return pd_mat, gd_mat, pos_sub

    def _calc_ind_sstar(pd_matrix, gd_matrix, positions_sub, match_bonus,
                       max_mismatch, mismatch_penalty):
        invalid = pd_matrix < 10
        scores = pd_matrix.copy()
        scores[invalid] = -np.inf
        scores[(~invalid) & (gd_matrix == 0)] += match_bonus
        scores[(~invalid) & (gd_matrix > 0) & (gd_matrix <= max_mismatch)] = mismatch_penalty
        scores[(~invalid) & (gd_matrix > max_mismatch)] = -np.inf

        n = positions_sub.size
        if n == 0:
            return 0.0

        max_scores = np.full(n, -np.inf, dtype=np.float64)
        for j in range(n):
            best = -np.inf
            for i in range(j):
                best = max(best, max_scores[i] + scores[j, i], scores[j, i])
            max_scores[j] = best

        if np.isneginf(max_scores).all():
            return 0.0
        return float(np.nanmax(max_scores))

    def sstar(tgt_gt: np.ndarray, pos: np.ndarray) -> List[float]:
        _, n_samples = tgt_gt.shape
        sstar_scores = []
        for i in range(n_samples):
            pdm, gdm, pos_sub = _create_matrices(tgt_gt, pos, i)
            score = _calc_ind_sstar(pdm, gdm, pos_sub, match_bonus,
                                  max_mismatch, mismatch_penalty)
            if score == -np.inf:
                score = 0.0
            sstar_scores.append(score)
        return sstar_scores

    def get_sstar_stats(tgt: str, ref: str) -> Tuple[float, float]:
        """Calculate mean and std of sstar for a specific tgt-ref combination"""
        # Identify target and reference samples
        tgt_samples = [i for i, s in enumerate(samples) if s.startswith(tgt + "_")]
        ref_samples = [i for i, s in enumerate(samples) if s.startswith(ref + "_")]

        # Combine samples
        selected_samples = tgt_samples + ref_samples
        all_sstar_scores = []

        # Process in windows
        sequence_length = int(ts.sequence_length)
        for left in range(0, sequence_length, window_size):
            right = min(left + window_size, sequence_length)

            # Get variants in window
            variants = list(ts.variants(samples=selected_samples, left=left, right=right))
            if not variants:
                continue

            # Get positions and genotypes
            pos = np.array([v.site.position for v in variants])
            gt = np.vstack([v.genotypes for v in variants])

            # Identify columns
            tgt_cols = list(range(len(tgt_samples)))
            ref_cols = list(range(len(tgt_samples), len(tgt_samples) + len(ref_samples)))

            # Filter positions
            if ref_cols:
                ref_gt = gt[:, ref_cols]
                keep_positions = ~np.any(ref_gt == 1, axis=1)
                gt = gt[keep_positions]
                pos = pos[keep_positions]

            if len(pos) == 0:
                continue

            # Calculate sstar and collect scores
            all_sstar_scores.extend(sstar(gt[:, tgt_cols], pos))

        # Calculate statistics
        if not all_sstar_scores:
            return (0.0, 0.0)
        return (np.mean(all_sstar_scores), np.std(all_sstar_scores))

    # Calculate statistics for all combinations and format results
    results = []
    for tgt, ref in combinations:
        mean, std = get_sstar_stats(tgt, ref)
        set_name = f"{tgt}-{ref}"
        results.append({"set": f"{set_name}-mean", "val": mean})
        results.append({"set": f"{set_name}-sd", "val": std})

    # Return as DataFrame
    return pd.DataFrame(results)
)"
  )
  reticulate::py$calc_sstar_summary(ts)
}


## load real data
tstest<-ts_read(file="~/demog-afr/preprocess/test_red.ts")
tsreal<-ts_read(file="/lisc/data/scratch/admixlab/mk_data/san/all_chroms_subset_WSS_subsel.trees")

load(file="~/demog-afr/preprocess/test_red.Robject")

metaset<-attr(tsave,"metadata")
metaset$sample_names<-names
metaset$subset_names<-names
metaset$sampling$name<-names
metaset$sampling$pop<-popv
metaset$arguments$SEQUENCE_LENGTH<-tsreal$sequence_length

attr(tsreal,"raw_individuals")$sampled<-TRUE
attr(tsreal,"nodes")$name<-NA
attr(tsreal,"nodes")$name<-c(hapsamp,rep(NA,nrow(attr(tsreal,"nodes"))-length(hapsamp)))
attr(tsreal,"nodes")$sampled<-c(rep(TRUE,length(hapsamp)),rep(FALSE,nrow(attr(tsreal,"nodes"))-length(hapsamp)))
attr(tsreal,"model")<-attr(tsave,"model")
attr(tsreal,"mutated")<-attr(tsave,"mutated")
attr(tsreal,"metadata")<-metaset
tsreal$nsamples<-metaset$sample_names

## REAL observed
observed<-list(diversity  = compute_diversity(tsreal),
               divergence  = compute_divergence(tsreal),
               f2  = compute_f2(tsreal),
               fst  = compute_fst(tsreal),
               td  = compute_td(tsreal),
               segs  = compute_segs(tsreal),
               afss  = compute_afss(tsreal),
               private  = compute_private(tsreal),
               sstar  = compute_sstar(tsreal) )

save(observed,file="/lisc/data/scratch/admixlab/mk_data/san/real_red_data_stats.Robject")
print("done")
print(Sys.time())

q()




