iter=as.numeric(unlist(strsplit(as.character(unlist((commandArgs(TRUE)))),split=" ")))
print(iter)

library(demografr)
library(slendr)
init_env(uv=T)
.libPaths( c("/lisc/data/scratch/admixlab/mk_data/rlib/", .libPaths()) )  
#library(tibble, quietly = TRUE)
library(e1071, quietly = TRUE)
library(GenomicRanges, quietly = TRUE)
library(reticulate, quietly = TRUE)
py_require("scipy")
library(future)
plan(multisession, workers = 10)

## sequence and window length
#slen=10000000
winsize=50000

## names of combinations
idx_pairs <- combn(seq_along(1:4), 2, simplify = FALSE)
idx_nam<-paste(do.call(rbind,idx_pairs)[,1],do.call(rbind,idx_pairs)[,2],sep="-")
uninam=c("eRHG","KS","RHGn","wRHG")
combnam<-c("eRHG-KS","eRHG-RHGn","eRHG-wRHG","KS-RHGn","KS-wRHG","RHGn-wRGH")


## the model
model_ghost <- function(Ne_Ghost, Ne_KS, Ne_RuN, Ne_ARHG, Ne_RHGn,Ne_wRHG, Ne_eRHG,T_AGhost, T_KS, T_ARHG, T_RHG,gf_KR, gf_RK, gf_KN, gf_NK, gf_KA, gf_AK, gf_NA, gf_AN, gf_KN2, gf_NK2,gf_WE, gf_EW, gf_EN, gf_NE, gf_WN, gf_NW,gf_Ghost) {
  A <- population("A", time = 50001,    N = 18507,remove = T_KS-100)
  Ghost <- population("Ghost", time = T_AGhost, N = Ne_Ghost, parent = A)
  AH <- population("AH", time = T_AGhost, N = 18507, parent = A)
  KS <- population("KS", time = T_KS, N = Ne_KS, parent = AH)
  RuN <- population("ARHG_RHGn", time = T_KS, N = Ne_RuN, parent = AH)
  ARHG <- population("ARHG", time = T_ARHG, N = Ne_ARHG, parent = RuN)
  RHGn <- population("RHGn", time = T_ARHG, N = Ne_RHGn, parent = RuN)
  wRHG <- population("wRHG", time = T_RHG, N = Ne_wRHG, parent = ARHG)
  eRHG <- population("eRHG", time = T_RHG, N = Ne_eRHG, parent = ARHG)
  
  gf<- list( gene_flow(from = KS, to = RuN, start = 9679, end = 9678, proportion = gf_KR),
   gene_flow(from = RuN, to = KS, start = 9678, end = 9677, proportion = gf_RK),
   gene_flow(from = KS, to = RHGn, start = 4490, end = 4489, proportion = gf_KN),
   gene_flow(from = RHGn, to = KS, start = 4491, end = 4490, proportion = gf_NK),
   gene_flow(from = KS, to = ARHG, start = 4408, end = 4407, proportion = gf_KA),
   gene_flow(from = ARHG, to = KS, start = 4407, end = 4406, proportion = gf_AK),
   gene_flow(from = RHGn, to = ARHG, start = 2756, end = 2755, proportion = gf_NA),
   gene_flow(from = ARHG, to = RHGn, start = 2757, end = 2756, proportion = gf_AN),
   gene_flow(from = KS, to = RHGn, start = 1383, end = 1382, proportion = gf_KN2),
   gene_flow(from = RHGn, to = KS, start = 1382, end = 1381, proportion = gf_NK2),
   gene_flow(from = wRHG, to = eRHG, start = 411, end = 410, proportion = gf_WE),
   gene_flow(from = eRHG, to = wRHG, start = 410, end = 409, proportion = gf_EW),
   gene_flow(from = eRHG, to = RHGn, start = 393, end = 392, proportion = gf_EN),
   gene_flow(from = RHGn, to = eRHG, start = 394, end = 393, proportion = gf_NE),
   gene_flow(from = wRHG, to = RHGn, start = 302, end = 301, proportion = gf_WN),
   gene_flow(from = RHGn, to = wRHG, start = 303, end = 302, proportion = gf_NW),
   gene_flow(from = Ghost, to = KS, start = 303, end = 302, proportion = gf_Ghost))
  
model <- compile_model(
  populations = list(A, Ghost, AH, RuN, ARHG,KS, RHGn, wRHG, eRHG),
  gene_flow = gf,
  generation_time = 1)

samples <- schedule_sampling(
  model, times = 0,
  list(KS, 25), list(RHGn, 50), list(eRHG, 11), list(wRHG, 17),
  strict = TRUE)

  return(list(model, samples))
}


model_noghost <- function(Ne_KS, Ne_RuN, Ne_ARHG, Ne_RHGn,Ne_wRHG, Ne_eRHG,T_KS, T_ARHG, T_RHG,gf_KR, gf_RK, gf_KN, gf_NK, gf_KA, gf_AK, gf_NA, gf_AN, gf_KN2, gf_NK2,gf_WE, gf_EW, gf_EN, gf_NE, gf_WN, gf_NW) {
  A <- population("A", time = 50001,    N = 18507,remove = T_KS-100)
  Ghost <- population("Ghost", time = 20000, N = 100, parent = A)
  AH <- population("AH", time = 18508, N = 18507, parent = A)
  KS <- population("KS", time = T_KS, N = Ne_KS, parent = AH)
  RuN <- population("ARHG_RHGn", time = T_KS, N = Ne_RuN, parent = AH)
  ARHG <- population("ARHG", time = T_ARHG, N = Ne_ARHG, parent = RuN)
  RHGn <- population("RHGn", time = T_ARHG, N = Ne_RHGn, parent = RuN)
  wRHG <- population("wRHG", time = T_RHG, N = Ne_wRHG, parent = ARHG)
  eRHG <- population("eRHG", time = T_RHG, N = Ne_eRHG, parent = ARHG)
  
  gf<- list( gene_flow(from = KS, to = RuN, start = 9679, end = 9678, proportion = gf_KR),
             gene_flow(from = RuN, to = KS, start = 9678, end = 9677, proportion = gf_RK),
             gene_flow(from = KS, to = RHGn, start = 4490, end = 4489, proportion = gf_KN),
             gene_flow(from = RHGn, to = KS, start = 4491, end = 4490, proportion = gf_NK),
             gene_flow(from = KS, to = ARHG, start = 4408, end = 4407, proportion = gf_KA),
             gene_flow(from = ARHG, to = KS, start = 4407, end = 4406, proportion = gf_AK),
             gene_flow(from = RHGn, to = ARHG, start = 2756, end = 2755, proportion = gf_NA),
             gene_flow(from = ARHG, to = RHGn, start = 2757, end = 2756, proportion = gf_AN),
             gene_flow(from = KS, to = RHGn, start = 1383, end = 1382, proportion = gf_KN2),
             gene_flow(from = RHGn, to = KS, start = 1382, end = 1381, proportion = gf_NK2),
             gene_flow(from = wRHG, to = eRHG, start = 411, end = 410, proportion = gf_WE),
             gene_flow(from = eRHG, to = wRHG, start = 410, end = 409, proportion = gf_EW),
             gene_flow(from = eRHG, to = RHGn, start = 393, end = 392, proportion = gf_EN),
             gene_flow(from = RHGn, to = eRHG, start = 394, end = 393, proportion = gf_NE),
             gene_flow(from = wRHG, to = RHGn, start = 302, end = 301, proportion = gf_WN),
             gene_flow(from = RHGn, to = wRHG, start = 303, end = 302, proportion = gf_NW) )
  
  model <- compile_model(
    populations = list(A, Ghost, AH, RuN, ARHG,KS, RHGn, wRHG, eRHG),
    gene_flow = gf,
    generation_time = 1)
  
  samples <- schedule_sampling(
    model, times = 0,
    list(KS, 25), list(RHGn, 50), list(eRHG, 11), list(wRHG, 17),
    strict = TRUE)
  
  return(list(model, samples))
}


## the priors  
## updated with ABC inference
priors_ghost <- list(
  Ne_Ghost  ~ runif(1675.419,  93944.04),  Ne_KS  ~ runif(55100.99, 97401.38),
  Ne_RuN  ~ runif(20921.35, 66213.48),  Ne_ARHG  ~ runif( 18934.54, 61591.38),
  Ne_RHGn  ~ runif( 1830.665, 20672.36),  Ne_wRHG  ~ runif(4815.942, 12309.55),
  Ne_eRHG  ~ runif(8732.584, 19156.39),  
  
  T_AGhost  ~ runif(37108.11,    57766.51),  T_KS  ~ runif(11426.21, 15007.2),
  T_ARHG  ~ runif(4284.651, 8549.017),  T_RHG  ~ runif( 223.4391, 1121.783),

  gf_KR ~ runif(0.0502588, 0.4196633),   gf_RK ~ runif(0, 0.280166),
  gf_KN ~ runif(0, 0.3396901),   gf_NK ~ runif( 0.07924232, 0.436566),
  gf_KA ~ runif(0, 0.2951918),  gf_AK ~ runif(0.01832189 , 0.3903628),
  gf_NA ~ runif(0.07977693, 0.3472172),  gf_AN ~ runif(0, 0.2392889),
  gf_KN2 ~ runif(0, 0.2152445), gf_NK2 ~ runif(0.3329192, 0.5162124),
  gf_WE ~ runif(0.09450591, 0.4730481), gf_EW ~ runif(0.0427776, 0.4205468),
  gf_EN ~ runif(0.08667002, 0.4402068),  gf_NE ~ runif(0.05114712, 0.3261772),
  gf_WN ~ runif(0.01925199 , 0.3283356),   gf_NW ~ runif(0.2205767, 0.5631253),
  gf_Ghost ~ runif(0.01016583, 0.05817588)
)



## priors for no ghost model (updated with ABC inference)
priors_noghost <- list(
  Ne_KS  ~ runif(55100.99, 97401.38),
  Ne_RuN  ~ runif(20921.35, 66213.48),  Ne_ARHG  ~ runif( 18934.54, 61591.38),
  Ne_RHGn  ~ runif( 1830.665, 20672.36),  Ne_wRHG  ~ runif(4815.942, 12309.55),
  Ne_eRHG  ~ runif(8732.584, 19156.39),  
  
  T_KS  ~ runif(11426.21, 15007.2),
  T_ARHG  ~ runif(4284.651, 8549.017),  T_RHG  ~ runif( 223.4391, 1121.783),
  
  gf_KR ~ runif(0.0502588, 0.4196633),   gf_RK ~ runif(0, 0.280166),
  gf_KN ~ runif(0, 0.3396901),   gf_NK ~ runif( 0.07924232, 0.436566),
  gf_KA ~ runif(0, 0.2951918),  gf_AK ~ runif(0.01832189 , 0.3903628),
  gf_NA ~ runif(0.07977693, 0.3472172),  gf_AN ~ runif(0, 0.2392889),
  gf_KN2 ~ runif(0, 0.2152445), gf_NK2 ~ runif(0.3329192, 0.5162124),
  gf_WE ~ runif(0.09450591, 0.4730481), gf_EW ~ runif(0.0427776, 0.4205468),
  gf_EN ~ runif(0.08667002, 0.4402068),  gf_NE ~ runif(0.05114712, 0.3261772),
  gf_WN ~ runif(0.01925199 , 0.3283356),   gf_NW ~ runif(0.2205767, 0.5631253)
)



## the statistics
# nucleotide diversity in each population in windows (median & mad)
compute_diversity <- function(ts) {
  samples <- ts_names(ts, split = "pop")
  slen<-ts$sequence_length
  winsize=50000
  tm<-ts_diversity(ts, sample_sets = samples,windows=seq(winsize,slen-winsize,winsize))
  tm$med<-unlist(lapply(tm$diversity,median))*10000
  tm$mad<-unlist(lapply(tm$diversity,mad))*10000
  tm<-data.frame(set=paste(rep(tm$set,2),c(rep("med",4),rep("mad",4)),sep="_"),val=c(tm$med,tm$mad))
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
  tm<-data.frame(set=paste(rep(paste(tm$x,tm$y,sep="-"),2),c(rep("med",6),rep("mad",6)),sep="_"),val=c(tm$med,tm$mad))
  return(tm)
}

# pairwise f2-statistics genome-wide
compute_f2 <- function(ts) {
  calcf2<-function(pops,ts) { ts_f2(ts,A=pops[[1]], B=pops[[2]]) }
  samples <- ts_names(ts, split = "pop")
  combnam<-c("eRHG-KS","eRHG-RHGn","eRHG-wRHG","KS-RHGn","KS-wRHG","RHGn-wRGH")
  idx_pairs <- combn(seq_along(1:4), 2, simplify = FALSE)
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
  tm$med<-unlist(lapply(tm$Fst,median))*100
  tm$mad<-unlist(lapply(tm$Fst,mad))*100
  tm<-data.frame(set=paste(rep(paste(tm$x,tm$y,sep="-"),2),c(rep("med",6),rep("mad",6)),sep="_"),val=c(tm$med,tm$mad))
  return(tm)
}

# population-wise Tajima's D in windows (median & mad)
compute_td <- function(ts) {
  samples <- ts_names(ts, split = "pop")
  slen<-ts$sequence_length
  winsize=50000
  tm<-ts_tajima(ts, sample_sets = samples,windows=seq(winsize,slen-winsize,winsize))
  tm$med<-unlist(lapply(tm$D,median))
  tm$mad<-unlist(lapply(tm$D,mad))
  tm<-data.frame(set=paste(rep(tm$set,2),c(rep("med",4),rep("mad",4)),sep="_"),val=c(tm$med,tm$mad))
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
  tm<-data.frame(set=paste(rep(tm$set,2),c(rep("med",4),rep("mad",4)),sep="_"),val=c(tm$med,tm$mad))
  return(tm)
}

# AFS-based statistics per population, in subsets of 11: median & mad of AFS, pairwise correlation between populations, alleles at frequency 1 per population, fixed alleles per population
# 11 because lowest number of individuals (eRHG)
compute_afss <- function(ts) {
  library(e1071, quietly = TRUE)
  pstat<-function(input) { c(median(input)*1000,mad(input)*1000, kurtosis(input)) }
  bstat<-function(pr,afs) { cov(afs[[pr[1]]],afs[[pr[2]]],method="spearman") }
  samples <- ts_names(ts, split = "pop")
  smples <- lapply(samples, head, 11)
  names(smples)<-names(samples)
  afsr<-list(); afsl<-list()
  for (i in (1:length(smples))) {afsr[[i]]<-ts_afs(ts,sample_sets=list(smples[[i]])); afsl[[i]]<-afsr[[i]]/sum(afsr[[i]]) }
  names(afsl)<-names(samples)
  perp<-do.call(rbind,lapply(afsl,pstat))  
  colnames(perp)<-c("med","mad","kurt")
  idx_pairs <- combn(seq_along(1:4), 2, simplify = FALSE)
  idx_nam<-paste(do.call(rbind,idx_pairs)[,1],do.call(rbind,idx_pairs)[,2],sep="-")
  pcov<-unlist(lapply(idx_pairs,bstat,afs=afsl))
  names(pcov)<-idx_nam
  afspe<-do.call(cbind,afsr)[c(2,length(afsr[[1]])),]/c(1000,100)
  rownames(afspe)<-c("singleton","fixed")
  colnames(afspe)<-names(samples)
  retva<-unlist(list(perp,pcov,afspe))
  retva<-data.frame(names=c(paste("med_",names(samples),sep=""),paste("mad_",names(samples),sep=""),paste("kurt_",1:4,sep=""),paste("cov_",1:6,sep=""),paste("singl_",names(samples),sep=""),paste("fixed_",names(samples),sep="")),value=retva)
  return(retva)
}

# mean and sd private SNPs for KSvsRHGn wRHGvsRHGn eRHGvsRHGn RHGnvsKS in windows
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
  combinset<-list(list(2,3),list(4,3),list(1,3),list(3,2))
#  tsm<-ts_mutate(ts,mutation_rate=1e-8)
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
    samples = ts.metadata["slendr"]["sample_names"]

    # Define all target-reference combinations to analyze
    combinations = [
        ("KS", "RHGn"),
        ("wRHG", "RHGn"),
        ("eRHG", "RHGn"),
        ("RHGn", "KS")
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


functions <- list(
  diversity  = compute_diversity,
  divergence  = compute_divergence,
  f2  = compute_f2,
  fst  = compute_fst,
  td  = compute_td,
  segs  = compute_segs,
  afss  = compute_afss,
  private  = compute_private,
  sstar = compute_sstar
)

## observed
load(file="/lisc/data/scratch/admixlab/mk_data/san/real_data_stats.Robject")

#validate_abc(model_ghost, priors, functions, observed,
#             sequence_length = 500000, recombination_rate = 1e-8)


## simulate both models

data_ghost <- simulate_abc(
  model_ghost, priors_ghost, functions, observed, iterations = 10,
  sequence_length = 100000000, recombination_rate = 1e-8, mutation_rate = 1.29e-8
)

print("ghost")

data_noghost <- simulate_abc(
  model_noghost, priors_noghost, functions, observed, iterations = 10,
  sequence_length = 100000000, recombination_rate = 1e-8, mutation_rate = 1.29e-8
)
print("noghost")

save(data_ghost,data_noghost,file=paste("/lisc/data/scratch/admixlab/mk_data/san/simul_nog/simul_nog_",iter,".Robject",sep=""))

print("done")
print(Sys.time())

q()



library("dplyr")
library(scales)
library(abc)

load(file="/lisc/data/scratch/admixlab/mk_data/san/real_data_stats.Robject")

allgo<-list();allno<-list()
it<-1; mis=c()
for (iter in c(1:2000)) {
  tt<-try(load(file=paste("/lisc/data/scratch/admixlab/mk_data/san/simul_nog/simul_nog_",iter,".Robject",sep="")))
  if (inherits(tt,"try-error")) { mis=c(mis,iter);next }
  data_ghost$observed<-observed
  data_noghost$observed<-observed
  allgo[[it]]<-data_ghost
  allno[[it]]<-data_noghost
  it=it+1
}

allgo<-combine_data(allgo)
allno<-combine_data(allno)

# select stats that make sense & normalize
input1<-allgo
input2<-allno
  sinp<-input1
  observed<-input1$observed
  patrn<-c("diversity_.*_med","divergence_.*_med","f2","fst_.*_med","td_.*_med","segs_.*_med","afss_med_|afss_mad_|afss_cov_|afss_fixed_","private.*_mean","sstar_.*-mean")
  selc<-list();slc<-list()
  for (j in (1:length(observed))) { selc[[j]]<-colnames(input1$simulated)[grep(patrn[j],colnames(input1$simulated))];slc[[j]]<-grep(patrn[j],colnames(input1$simulated)) }

  nobserved<-list()
  alcol<-colnames(input1$simulated)
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
  normvals <- as.data.frame(rbind(input1$simulated[,unlist(selc)],input2$simulated[,unlist(selc)],obsv)) %>%
    mutate(across(where(is.numeric), ~ rescale(.x, to = c(0, 1), na.rm = TRUE)))

  simnorm1<-normvals[-nrow(normvals),][1:nrow(input1$simulated),]
  simnorm2<-normvals[-nrow(normvals),][-c(1:nrow(input1$simulated)),]
  obsv<-unlist(normvals[nrow(normvals),])

  fosn<-list(data.frame(set=osn[1:4],val=obsv[1:4]),data.frame(set=osn[5:10],val=obsv[5:10]), data.frame(pop=osn[11:16],f2=obsv[11:16]),
             data.frame(set=osn[17:22],val=obsv[17:22]),data.frame(set=osn[23:26],val=obsv[23:26]),
             data.frame(set=osn[27:30],val=obsv[27:30]),data.frame(names=osn[31:48],value=obsv[31:48]),
             data.frame(set=osn[49:52],val=obsv[49:52]),data.frame(set=osn[53:56],val=obsv[53:56]))
  names(fosn)<-names(input1$observed)

  input1$observed<-fosn
  input1$simulated<-simnorm1
  input2$observed<-fosn
  input2$simulated<-simnorm2

#gg_norm<-normfun(allgo,allno)
allgo_norm<-input1
allno_norm<-input2

abcG <- run_abc(allgo_norm, engine = "abc", tol = 0.05, method = "neuralnet",numnet=100)
abcN <- run_abc(allno_norm, engine = "abc", tol = 0.05, method = "neuralnet",numnet=100)

## the cross validation needed adjustment using the original abc package
#models <- list(abcG, abcN)
sumstat<-rbind(abcG$ss,abcN$ss)
colnames(sumstat)<-gsub("-","_",colnames(sumstat))
index<-c(rep("Ghost",1000),rep("NoGhost",1000))

cv_sel<-cv4postpr(index, sumstat, nval=100, tols=0.05, method="neuralnet")

#plot(cv_sel)
summary(cv_sel)

asumstat<-as.data.frame(rbind(allgo_norm$simulated,allno_norm$simulated))
colnames(asumstat)<-gsub("-","_",colnames(asumstat))
aindex<-c(rep("Ghost",20000),rep("NoGhost",20000))
target<-allgo_norm$observed[[1]][,2];tn<-allgo_norm$observed[[1]][,1]
for (j in (2:length(allgo_norm$observed))) { target<-c(target,allgo_norm$observed[[j]][,2]);tn=c(tn,allgo_norm$observed[[j]][,1])}


post_sel<-postpr(target, index=aindex, sumstat=asumstat, tol=0.05, method="neuralnet",nval=100)
summary(post_sel)

save(cv_sel,post_sel,file="/lisc/data/scratch/admixlab/mk_data/san/simul_nog/nog_selection.Robject")


