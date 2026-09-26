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
plan(multisession, workers = 8) # use 8 CPUs

## sequence and window length
#slen=10000000
winsize=50000

## names of combinations
idx_pairs <- combn(seq_along(1:4), 2, simplify = FALSE)
idx_nam<-paste(do.call(rbind,idx_pairs)[,1],do.call(rbind,idx_pairs)[,2],sep="-")
uninam=c("eRHG","KS","RHGn","wRHG")
combnam<-c("eRHG-KS","eRHG-RHGn","eRHG-wRHG","KS-RHGn","KS-wRHG","RHGn-wRGH")


## the model
model <- function(Ne_Ghost, Ne_KS, Ne_RuN, Ne_ARHG, Ne_RHGn,Ne_wRHG, Ne_eRHG,T_AGhost, T_KS, T_ARHG, T_RHG,gf_KR, gf_RK, gf_KN, gf_NK, gf_KA, gf_AK, gf_NA, gf_AN, gf_KN2, gf_NK2,gf_WE, gf_EW, gf_EN, gf_NE, gf_WN, gf_NW,gf_Ghost) {
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

## writing test data
#model2<-model(Ne_Ghost=1, Ne_KS=1, Ne_RuN=1, Ne_ARHG=1, Ne_RHGn=1,Ne_wRHG=1, Ne_eRHG=1,T_AGhost=50000, T_KS=40000, T_ARHG=30000, T_RHG=20000,gf_KR=0, gf_RK=0, gf_KN=0, gf_NK=0, gf_KA=0, gf_AK=0, gf_NA=0, gf_AN=0, gf_KN2=0, gf_NK2=0,gf_WE=0, gf_EW=0, gf_EN=0, gf_NE=0, gf_WN=0, gf_NW=0,gf_Ghost=0)
#slen=100000;ts <- simulate_model(model, priors, sequence_length = slen, recombination_rate = 1e-8,mutation_rate=1.29e-8)
#tsave=ts;save(tsave,file="~/demog-afr/preprocess/test.Robject");ts_write(tsave,file="~/demog-afr/preprocess/test.ts")


## the priors  
priors <- list(
  Ne_Ghost  ~ runif(1000,  100000),
  Ne_KS  ~ runif(86572-25000, 86572+25000),
  Ne_RuN  ~ runif(46968-25000, 46968+25000),
  Ne_ARHG  ~ runif(57397-25000, 57397+25000),
  Ne_RHGn  ~ runif(18022-10000, 18022+10000),
  Ne_wRHG  ~ runif(9461-5000, 9461+5000),
  Ne_eRHG  ~ runif(12592-6000, 12592+6000),
  
  T_AGhost  ~ runif(20000,    60000),
  T_KS  ~ runif(9679, 12117+2500),
  T_ARHG  ~ runif(4409, 9600),
  T_RHG  ~ runif(412, 892+500),

  gf_KR ~ runif(0, 0.4),
  gf_RK ~ runif(0, 0.3),
  gf_KN ~ runif(0, 0.4),
  gf_NK ~ runif(0, 0.4),
  gf_KA ~ runif(0, 0.3),
  gf_AK ~ runif(0, 0.4),
  gf_NA ~ runif(0, 0.3),
  gf_AN ~ runif(0, 0.3),
  gf_KN2 ~ runif(0, 0.3),
  gf_NK2 ~ runif(0, 0.3),
  gf_WE ~ runif(0, 0.4),
  gf_EW ~ runif(0, 0.4),
  gf_EN ~ runif(0, 0.4),
  gf_NE ~ runif(0, 0.3),
  gf_WN ~ runif(0, 0.3),
  gf_NW ~ runif(0, 0.4),
  gf_Ghost ~ runif(0, 0.2)
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

## example observed
# to be replaced with real data later
observed<-list(diversity=data.frame(set=paste(rep(uninam,2),c(rep("med",4),rep("mad",4)),sep="_"),val=c(3.1,4.4,3.2,7.1,0.2,0.4,0.3,0.5)),
               divergence=data.frame(set=paste(rep(combnam,2),c(rep("med",6),rep("mad",6)),sep="_"),val=c(9.1,7.2,5.1,9.1,9.5,7.1,0.01,0.04,0.3,0.3,0.2,0.1)),
               f2=data.frame(pop=combnam,f2=c(0.7,0.8,0.4,0.3,0.5,0.7)),
               fst=data.frame(set=paste(rep(combnam,2),c(rep("med",6),rep("mad",6)),sep="_"),val=c(3.7,4.1,1.3,3.1,4.2,4.1,1.1,1.3,0.9,0.8,0.6,1.3)),
               td=data.frame(set=paste(rep(uninam,2),c(rep("med",4),rep("mad",4)),sep="_"),val=c(-0.3,-1.1,0.1,-0.22,0.2,0.3,0.4,0.1)),
               segs=data.frame(set=paste(rep(uninam,2),c(rep("med",4),rep("mad",4)),sep="_"),val=c(1.2,3.1,2.1,3.22,0.2,0.3,0.4,0.1)),
               afss=data.frame(names=c(paste("med_",uninam,sep=""),paste("mad_",uninam,sep=""),paste("kurt_",1:4,sep=""),paste("cov_",1:6,sep=""),paste("singl_",uninam,sep=""),paste("fixed_",uninam,sep="")),value=unlist(list(data.frame(med=c(7.5,4.5,6.8,7.5),mad=c(3.6,3.2,5.5,6.5),kurt=c(12.1,11.3,14.3,12.3)),c(36.2,39.5,38.4,35.1,37.1,33.1),t(data.frame(singleton=c(15.1,20.1,13.1,14.1),fixed=c(10.1,5.3,9.1,8.7)))))),
               private=data.frame(set=paste(c(rep("KS",2),rep("wRHG",2),rep("eRHG",2),rep("RHGn",2)),rep(c("mean","sd"),4),sep="_"),val=c(5.5,8.1,7.3,4.1,6.4,5.4,7.1,5.5)),
               sstar=data.frame(set=c("KS-RHGn-mean","KS-RHGn-sd","wRHG-RHGn-mean","wRHG-RHGn-sd","eRHG-RHGn-mean","eRHG-RHGn-sd","RHGn-KS-mean","RHGn-KS-sd"),val=c(5400,3300,2500,4000,13000,12000,16000,13000))
                )
               

#validate_abc(model, priors, functions, observed,
#             sequence_length = slen, recombination_rate = 1e-8)


data <- simulate_abc(
  model, priors, functions, observed, iterations = 8,
  sequence_length = 250000000, recombination_rate = 1e-8, mutation_rate = 1.29e-8
)

save(data,file=paste("/lisc/data/scratch/admixlab/mk_data/san/simul/simul_",iter,".Robject",sep=""))

print("done")
print(Sys.time())

q()

