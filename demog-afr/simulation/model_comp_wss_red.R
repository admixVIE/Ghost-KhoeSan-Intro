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
idx_pairs <- combn(seq_along(1:2), 2, simplify = FALSE)
idx_nam<-paste(do.call(rbind,idx_pairs)[,1],do.call(rbind,idx_pairs)[,2],sep="-")
uninam=c("Nama","MSL")
combnam<-c("Nama-MSL","MSL-Nama")

############################################
##### standard demography with ghost
############################################

## the model (ABC updated)
model_ghost <- function(Ne_Ghost, Ne_KS, Ne_RuN, Ne_ARHG, Ne_RHGn,Ne_wRHG, Ne_eRHG,T_AGhost, T_KS, T_ARHG, T_RHG,gf_KR, gf_RK, gf_KN, gf_NK, gf_KA, gf_AK, gf_NA, gf_AN, gf_KN2, gf_NK2,gf_WE, gf_EW, gf_EN, gf_NE, gf_WN, gf_NW,gf_Ghost) {
  A <- population("A", time = 50001,    N = 18507,remove = T_KS-100)
  Ghost <- population("Ghost", time = T_AGhost, N = Ne_Ghost, parent = A)
  AH <- population("AH", time = T_AGhost, N = 18507, parent = A)
  Nama <- population("Nama", time = T_KS, N = Ne_KS, parent = AH)
  RuN <- population("ARHG_RHGn", time = T_KS, N = Ne_RuN, parent = AH)
  ARHG <- population("ARHG", time = T_ARHG, N = Ne_ARHG, parent = RuN)
  MSL <- population("MSL", time = T_ARHG, N = Ne_RHGn, parent = RuN)
  wRHG <- population("wRHG", time = T_RHG, N = Ne_wRHG, parent = ARHG)
  eRHG <- population("eRHG", time = T_RHG, N = Ne_eRHG, parent = ARHG)
  
  gf<- list( gene_flow(from = Nama, to = RuN, start = 9679, end = 9678, proportion = gf_KR),
             gene_flow(from = RuN, to = Nama, start = 9678, end = 9677, proportion = gf_RK),
             gene_flow(from = Nama, to = MSL, start = 4490, end = 4489, proportion = gf_KN),
             gene_flow(from = MSL, to = Nama, start = 4491, end = 4490, proportion = gf_NK),
             gene_flow(from = Nama, to = ARHG, start = 4408, end = 4407, proportion = gf_KA),
             gene_flow(from = ARHG, to = Nama, start = 4407, end = 4406, proportion = gf_AK),
             gene_flow(from = MSL, to = ARHG, start = 2756, end = 2755, proportion = gf_NA),
             gene_flow(from = ARHG, to = MSL, start = 2757, end = 2756, proportion = gf_AN),
             gene_flow(from = Nama, to = MSL, start = 1383, end = 1382, proportion = gf_KN2),
             gene_flow(from = MSL, to = Nama, start = 1382, end = 1381, proportion = gf_NK2),
             gene_flow(from = wRHG, to = eRHG, start = 411, end = 410, proportion = gf_WE),
             gene_flow(from = eRHG, to = wRHG, start = 410, end = 409, proportion = gf_EW),
             gene_flow(from = eRHG, to = MSL, start = 393, end = 392, proportion = gf_EN),
             gene_flow(from = MSL, to = eRHG, start = 394, end = 393, proportion = gf_NE),
             gene_flow(from = wRHG, to = MSL, start = 302, end = 301, proportion = gf_WN),
             gene_flow(from = MSL, to = wRHG, start = 303, end = 302, proportion = gf_NW),
             gene_flow(from = Ghost, to = Nama, start = 303, end = 302, proportion = gf_Ghost))
  
  model <- compile_model(
    populations = list(A, Ghost, AH, RuN, ARHG,Nama, MSL, wRHG, eRHG),
    gene_flow = gf,
    generation_time = 1)
  
  samples <- schedule_sampling(
    model, times = 0,
    list(Nama, 5), list(MSL, 3),
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


## writing test data
#model2<-model_ghost(Ne_Ghost=1, Ne_KS=1, Ne_RuN=1, Ne_ARHG=1, Ne_RHGn=1,Ne_wRHG=1, Ne_eRHG=1,T_AGhost=50000, T_KS=40000, T_ARHG=30000, T_RHG=20000,gf_KR=0, gf_RK=0, gf_KN=0, gf_NK=0, gf_KA=0, gf_AK=0, gf_NA=0, gf_AN=0, gf_KN2=0, gf_NK2=0,gf_WE=0, gf_EW=0, gf_EN=0, gf_NE=0, gf_WN=0, gf_NW=0,gf_Ghost=0)
#slen=100000;ts <- simulate_model(model_ghost, priors_ghost, sequence_length = slen, recombination_rate = 1e-8,mutation_rate=1.29e-8)
#tsave=ts;save(tsave,file="~/demog-afr/preprocess/test_red.Robject");ts_write(tsave,file="~/demog-afr/preprocess/test_red.ts")


############################################
##### standard demography WITHOUT ghost
############################################

# the model
model_noghost <- function(Ne_KS, Ne_RuN, Ne_ARHG, Ne_RHGn,Ne_wRHG, Ne_eRHG,T_KS, T_ARHG, T_RHG,gf_KR, gf_RK, gf_KN, gf_NK, gf_KA, gf_AK, gf_NA, gf_AN, gf_KN2, gf_NK2,gf_WE, gf_EW, gf_EN, gf_NE, gf_WN, gf_NW) {
  A <- population("A", time = 50001,    N = 18507,remove = T_KS-100)
  Ghost <- population("Ghost", time = 18508, N = 100, parent = A)
  AH <- population("AH", time = 18508, N = 18507, parent = A)
  Nama <- population("Nama", time = T_KS, N = Ne_KS, parent = AH)
  RuN <- population("ARHG_RHGn", time = T_KS, N = Ne_RuN, parent = AH)
  ARHG <- population("ARHG", time = T_ARHG, N = Ne_ARHG, parent = RuN)
  MSL <- population("MSL", time = T_ARHG, N = Ne_RHGn, parent = RuN)
  wRHG <- population("wRHG", time = T_RHG, N = Ne_wRHG, parent = ARHG)
  eRHG <- population("eRHG", time = T_RHG, N = Ne_eRHG, parent = ARHG)
  
  gf<- list( gene_flow(from = Nama, to = RuN, start = 9679, end = 9678, proportion = gf_KR),
             gene_flow(from = RuN, to = Nama, start = 9678, end = 9677, proportion = gf_RK),
             gene_flow(from = Nama, to = MSL, start = 4490, end = 4489, proportion = gf_KN),
             gene_flow(from = MSL, to = Nama, start = 4491, end = 4490, proportion = gf_NK),
             gene_flow(from = Nama, to = ARHG, start = 4408, end = 4407, proportion = gf_KA),
             gene_flow(from = ARHG, to = Nama, start = 4407, end = 4406, proportion = gf_AK),
             gene_flow(from = MSL, to = ARHG, start = 2756, end = 2755, proportion = gf_NA),
             gene_flow(from = ARHG, to = MSL, start = 2757, end = 2756, proportion = gf_AN),
             gene_flow(from = Nama, to = MSL, start = 1383, end = 1382, proportion = gf_KN2),
             gene_flow(from = MSL, to = Nama, start = 1382, end = 1381, proportion = gf_NK2),
             gene_flow(from = wRHG, to = eRHG, start = 411, end = 410, proportion = gf_WE),
             gene_flow(from = eRHG, to = wRHG, start = 410, end = 409, proportion = gf_EW),
             gene_flow(from = eRHG, to = MSL, start = 393, end = 392, proportion = gf_EN),
             gene_flow(from = MSL, to = eRHG, start = 394, end = 393, proportion = gf_NE),
             gene_flow(from = wRHG, to = MSL, start = 302, end = 301, proportion = gf_WN),
             gene_flow(from = MSL, to = wRHG, start = 303, end = 302, proportion = gf_NW) )
  
  model <- compile_model(
    populations = list(A, Ghost, AH, RuN, ARHG,Nama, MSL, wRHG, eRHG),
    gene_flow = gf,
    generation_time = 1)
  
  samples <- schedule_sampling(
    model, times = 0,
    list(Nama, 5), list(MSL, 3),
    strict = TRUE)
  
  return(list(model, samples))
}

## the priors
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


############################################
############################################
# WSS models (converted from Demes YAML)
############################################
############################################

## function to convert migration rates to proportions to adhere to slendr definitions through msprime
### that is, since migration_rate is getting larger than 1, adjust by calculating the proportions first
confu <- function(m, L) -expm1(L * log1p(-m))

############################################
# WSS continuous-migration model
############################################

model_wss_cont <- function(t_anc_end,t_stem1_end, Nanc, Nstem1, Nstem2, Nvin, Ngbr, Nnam, Ngum, Nmsl, Nep, Nnam2, Ngum2, Nmsl2, Ngbr2, gf_stem, gf_nam, gf_msl, gf_gum, gf_msl2, gf_gum2, gf_gum3, gf_gbr,gf_ep,gf_ep2,gf_nam2,gf_nam3) {
  confu <- function(m, L) -expm1(L * log1p(-m))
  t_anc_end   <- t_anc_end/29 
  t_stem1_end <- t_stem1_end/29 
  t_stem1_V   <- 550000/29
  t_stem2_end <- 5000/29     
  t_vind_end  <- 50000/29    
  t_NI_start  <- 80000/29    
  t_NI_end    <- 45000/29    
  t_msl_start <- 60000/29    
  t_gbr_start <- 50000/29    
  t_ep_start  <- 12000/29    
  t_nama_shift <- 261/29     
  
  anc <- population("anc", time = t_anc_end+1, N = Nanc,  remove= t_anc_end-1)
  stem1 <- population("stem1", time = t_anc_end, N = Nanc, parent = anc , remove= t_stem1_end )
  stem2 <- population("stem2", time = t_anc_end, N = Nstem2, parent = anc, remove = t_stem2_end )
  Vindija <- population("Vindija", time = t_stem1_V, N = Nvin, parent = stem1, remove = t_vind_end )
  NI <- population("NI", time = t_NI_start, N = Nvin, parent = Vindija,remove = t_NI_end )
  Nama <- population("Nama", time = t_stem1_end+1, N = Nnam, parent = stem1)
  Gumuz <- population("Gumuz", time = t_stem1_end+1, N = Ngum, parent = stem1 )
  MSL <- population("MSL", time = t_msl_start, N = Nmsl, parent = Gumuz )
  GBR <- population("GBR", time = t_gbr_start, N = Ngbr, parent = Gumuz )
  EP <- population( "EP", time = t_ep_start, N = Nep, parent = Gumuz )
  
  stem1<-slendr::resize(stem1,"step", time = t_stem1_V, N = Nstem1)
  Nama<-slendr::resize(Nama,"step", time = t_nama_shift, N = Nnam2)
  Gumuz<-slendr::resize(Gumuz,"step", time = t_stem2_end, N = Ngum2)
  MSL<-slendr::resize(MSL, "step",  time = t_stem2_end, N = Nmsl2)
  GBR<-slendr::resize(GBR, Ngbr2, "exponential", t_vind_end-1, end = 0)

  gf <- list(
    gene_flow(from = stem1, to = stem2, start = t_anc_end, end = t_stem1_end,proportion = confu(gf_stem, t_stem1_V-t_stem1_end)),
    gene_flow(from = stem2, to = stem1, start = t_anc_end, end = t_stem1_end, proportion = confu(gf_stem, t_stem1_V-t_stem1_end)),
    gene_flow(from = stem2, to = Nama, start = t_stem1_end, end = t_stem2_end, proportion = confu(gf_nam, t_stem1_end-t_stem2_end)),
    gene_flow(from = Nama,  to = stem2, start = t_stem1_end, end = t_stem2_end, proportion = confu(gf_nam, t_stem1_end-t_stem2_end)),
    gene_flow(from = stem2, to = MSL, start = t_msl_start, end = t_stem2_end, proportion = confu(gf_msl, t_msl_start-t_stem2_end)),
    gene_flow(from = MSL,   to = stem2, start = t_msl_start, end = t_stem2_end, proportion = confu(gf_msl, t_msl_start-t_stem2_end)),
    gene_flow(from = stem2, to = Gumuz, start = t_stem1_end, end = t_stem2_end, proportion = confu(gf_gum, t_stem1_end-t_stem2_end)),
    gene_flow(from = Gumuz, to = stem2, start = t_stem1_end, end = t_stem2_end, proportion = confu(gf_gum, t_stem1_end-t_stem2_end)),
    gene_flow(from = stem2, to = EP, start = t_ep_start, end = t_stem2_end, proportion = confu(gf_gum, t_ep_start-t_stem2_end)),
    gene_flow(from = EP,    to = stem2, start = t_ep_start, end = t_stem2_end, proportion = confu(gf_gum, t_ep_start-t_stem2_end)),
    gene_flow(from = Nama, to = MSL, start = t_msl_start, end = 0, proportion = confu(gf_msl2, t_msl_start)),
    gene_flow(from = MSL,  to = Nama, start = t_msl_start, end = 0, proportion = confu(gf_msl2, t_msl_start)),
    gene_flow(from = Nama,  to = Gumuz, start = t_stem1_end, end = 0, proportion = confu(gf_gum2, t_stem1_end)),
    gene_flow(from = Gumuz, to = Nama,  start = t_stem1_end, end = 0, proportion = confu(gf_gum2, t_stem1_end)),
    gene_flow(from = Nama, to = EP, start = t_ep_start, end = 0,proportion = confu(gf_gum2, t_ep_start)),
    gene_flow(from = EP,   to = Nama, start = t_ep_start, end = 0, proportion = confu(gf_gum2, t_ep_start)),
    gene_flow(from = MSL,   to = Gumuz, start = t_msl_start, end = 0,proportion = confu(gf_gum3, t_msl_start)),
    gene_flow(from = Gumuz, to = MSL,   start = t_msl_start, end = 0,proportion = confu(gf_gum3, t_msl_start)),
    gene_flow(from = MSL, to = EP, start = t_ep_start, end = 0, proportion = confu(gf_gum3, t_ep_start)),
    gene_flow(from = EP,  to = MSL, start = t_ep_start, end = 0, proportion = confu(gf_gum3, t_ep_start)),
    gene_flow(from = Gumuz, to = GBR, start = t_gbr_start-1, end = 0,proportion = confu(gf_gbr, t_gbr_start-1)),
    gene_flow(from = GBR,   to = Gumuz, start = t_gbr_start-1, end = 0, proportion = confu(gf_gbr, t_gbr_start-1)),
    gene_flow(from = GBR, to = EP, start = t_ep_start, end = 0,proportion = confu(gf_gbr, t_ep_start)),
    gene_flow(from = EP,  to = GBR, start = t_ep_start, end = 0, proportion = confu(gf_gbr, t_ep_start)),
    gene_flow(from = Gumuz, to = EP, start = t_ep_start, end = 0, proportion = confu(gf_ep, t_ep_start)),
    gene_flow(from = EP,    to = Gumuz, start = t_ep_start, end = 0, proportion = confu(gf_ep, t_ep_start)) )
  
  gf <- c(gf,list(
            gene_flow(from = NI, to = GBR, start = t_NI_end+1, end = t_NI_end, proportion = 0.015),
            gene_flow(from = EP, to = Nama, start = 2000/29, end = (2000/29) - 1, proportion = gf_nam2),
            gene_flow(from = GBR, to = Gumuz, start = t_ep_start, end = t_ep_start - 1, proportion = gf_ep2),
            gene_flow(from = GBR, to = Nama, start = 290/29, end = (290/29) - 1, proportion = gf_nam3))
            )
  
  model <- compile_model(populations = list(anc, stem1, stem2, Vindija, NI, Nama, Gumuz, MSL, GBR, EP),gene_flow = gf,generation_time = 1)
  
  samples <- schedule_sampling(model, times = 0,
                               list(Nama, 5), list(MSL, 3),
                               strict = TRUE)
  
  return(list(model, samples))
}

#plot_model(model,log=T)

## priors for WSS continuous
priors_wss_cont <- list(
  t_anc_end ~ runif(1199000, 1382000),t_stem1_end ~ runif(127800,139900),Nanc ~ runif(4180, 6530),Nstem1 ~ runif(7130, 8400),
  Nstem2 ~ runif(11330, 13580),Nvin ~ runif(2380,2750),Nnam ~ runif(11160, 12390),Ngum ~ runif(3560,3810),
  Nmsl ~ runif(9500,10210),Ngbr ~ runif(940,981),Nep ~ runif(12800,13290),Nnam2 ~ runif(215,226),
  Ngum2 ~ runif(3560,3810),Nmsl2 ~ runif(25930,28930),Ngbr2 ~ runif(11490,12090),
  
  gf_stem ~ runif(6.14e-05,  7.83e-05),gf_nam ~ runif(4.5e-05,  7.04e-05),gf_msl ~ runif(13.8e-05,  18.9e-05),gf_gum ~ runif(2.18e-05,  3.65e-05),
  gf_msl2 ~ runif(0.41e-05,  1.43e-05),gf_gum2 ~ runif(3.92e-05,  4.56e-05),gf_gum3 ~ runif(21.2e-05,  21.8e-05),gf_gbr ~ runif(3.77e-05,  4.33e-05),
  gf_ep ~ runif(33e-05,  36.3e-05),gf_ep2 ~ runif(0.639,  0.656),gf_nam2 ~ runif(0.245,  0.263),gf_nam3 ~ runif(0.151,  0.160)
)


#data2 <- simulate_abc(model_wss_cont, priors_wss_cont, functions, observed, iterations = 8,sequence_length = 250000, recombination_rate = 1e-8, mutation_rate = 1.29e-8 )
  
############################################
##### WSS with pulses model
############################################

model_wss_puls <- function(t_anc_end,t_stem1_shrink_end, Nanc, Nstem1shrink, t_nama_start, t_gumuz_start, t_stem2_end, Nstem2, Nvin, Nnam, Ngum, Nea,Nmsl, Ngbr, Nep, Nnam2, Ngum2, Nmsl2, Ngbr2, gf_stem, gf_nam_msl, gf_nam_gum, gf_msl_ep, gf_gum_gbr,gf_gum_ep, gf_gum_merge, gf_nama_merge,  gf_ep2, gf_msl_pulse, gf_nam2, gf_nam3) {
  confu <- function(m, L) -expm1(L * log1p(-m))
  t_anc_end           <- t_anc_end / 29
  t_stem1_shrink_end  <- t_stem1_shrink_end / 29
  t_stem2_end         <- t_stem2_end / 29
  t_stem1_V           <- 550000 / 29
  t_vind_end          <- 50000 / 29
  t_NI_start          <- 80000 / 29
  t_NI_end            <- 45000 / 29
  t_msl_start         <- 60000 / 29
  t_gbr_start         <- 50000 / 29
  t_ep_start          <- 12000 / 29
  t_gumuz_start       <- t_gumuz_start / 29
  t_nama_start        <- t_nama_start / 29
  t_gum_shift_end     <- 5000 / 29
  t_nama_shift        <- 261 / 29
  Nstem1E <- 9104.777532154083
  Nstem1S <- 13436.752132969787

  anc   <- population("anc",   time = t_anc_end + 1, N = Nanc, remove = t_anc_end - 1)
  stem1 <- population("stem1", time = t_anc_end, N = Nanc, parent = anc, remove = t_stem1_shrink_end)
  stem2 <- population("stem2", time = t_anc_end,  N = Nstem2, parent = anc, remove = t_stem2_end)
  stem1E <- population("stem1E", time = t_stem1_V, N = Nstem1E, parent = stem1, remove = gf_gum_merge)
  stem1S <- population("stem1S", time = t_stem1_V,  N = Nstem1S, parent = stem1, remove = gf_nama_merge)
  Vindija <- population("Vindija", time = t_stem1_V - 1, N = Nvin, parent = stem1, remove = t_vind_end)
  NI      <- population("NI",      time = t_NI_start,    N = Nvin, parent = Vindija, remove = t_NI_end)
  Gumuz <- population("Gumuz", time = t_gumuz_start, N = Ngum, parent = stem2)
  Nama  <- population("Nama",  time = t_nama_start,  N = Nnam, parent = stem2)
  MSL   <- population("MSL",   time = t_msl_start,   N = Nmsl, parent = Gumuz)
  GBR   <- population("GBR",   time = t_gbr_start,   N = Ngbr, parent = Gumuz)
  EP    <- population("EP",    time = t_ep_start,    N = Nep,  parent = Gumuz)
  
  stem1<-slendr::resize(stem1, "step", time = t_stem1_V, N = Nstem1shrink)
  Nama<-slendr::resize(Nama,  "step", time = t_nama_shift,       N = Nnam2)
  MSL<-slendr::resize(MSL,   "step", time = t_gum_shift_end,    N = Nmsl2)
  GBR<-slendr::resize(GBR,   Ngbr2,  "exponential",             t_vind_end - 1, end = 0)
  
  gf <- list(
    gene_flow(from = stem1, to = stem2, start = t_stem1_V,          end = t_stem1_shrink_end, proportion = confu(gf_stem, t_stem1_V - t_stem1_shrink_end)),
    gene_flow(from = stem2, to = stem1, start = t_stem1_V,          end = t_stem1_shrink_end,proportion = confu(gf_stem, t_stem1_V - t_stem1_shrink_end)),
    gene_flow(from = Nama,  to = MSL,   start = t_msl_start,        end = 0,proportion = confu(gf_nam_msl, t_msl_start)),
    gene_flow(from = MSL,   to = Nama,  start = t_msl_start,        end = 0,proportion = confu(gf_nam_msl, t_msl_start)),
    gene_flow(from = Nama,  to = Gumuz, start = t_gumuz_start,       end = 0,proportion = confu(gf_nam_gum, t_gumuz_start)),
    gene_flow(from = Gumuz, to = Nama,  start = t_gumuz_start,       end = 0,proportion = confu(gf_nam_gum, t_gumuz_start)),
    gene_flow(from = Nama,  to = EP,    start = t_ep_start,         end = 0,proportion = confu(gf_nam_gum, t_ep_start)),
    gene_flow(from = EP,    to = Nama,  start = t_ep_start,         end = 0,proportion = confu(gf_nam_gum, t_ep_start)),
    gene_flow(from = MSL,   to = Gumuz, start = t_msl_start,        end = 0,proportion = confu(gf_msl_ep, t_msl_start)),
    gene_flow(from = Gumuz, to = MSL,   start = t_msl_start,        end = 0,proportion = confu(gf_msl_ep, t_msl_start)),
    gene_flow(from = MSL,   to = EP,    start = t_ep_start,         end = 0,proportion = confu(gf_msl_ep, t_ep_start)),
    gene_flow(from = EP,    to = MSL,   start = t_ep_start,         end = 0,proportion = confu(gf_msl_ep, t_ep_start)),
    gene_flow(from = Gumuz, to = GBR,   start = t_gbr_start - 1,    end = 0,proportion = confu(gf_gum_gbr, t_gbr_start - 1)),
    gene_flow(from = GBR,   to = Gumuz, start = t_gbr_start - 1,    end = 0,proportion = confu(gf_gum_gbr, t_gbr_start - 1)),
    gene_flow(from = GBR,   to = EP,    start = t_ep_start,         end = 0,proportion = confu(gf_gum_gbr, t_ep_start)),
    gene_flow(from = EP,    to = GBR,   start = t_ep_start,         end = 0,proportion = confu(gf_gum_gbr, t_ep_start)),
    gene_flow(from = Gumuz, to = EP,    start = t_ep_start,         end = 0,proportion = confu(gf_gum_ep, t_ep_start)),
    gene_flow(from = EP,    to = Gumuz, start = t_ep_start,         end = 0,proportion = confu(gf_gum_ep, t_ep_start))
  )
  
  gf <- c(gf, list(
    gene_flow(from = stem1E, to = Gumuz, start = t_gumuz_start, end = t_gumuz_start - 1,proportion = gf_gum_merge),
    gene_flow(from = stem1S, to = Nama,  start = t_nama_start -1 ,  end = t_nama_start -2,proportion = gf_nama_merge),
    gene_flow(from = GBR,    to = EP,    start = t_ep_start,    end = t_ep_start - 1,proportion = gf_ep2),
    gene_flow(from = NI,     to = GBR,   start = t_NI_end + 1,  end = t_NI_end,proportion = 0.015),
    gene_flow(from = stem2,  to = MSL,   start = t_stem2_end+2,   end = t_stem2_end + 1,proportion = gf_msl_pulse),
    gene_flow(from = EP,     to = Nama,  start = 2000 / 29,     end = (2000 / 29) - 1,proportion = gf_nam2),
    gene_flow(from = GBR,    to = Nama,  start = 290 / 29,      end = (290 / 29) - 1,proportion = gf_nam3)
  ))
  
  model <- compile_model(
    populations = list(anc, stem1, stem2, stem1E, stem1S, Vindija, NI, Gumuz, Nama, MSL, GBR, EP), gene_flow = gf,
    generation_time = 1
  )
  
  samples <- schedule_sampling(
    model, times = 0,
    list(Nama, 5), list(MSL, 3),
    strict = TRUE
  )
  
  return(list(model, samples))
}

#plot_model(model,log=T)

## priors for WSS with pulses
priors_wss_puls <- list(
  t_anc_end ~ runif(695000,1785000),t_stem1_shrink_end ~ runif(275600,478300),Nanc ~ runif(6410,11910),Nstem1shrink ~ runif(100,851),
  t_nama_start ~ runif(101300,124600),t_gumuz_start ~ runif(95400,111000), t_stem2_end ~ runif(22210,26300),
  Nstem2 ~ runif(22130,26340),Nvin ~ runif(2100,2730),Nnam ~ runif(11500,13560),Ngum ~ runif(3270,3570),Nea ~ runif(7970,8860),
  Nmsl ~ runif(9690,11290),Ngbr ~ runif(928,1002),Nep ~ runif(12780,13500),Nnam2 ~ runif(214,231),
  Nmsl2 ~ runif(26820,32280),Ngbr2 ~ runif(11850,12600), gf_stem ~ runif(3.95e-05,  12.3e-05),
  
  gf_nam_msl ~ runif(0.66e-05,  1.57e-05), gf_nam_gum ~ runif(4.45e-05,  5.52e-05), gf_msl_ep ~ runif(20.0e-05,  21.7e-05),
  gf_gum_gbr ~ runif(3.25e-05,  4.21e-05),gf_gum_ep ~ runif(33.7e-05,  40.9e-05),

  gf_gum_merge ~ runif(0.636, 0.667), gf_nama_merge ~ runif(0.242, 0.306),  gf_ep2 ~ runif(0.636,0.667), gf_msl_pulse ~ runif(0.182,0.211),
  gf_nam2 ~ runif(0.237,0.267), gf_nam3 ~ runif(0.154,0.167)
)

#data2 <- simulate_abc(model_wss_puls, priors_wss_puls, functions, observed, iterations = 8,sequence_length = 250000, recombination_rate = 1e-8, mutation_rate = 1.29e-8 )



############################################
### FUNCTIONS
############################################

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
  tm$med<-unlist(lapply(tm$Fst,median))*100
  tm$mad<-unlist(lapply(tm$Fst,mad))*100
  tm<-data.frame(set=paste(rep(paste(tm$x,tm$y,sep="-"),2),c(rep("med",1),rep("mad",1)),sep="_"),val=c(tm$med,tm$mad))
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

# mean and sd private SNPs in windows
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

## load observed
load(file="/lisc/data/scratch/admixlab/mk_data/san/real_red_data_stats.Robject")

# simulate data for all 4 models
data_wss_cont <- simulate_abc(
  model_wss_cont, priors_wss_cont, functions, observed, iterations = 10,
  sequence_length = 250000000, recombination_rate = 1e-8, mutation_rate = 1.29e-8
)
print("wsscont")

data_wss_puls <- simulate_abc(
  model_wss_puls, priors_wss_puls, functions, observed, iterations = 10,
  sequence_length = 250000000, recombination_rate = 1e-8, mutation_rate = 1.29e-8
)
print("wsspuls")

data_ghost <- simulate_abc(
  model_ghost, priors_ghost, functions, observed, iterations = 10,
  sequence_length = 250000000, recombination_rate = 1e-8, mutation_rate = 1.29e-8
)
print("ghost")

data_noghost <- simulate_abc(
  model_noghost, priors_noghost, functions, observed, iterations = 10,
  sequence_length = 250000000, recombination_rate = 1e-8, mutation_rate = 1.29e-8
)
print("noghost")

save(data_ghost,data_noghost,data_wss_cont,data_wss_puls,file=paste("/lisc/data/scratch/admixlab/mk_data/san/simul_red/simul_red_",iter,".Robject",sep=""))

print("done")
print(Sys.time())

q()




load(file=paste("/lisc/data/scratch/admixlab/mk_data/san/simul_red/simul_red_",1,".Robject",sep=""))
obsnam<-c()
for (j in (1:length(data_ghost$observed))) { obsnam<-c(obsnam,paste(names(data_ghost$observed)[j],"_",data_ghost$observed[[j]][,1],sep="")) }
observed<-data_ghost$observed

allgo<-list();allno<-list()
wss1<-list();wss2<-list()
it<-1; mis=c()
for (iter in c(1:2500)) {
  tt<-try(load(file=paste("/lisc/data/scratch/admixlab/mk_data/san/simul_red/simul_red_",iter,".Robject",sep="")))
  if (inherits(tt,"try-error")) { mis=c(mis,iter);next }
  data_ghost$observed<-observed;colnames(data_ghost$simulated)<-obsnam
  data_noghost$observed<-observed;colnames(data_noghost$simulated)<-obsnam
  allgo[[it]]<-data_ghost
  allno[[it]]<-data_noghost
  data_wss_puls$observed<-observed;colnames(data_wss_puls$simulated)<-obsnam
  data_wss_cont$observed<-observed;colnames(data_wss_cont$simulated)<-obsnam
  wss1[[it]]<-data_wss_puls
  wss2[[it]]<-data_wss_cont
  it=it+1
}

allgo<-combine_data(allgo)
allno<-combine_data(allno)
wss1<-combine_data(wss1)
wss2<-combine_data(wss2)

# select stats that make sense & normalize
input1<-allgo
input2<-allno
input3<-wss1
input4<-wss2
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
## normalize everything at once for all scenarios
normvals <- as.data.frame(rbind(input1$simulated[,unlist(selc)],input2$simulated[,unlist(selc)],wss1$simulated[,unlist(selc)],wss2$simulated[,unlist(selc)],obsv)) %>%  mutate(across(where(is.numeric), ~ rescale(.x, to = c(0, 1), na.rm = TRUE)))

sval<-nrow(input1$simulated)
simnorm1<-normvals[-nrow(normvals),][1:sval,]
simnorm2<-normvals[-nrow(normvals),][c(1:sval)+sval,]
simnorm3<-normvals[-nrow(normvals),][c(1:sval)+(sval*2),]
simnorm4<-normvals[-nrow(normvals),][c(1:sval)+(sval*3),]
obsv<-unlist(normvals[nrow(normvals),])

fosn<-list(data.frame(set=osn[1:2],val=obsv[1:2]),data.frame(set=osn[3],val=obsv[3]), data.frame(pop=osn[4],f2=obsv[4]),
           data.frame(set=osn[5],val=obsv[5]),data.frame(set=osn[6:7],val=obsv[6:7]),
           data.frame(set=osn[8:9],val=obsv[8:9]),data.frame(names=osn[10:16],value=obsv[10:16]),
           data.frame(set=osn[17:18],val=obsv[17:18]),data.frame(set=osn[19:20],val=obsv[19:20]))
names(fosn)<-names(input1$observed)

input1$observed<-fosn
input1$simulated<-simnorm1
input2$observed<-fosn
input2$simulated<-simnorm2
wss1$observed<-fosn
wss1$simulated<-simnorm3
wss2$observed<-fosn
wss2$simulated<-simnorm4

# make sure no NAs are there(or remove them?)
allgo_norm<-input1
allgo_norm$simulated[is.na(allgo_norm$simulated)]<-0
allno_norm<-input2
allno_norm$simulated[is.na(allno_norm$simulated)]<-0
wss1_norm<-wss1
wss1_norm$simulated[is.na(wss1_norm$simulated)]<-0
wss2_norm<-wss2
wss2_norm$simulated[is.na(wss2_norm$simulated)]<-0

abcG <- run_abc(allgo_norm, engine = "abc", tol = 0.05, method = "neuralnet",numnet=100)
abcN <- run_abc(allno_norm, engine = "abc", tol = 0.05, method = "neuralnet",numnet=100)
abcP <- run_abc(wss1_norm, engine = "abc", tol = 0.05, method = "neuralnet",numnet=100)
abcC <- run_abc(wss2_norm, engine = "abc", tol = 0.05, method = "neuralnet",numnet=100)

## the cross validation needed adjustment using the original abc package
#models <- list(abcG, abcN)
library(abc)
sumstat<-rbind(abcG$ss,abcN$ss,abcP$ss,abcC$ss)
colnames(sumstat)<-gsub("-","_",colnames(sumstat))
index<-c(rep("Ghost",1250 ),rep("NoGhost",1250 ),rep("WSS_puls",1250 ),rep("WSS_cont",1250 ))

cv_sel<-cv4postpr(index, sumstat, nval=100, tols=0.1, method="neuralnet")

summary(cv_sel)
plot(cv_sel)

asumstat<-as.data.frame(rbind(allgo_norm$simulated,allno_norm$simulated,wss1_norm$simulated,wss2_norm$simulated))
#asumstat<-as.data.frame(rbind(allgo_norm$simulated,wss1_norm$simulated,wss2_norm$simulated))
colnames(asumstat)<-gsub("-","_",colnames(asumstat))
aindex<-c(rep("Ghost",25000 ),rep("NoGhost",25000 ),rep("WSS_puls",25000 ),rep("WSS_cont",25000 ))
target<-allgo_norm$observed[[1]][,2];tn<-allgo_norm$observed[[1]][,1]
for (j in (2:length(allgo_norm$observed))) { target<-c(target,allgo_norm$observed[[j]][,2]);tn=c(tn,allgo_norm$observed[[j]][,1])}


post_sel<-postpr(target, index=aindex, sumstat=asumstat, tol=0.1, method="neuralnet",nval=100)
summary(post_sel)

save(cv_sel,post_sel,abcG, abcN, abcP, abcC,file="/lisc/data/scratch/admixlab/mk_data/san/simul_red/red_selection.Robject")

