# demog-afr
modeling for the KS paper

## Pre-processing real data

First, the VCF files were subset using `bcftools` and phased using `SHAPEIT5` (AR add).


Then, `tsinfer` and `tsdate` were applied to create tree sequence format data (JH add).


Filtering for high-quality sites across individuals was performed on the all-sites VCF files, subsetting to individuals of interest, using `bcftools`, as shown in the script `preprocessing/aks_filtering.sh`.


The 5000 best windows across the genome were then identified with the script `preprocessing/real_prep.R` and saved in `bed` format.


The tree sequences were then intersected with these `bed` files, concatenated across 22 autosomal chromosomes and saved in a single pseudo-continuous tree sequence object resembling the simulated data of 250 Mbp in `slendr`. This procedure is based on custom python code reported in `preprocessing/ts_modify_final.ipynb`. Finally, population and individual labels were updated to match the naming conventions in the simulations.


This object was loaded in `R` and summary statistics were calculated as for the simulated data, using the same functions and yielding the same order of statistics. Finally, the `slendr`-specific metadata was updated based on a simulated tree seuqence. These processes are shown in `preprocessing/real_data.R`. 


## ABC inference

We created ranges of priors starting from the parameters of the "B26_Deep_Merged_Half.yaml" model. Here, large ranges were given for effective population sizes, divergence times and gene flow proportions. Gene flow pulse times were kept from the original model to reduce the number of parameters. A model with ghost introgression was built using the `slendr` definitions and `demografr` framework. Then, the following summary statistics were calculated: diversity (windowed), pairwise divergence (windowed), pairwise f2 (genome-wide), pairwise FST (windowed), population-wise Tajima's D (windowed), allele-frequency-based statistics and pairwise correlation (genome-wide), relative pairwise private SNPs (windowed), pairwise S* values (windowed). `msprime` backward simulations of 250 Mbp (5000 windows) were performed in batches of 8. The script for these procedures is in `simulation/demog_sim.R`, with `simulation/demog_run.sh` to launch 3125 jobs (25000 iterations).


The inference of best parameters for this model was performed using the `demografr` framework. All simulated iterations were collected, merged with the statistics from real data, and a subset of statistics selected (real values not outside the distribution of simulations for all combinations, sufficient variance of simulated data), ending with 56 informative statistics. These were then normalized to values between 0 and 1 to provide input within the same ranges the ABC inference. ABC with neural networks (`tol=0.05, numnet=100`) was performed. These steps are provided in `abc_run.R`.


## Model selection

Model selection between ghost and no-ghost scenarios was performed as shown in `simulation/model_comp.R`. Here, we used as priors the posteriors of the ABC inference, where no-ghost scenarios simply do not contain gene flow from the ghost population. 20000 iterations of 100 Mbp were performed. Given the simplicity of this comparison, we used `tol=0.05` for cross-validation and estimating posterior model probabilities.


Model selection including WSS models is shown in `simulation/model_comp_wss.R`. Here, no RHG populations exist, and the analysis is restricted to Nama (representing Khoe-San) and MSL (representing RHGn populations). All parameters involving RHG populations were still included in the ghost/no-ghost models, only no sampling of these populations was performed. For priors of the WSS models, we used the posterior ranges from the original publication (Table S5 for continuous migration; Table S7 for migration pulses; Table S2 for fixed parameters; and the original `demes` files for missing parameters). Here, we performed 20000 iterations of 100 Mbp for each of the four models, and calculated the summary statistics available for a 2-population setting. Since some models were fitting poorly, we used `tol=0.1` for cross-validation and estimating posterior model probabilities.


Finally, we performed a model selection only based on three MSL and five Nama individuals, subsetting the original data accordingly, as using all Khoe-San and RHGn individuals is not exactly as presented in the original study. Adjustments of statistics were made according to the smaller number of individuals. Again, `tol=0.1` was used.


## Plots

Plots of diagnostic and validation steps are provided, as documented in `abc_mainfig.R`, and shown in `plots/`:

- `ABC_composite_figure.pdf`: composite figure for the main text
- `SI_errors.pdf`: prediction error for all parameters of the best model
- `SI_posteriors.pdf`: ABC prior and posterior distributions for all parameters of the best model



