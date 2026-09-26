#!/bin/bash
#
#SBATCH --job-name=geno_proc
#SBATCH --cpus-per-task=8
#SBATCH --mem=25GB
#SBATCH --time=22:00:00
#SBATCH --output=./logs/filt_%j.out
#SBATCH --error=./logs/filt_%j.err
#SBATCH --export=NONE
#SBATCH --array=1-23

kdir=/lisc/data/scratch/admixlab/mk_data/tmp/ 
basdir=/lisc/data/scratch/admixlab/
wdir=$basdir/mk_data/san
module load HTSlib BCFtools 

ID=$SLURM_ARRAY_TASK_ID

###################################################

# normal: minimal coverage filter (all individuals need to have data), max. two alleles, no indels
echo "normal"
bcftools view -S $wdir/inds_to_subset.txt $kdir/25KS.48RHG.74comp.HCBP.$ID.recalSNP99.9.recalINDEL99.0.vcf.gz | bcftools annotate -x INFO | bcftools view -a | bcftools view -V indels | bcftools view -M2 | bcftools filter -e '(FMT/DP[*]<3 | FMT/DP[*]>250)' |  bcftools query -f "%POS\n" | bgzip > $wdir/$ID.len.txt.gz
# strict: stringent coverage, max. 10 individuals with missing genotypes (despite data is there)
echo "strict"
bcftools view -S $wdir/inds_to_subset.txt $kdir/25KS.48RHG.74comp.HCBP.$ID.recalSNP99.9.recalINDEL99.0.vcf.gz | bcftools annotate -x INFO | bcftools view -a | bcftools view -V indels | bcftools view -M2 | bcftools filter -e '(FMT/DP[*]<5 | FMT/DP[*]>150)' |  bcftools filter -i 'COUNT(FMT/GT=="./.")<10' |  bcftools query -f "%POS\n" | bgzip > $wdir/$ID.str.txt.gz
# very strict: stringent coverage, max. 10 individuals with missing genotypes (despite data is there), allele imbalance, mapability
 echo "very strict"
bcftools view -S $wdir/inds_to_subset.txt $kdir/25KS.48RHG.74comp.HCBP.$ID.recalSNP99.9.recalINDEL99.0.vcf.gz | bcftools annotate -x INFO | bcftools view -a | bcftools view -V indels | bcftools view -M2 | bcftools filter -e '(FMT/DP[*]<5 | FMT/DP[*]>150)' | bcftools filter -e '(FMT/AD[*:0]/FMT/DP[*]<0.15 | FMT/AD[*:1]/FMT/DP[*]<0.15) & FMT/GT[*]=="het"' |  bcftools filter -i 'COUNT(FMT/GT=="./.")<10' | java -jar ~/jvarkit.jar vcfbigwig -T mapability -B $wdir/../k36.Umap.MultiTrackMappability.bw  | egrep "#|mapability=1" | bcftools query -f "%POS\n" | bgzip > $wdir/$ID.ver.txt.gz
echo "done"

exit




