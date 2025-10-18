#! /bin/bash
#
##SBATCH --mem-per-cpu 12000
#SBATCH --mem-per-cpu 6000
#SBATCH -c 8
##SBATCH -t 1200:00
#SBATCH -t 900:00
#SBATCH -J fit_latent
##SBATCH -J lat_cv_80
##SBATCH -a 1-10
#SBATCH --mail-user=oward@sfu.ca
#SBATCH --mail-type=ALL

module load gcc r/4.5.0

echo "Launching R"
date

Rscript Fit_2015_Stan.R
#Rscript Null_CV.R
#Rscript Latent_CV.R

echo "Completed"
date

# end of script
