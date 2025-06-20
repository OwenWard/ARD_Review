#! /bin/bash
#
#SBATCH --mem-per-cpu 12000
#SBATCH -c 4
#SBATCH -t 3000:00
#SBATCH --mail-user=oward@sfu.ca
##SBATCH --mail-type=ALL

module load gcc r/4.2.2

echo "Launching R"
date

Rscript Fit_2015_Stan.R

echo "Completed"
date

# end of script
