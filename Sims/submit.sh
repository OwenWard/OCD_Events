#!/bin/sh
#
#SBATCH -A tzsts
#SBATCH -J sims
#SBATCH -c 4
#SBATCH -t 2000:00
#SBATCH --mem-per-cpu 12gb
##SBATCH -a 1-10
#SBATCH --mail-user=ogw2103@columbia.edu
#SBATCH --mail-type=ALL

module load R
echo "Launching R"
date

R CMD BATCH --no-save --vanilla sim_fits.R routput_sims_long
echo "Completed"
date

#end of script
