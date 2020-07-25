#!/bin/sh
#
#SBATCH -A stats
#SBATCH -J sims
#SBATCH -c 4
#SBATCH -t 320:00
#SBATCH --mem-per-cpu 12gb
##SBATCH -a 1-10
#SBATCH --mail-user=ogw2103@columbia.edu
#SBATCH --mail-type=ALL

module load R/4.0.1
echo "Launching R"
date

R CMD BATCH --no-save --vanilla agg_sims.R routput_agg_sims
echo "Completed"
date

#end of script
