#!/bin/sh
#
#SBATCH -A stats
#SBATCH -J exp_dt
#SBATCH -c 4
#SBATCH -t 600:00
#SBATCH --mem-per-cpu 12gb
##SBATCH -a 1
#SBATCH --mail-user=ogw2103@columbia.edu
#SBATCH --mail-type=ALL
##SBATCH -o cmmhp-%A.%a.out
##SBATCH --error=cmmhp-%A.%a.err 

module load R/4.0.1
echo "Launching R"
date

##R CMD BATCH --no-save --vanilla diagnostics.R routput_diag_1_$SLURM_ARRAY_TASK_ID
##R CMD BATCH --no-save --vanilla Exp_1.R routput_exp_1_long
R CMD BATCH --no-save --vanilla Exp_2.R routput_exp_dt
echo "Completed"
date

#end of script
