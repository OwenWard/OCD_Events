#### Experiment 2, Oct 26th
#### to investigate the relationship between 
#### the window size and the clustering result,
#### for a fixed data set and initial B matrix

.libPaths("/moto/stats/users/ogw2103/rpackages")

library(here)

source(here("Experiments/", "utils.R"))

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid

dT_vec <- seq(from = 0.1, to = 2, by = 0.1)

### then load in a dataset and an initial B? once which works 
### well for a known dT