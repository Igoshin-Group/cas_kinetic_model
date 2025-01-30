#!/bin/bash
#SBATCH --job-name=off-target-model-fit
#SBATCH --partition=commons
#SBATCH --account=commons
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=01:00:00
#SBATCH --mail-user=adh8@rice.edu
#SBATCH --mail-type=FAIL

# Available from environment: SLURM_ARRAY_TASK_ID = ID of current task

# Other variables:
# RANK = on-target model rank (which on-target parameters to use)
# SEED = which seed to use for the RNG

# purge automatically loaded modules
module purge

# load MATLAB module
module load MATLAB/2021a

# change directory to code directory
cd $SHARED_SCRATCH/adh8/cas_kinetic_model/matlab/model_implementation

# run fitting code (should have one job per substrate)
srun matlab -nodisplay -r "init_path(), run_off_target_fit_cluster(${SEED}, ${SLURM_ARRAY_TASK_ID}, ${RANK}), quit()"
