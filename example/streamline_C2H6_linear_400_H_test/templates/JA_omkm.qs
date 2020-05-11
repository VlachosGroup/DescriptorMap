#!/bin/bash -l
#SBATCH --array=1-__n_jobs%__n_concurrent
#SBATCH --job-name="__job_name omkm run"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=__partition
#SBATCH --time=__run_time
#SBATCH --output="__log_pathomkm-%A_%a.out"
#SBATCH --export=NONE
#SBATCH --no-step-tmpdir
#SBATCH --mail-user=__mail-user
#SBATCH --mail-type=ALL

# Set up environment
export VALET_PATH=/work/ccei_biomass/sw/valet
vpkg_require openmkm

# Read folder file
cd "./omkm"
FOLDER_FILE='./folderlist.txt'
#Change to the job directory
FOLDER=$(sed -n "$SLURM_ARRAY_TASK_ID p" "$FOLDER_FILE")
cd "$FOLDER"
echo Running $FOLDER
# Run OMKM
omkm reactor.yaml thermo.xml