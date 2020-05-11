#!/bin/bash -l
#SBATCH --job-name="__job_name analyze"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=__partition
#SBATCH --time=__analyze_time
#SBATCH --output="__log_pathanalysis-%j.out"
#SBATCH --export=NONE
#SBATCH --no-step-tmpdir
#SBATCH --mail-user=__mail-user
#SBATCH --mail-type=ALL

# Set up environment
export VALET_PATH=/work/ccei_biomass/sw/valet
vpkg_require anaconda/5.2.0:python3
source activate descmap_dev

# Run scripts
cd "./analysis"
python analysis.py