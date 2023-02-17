#!/bin/bash -l
#SBATCH --job-name="__job_name setup"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=__partition
#SBATCH --time=__setup_time
#SBATCH --output="./log/setup-%j.out"
#SBATCH --export=NONE
#SBATCH --no-step-tmpdir
#SBATCH --mail-user=__mail-user
#SBATCH --mail-type=ALL

# Set up environment
export VALET_PATH=/work/ccei_biomass/sw/valet
vpkg_require anaconda/5.2.0:python3
source activate descmap_dev

# Run scripts
cd "./setup"
python setup.py