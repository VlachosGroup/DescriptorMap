#!/bin/bash -l
#SBATCH --job-name="DescMap analyze"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=ccei_biomass
#SBATCH --time=24:00:00
#SBATCH --output="./log/analysis-%j.out"
#SBATCH --export=NONE
#SBATCH --no-step-tmpdir
# SBATCH --mail-user=xzong@udel.edu
#SBATCH --mail-type=ALL

# Set up environment
export VALET_PATH=/work/ccei_biomass/sw/valet
vpkg_require anaconda/5.2.0:python3
source activate descmap

# Run scripts
cd "./analysis"
python geometric_analysis.py
