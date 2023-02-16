#!/bin/bash -l
#SBATCH --array=1-__n_jobs%__n_concurrent
#SBATCH --job-name="__job_name mkm run"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=__partition
#SBATCH --time=__run_time
#SBATCH --output="./log/pathmkm-%A_%a.out"
#SBATCH --export=NONE
#SBATCH --no-step-tmpdir
# SBATCH --mail-user=__mail-user
#SBATCH --mail-type=ALL

# Set up environment
export VALET_PATH=/work/ccei_biomass/sw/valet

# Do standard OpenMP environment setup:
. /opt/shared/slurm/templates/libexec/openmp.sh

vpkg_require chemkin-reactor/2019.0.120:intel
vpkg_require graphviz

# Read folder file
cd "./omkm"
FOLDER_FILE='./folderlist.txt'
#Change to the job directory
FOLDER=$(sed -n "$SLURM_ARRAY_TASK_ID p" "$FOLDER_FILE")
cd "$FOLDER"
echo Running $FOLDER

# Do standard CHEMKIN reactor job environment setup:
. /work/ccei_biomass/sw/chemkin-reactor/jobenv-setup.sh

.  graphviz.sh

.  refineVisualization.sh

# Each CHEMKIN program should be executed via srun to ensure proper resource
# setup and tracking.

#
# Prep work:
#
CHEMKIN_PREP
rc=$?
if [ $rc -eq 0 ]; then
	#
	# Execute the model:
	#
	CHEMKIN_EXEC
	rc=$?
fi
exit $rc