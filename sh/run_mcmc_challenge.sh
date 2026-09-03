#!/bin/bash
#SBATCH --job-name=MCMC_DC2
#SBATCH --output=./projects/cpip_data_challenge_2/logs/%x_%a_%A.out
#SBATCH --error=./projects/cpip_data_challenge_2/logs/%x_%a_%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --ntasks-per-socket=2
#SBATCH --cpus-per-task=8
#SBATCH --time=4-00:00:00
#SBATCH --partition=standard
#SBATCH --account=timeifler

echo Job starting at `date` on node `hostname`

YAML=./projects/cpip_data_challenge_2/yamls/MCMC_${SLURM_ARRAY_TASK_ID}.yaml

# Clear the environment from any previously loaded modules
# module purge > /dev/null 2>&1
module load micromamba
source ~/.bashrc

cd $SLURM_SUBMIT_DIR
micromamba activate cocoa
sleep $(( SLURM_ARRAY_TASK_ID*3 ))
source start_cocoa.sh

export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mpirun -n ${SLURM_NTASKS} --oversubscribe --mca pml ob1 \
  --mca btl vader,tcp,self --bind-to core:overload-allowed \
  --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} \
  --mca mpi_yield_when_idle 1 \
  --report-bindings \
  cobaya-run ${YAML} -r
