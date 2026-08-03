#!/bin/bash
#SBATCH --job-name=MCMC_BPCA_VAL
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

# Clear the environment from any previously loaded modules
module purge > /dev/null 2>&1
source ~/.bashrc 

cd $SLURM_SUBMIT_DIR
conda activate cocoa
source start_cocoa

export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mpirun -n ${SLURM_NTASKS} --oversubscribe --mca pml ^ucx \
  --mca btl vader,tcp,self --bind-to core:overload-allowed \
  --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} \
  --mca mpi_yield_when_idle 1 \
  cobaya-run ./projects/cpip_data_challenge_2/yamls/BPCA_VAL_${SLURM_ARRAY_TASK_ID}.yaml -r
