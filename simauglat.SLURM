#!/bin/bash

# To submit this file, use: $sbatch simauglat.SLURM

#------------------------------------------------
# Configuration of this SLURM batch job

#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fweber144@protonmail.com
#SBATCH --job-name=auglat
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

#SBATCH --nodes=3
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=16

####SBA|TCH --mem-per-cpu=8000
#SBATCH --mem=24000

#SBATCH --time=15:00:00

#------------------------------------------------
# Configuration of modules

module load R/4.2.0 openmpi/gcc/4.1.1 cmake/3.20.4

#------------------------------------------------
# Configuration of MPI

export OMPI_MCA_mtl=^psm

#------------------------------------------------
# Write job info to output file

echo "------------------------------------------------"
echo "User ID: $USER"
echo "The ID of the job allocation: $SLURM_JOB_ID"
echo "Total number of nodes in the job's resource allocation: $SLURM_JOB_NUM_NODES"
echo "List of nodes allocated to the job: $SLURM_JOB_NODELIST"
echo "Job was started on nodes: $SLURM_NODELIST"
echo "ID of the nodes allocated: $SLURM_NODEID"
echo "Number of tasks (same as --ntasks): $SLURM_NTASKS"
echo "Number of tasks requested per node (only set if the --ntasks-per-node option is specified): $SLURM_NTASKS_PER_NODE"
echo "Number of tasks to be initiated on each node (values are comma separated and in the same order as SLURM_JOB_NODELIST [...]): $SLURM_TASKS_PER_NODE"
echo "Number of CPUs on the allocated node (important for the number of threads (OMP, OPENBLAS, MKL)): $SLURM_CPUS_ON_NODE"
echo "Number of CPUs requested per task (only set if the --cpus-per-task option is specified): $SLURM_CPUS_PER_TASK"
echo "Same as --mem-per-cpu: $SLURM_MEM_PER_CPU"
echo "Same as --mem: $SLURM_MEM_PER_NODE"
echo "Submit directory (the directory from which sbatch was invoked or, if applicable, the directory specified by the -D, --chdir option): $SLURM_SUBMIT_DIR"
echo "Home directory: $HOME"
### Not defined (but not harmful neither):
# echo "Scratch directory: $SCRATCH"
###
echo "------------------------------------------------"

#------------------------------------------------
# Job content

mpirun --map-by socket:NOOVERSUBSCRIBE -np 1 R --no-save < sim.R

#------------------------------------------------
# Exit

exit 0
