#!/bin/bash -l

# Include your account in your cluster.
#SBATCH --account=snic2021-3-15
# The name of the job in the queue
#SBATCH --job-name=LB-AMBER_s1
#SBATCH --partition main

# Output file names for stdout and stderr
#SBATCH --error=slurm.out
#SBATCH --output=slurm.err

# Add your email below.
# Receive e-mails when your job fails
#SBATCH --mail-user=sergiopc@kth.se
#SBATCH --mail-type=FAIL

# Time should slightly greater than the time for one full iteration
# if you want to do one iteration per slurm job (recomended).
# If you can allocate big chunks of time in your cluster you can put the time
# of N-iterations and add this number of interations in the variable
# `max_iteration=$((($iteration+1)))` bellow.
#SBATCH --time=0:30:00

# Total number of nodes and MPI tasks
# This number of nodes and tasks has been found to work well for 60-80k atoms in beskow (@DelemotteLab).
# You can of course adapt it to your HPC environment following the guidelines of the main README.md

# Number of nodes and number of MPI tasks per node
#SBATCH --nodes=4
# In slurm jargon tasks is like MPI-ranks
#SBATCH --ntasks-per-node=128

# Choose version of gromacs
ml PDC
ml GROMACS/2020.5-cpeCray-21.11
# ml Anaconda3/2021.05


time=0.45
cmd="srun  gmx_mpi mdrun -v -maxh $time -s topol.tpr  -pin on  -cpi state.cpt"
echo $cmd
$cmd
err=$?
if [   $err == 0 ]; then
if [ ! -f "confout.gro" ]; then
	sbatch gromacs_dardel.sh
fi
fi
