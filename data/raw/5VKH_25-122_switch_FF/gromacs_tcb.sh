#!/bin/bash

# Submit to the tcb partition
#SBATCH -p tcb

# The name of the job in the queue
#SBATCH -J switch_FF
# wall-clock time given to this job
#SBATCH -t 23:30:00

# Number of nodes and number of MPI processes per node
#SBATCH -N 1
# Request a GPU node and two GPUs (per node)
# Remove if you don't need a GPU
#SBATCH -C gpu --gres=gpu:2

# Output file names for stdout and stderr
#SBATCH -e job-%j.err -o job-%j.out

# Receive e-mails when your job starts and ends
#SBATCH --mail-user=sergio.perez.conesa@scilifelab.se --mail-type=FAIL

#module swap PrgEnv-cray PrgEnv-gnu


###YOU PROBABLY WANT TO EDIT STUFF BELOW HERE
module unload gromacs
module load gromacs/2020.1
time=23


cmd="gmx mdrun -nt 24 -v -maxh $time -s topol.tpr  -pin on -cpi state.cpt"
echo $cmd
$cmd
err=$?
if [   $err == 0 ]; then
if [ ! -f "confout.gro" ]; then
	sbatch gromacs_tcb.sh
fi
fi
