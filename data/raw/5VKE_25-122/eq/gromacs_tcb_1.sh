#!/bin/bash

# Submit to the tcb partition
#SBATCH -p tcb

# The name of the job in the queue
#SBATCH -J 5VKE_25-122_1
# wall-clock time given to this job
#SBATCH -t 23:30:00

# Number of nodes and number of MPI processes per node
#SBATCH -N 1 -n 8 
# Request a GPU node and two GPUs (per node)
# Remove if you don't need a GPU
#SBATCH -C gpu --gres=gpu:2 -x gpu04

# Output file names for stdout and stderr
#SBATCH -e job-%j.err -o job-%j.out

# Receive e-mails when your job starts and ends
#SBATCH --mail-user=sergio.perez.conesa@scilifelab.se --mail-type=FAIL

###YOU PROBABLY WANT TO EDIT STUFF BELOW HERE
module unload gromacs
module load gromacs/2020.5

var=1
time=23


var0=$((($var-1)))
if [ ! -f "step${var}.cpt" ]; then
	gmx grompp -maxwarn 1 -f ../../../../models/mdp/NPTres${var}.mdp  -c step${var0}.gro  -p ../topol.top  -o topol${var}.tpr -n ../index.ndx -pp topol_pp.top -r step${var0}.gro
    cmd="gmx mdrun -v -maxh $time -s topol${var}.tpr  -pin on -deffnm step$var "
	echo $cmd
	$cmd
	if [ ! -f "step${var}.gro" ]; then
		sbatch gromacs_tcb_${var}.sh
	else
		var=$((($var+1)))
		sbatch gromacs_tcb_${var}.sh
        fi 
else
    cmd="gmx mdrun -v -maxh $time -s topol${var}.tpr  -pin on -deffnm step$var -cpi step${var}.cpt"
    echo $cmd
    $cmd
	if [ ! -f "step${var}.gro" ]; then
		sbatch gromacs_tcb_${var}.sh
	else
		var=$((($var+1)))
		sbatch gromacs_tcb_${var}.sh
        fi 
fi
