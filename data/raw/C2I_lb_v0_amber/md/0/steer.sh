#!/bin/bash -l

# Include your allocation number
#SBATCH -A 2020-3-28
# Name your job
#SBATCH -J steerXXXX
# Total number of nodes and MPI tasks
#SBATCH --nodes 4 --ntasks-per-node=32
# length in hours
#SBATCH -t  00:30:00

export OMP_NUM_THREADS=2
module swap PrgEnv-cray PrgEnv-gnu
module add gromacs/2020.5


echo "New Sim $var">>beskow_queue

date >>beskow_queue
time=0.45


var=XXXX
var0=$((($var-1)))
if [ ! -f "state.cpt" ]; then
	gmx grompp -maxwarn 1 -f grompp.mdp  -c ../../$var0/restrained/confout.gro  -p ../../../../topology/topol.top  -o topol.tpr -n ../../../../topology/index.ndx -r ../../$var0/restrained/confout.gro
    cmd="srun  gmx_mpi mdrun -v -maxh $time -s topol.tpr  -pin on"
	echo $cmd
	$cmd
	if [ ! -f "confout.gro" ]; then
		sbatch steer.sh
	else
		var=$((($var+1)))
        cd ../../$var/restrained/
		sbatch steer.sh
        fi
else
    cmd="srun  gmx_mpi mdrun -v -maxh $time -s topol.tpr  -pin on -cpi state.cpt"
    echo $cmd
    $cmd
	if [ ! -f "confout.gro" ]; then
		sbatch steer.sh
	else
		var=$((($var+1)))
        cd ../../$var/restrained/
		sbatch steer.sh
        fi
fi
date >>beskow_queue
