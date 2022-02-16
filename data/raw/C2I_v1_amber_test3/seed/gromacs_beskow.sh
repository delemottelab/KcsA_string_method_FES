#!/bin/bash -l

# Include your allocation number
#SBATCH -A 2020-3-28
# Name your job
#SBATCH -J swarm
# Total number of nodes and MPI tasks
#SBATCH --nodes 4 --ntasks-per-node=32
# length in hours
#SBATCH -t  00:30:00

export OMP_NUM_THREADS=1
# APRUN_OPTIONS="-n 64 -d 1 -cc none"
module swap PrgEnv-cray PrgEnv-gnu
module add gromacs/2020.5


var=0
echo "New Sim $var">>beskow_queue

time=0.45


if [ ! -f "state.cpt" ]; then
	gmx grompp -maxwarn 1 -f swarms.mdp  -o topol.tpr -n ./index.ndx -pp topol_pp.top
    cmd="srun  gmx_mpi mdrun -v -maxh $time -s topol.tpr  -pin on -cpi state.cpt "
	echo $cmd
	$cmd
	if [ ! -f "confout.gro" ]; then
		sbatch gromacs_beskow.sh
        fi
else
    cmd="srun  gmx_mpi mdrun -v -maxh $time -s topol.tpr  -pin on  -cpi state.cpt"
    echo $cmd
    $cmd
	if [ ! -f "confout.gro" ]; then
		sbatch gromacs_beskow.sh
        fi
fi
