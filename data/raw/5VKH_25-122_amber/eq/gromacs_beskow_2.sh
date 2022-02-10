#!/bin/bash -l

# Include your allocation number
#SBATCH -A 2020-3-28 
# Name your job
#SBATCH -J AMBER_2
# Total number of nodes and MPI tasks
#SBATCH --nodes 4 --ntasks-per-node=32 
# length in hours
#SBATCH -t  00:30:00 

export OMP_NUM_THREADS=1
APRUN_OPTIONS="-n 64 -d 1 -cc none"
module swap PrgEnv-cray PrgEnv-gnu
module add gromacs/2020.5


var=2
echo "New Sim $var">>beskow_queue

date >>beskow_queue
time=0.45


var0=$((($var-1)))
if [ ! -f "step${var}.cpt" ]; then
	gmx grompp -maxwarn 1 -f NPTres${var}.mdp  -c step${var0}.gro  -p ../topol.top  -o topol${var}.tpr -n ../index.ndx -pp topol_pp.top -r step${var0}.gro
    cmd="srun  gmx_mpi mdrun -v -maxh $time -s topol${var}.tpr  -pin on -deffnm step$var "
	echo $cmd
	$cmd
	if [ ! -f "step${var}.gro" ]; then
		sbatch gromacs_beskow_${var}.sh
	else
		var=$((($var+1)))
		sbatch gromacs_beskow_${var}.sh
        fi 
else
    cmd="srun  gmx_mpi mdrun -v -maxh $time -s topol${var}.tpr  -pin on -deffnm step$var -cpi step${var}.cpt"
    echo $cmd
    $cmd
	if [ ! -f "step${var}.gro" ]; then
		sbatch gromacs_beskow_${var}.sh
	else
		var=$((($var+1)))
		sbatch gromacs_beskow_${var}.sh
        fi 
fi
date >>beskow_queue
