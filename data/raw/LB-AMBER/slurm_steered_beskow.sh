#!/bin/bash -l

# Include your account in your cluster.
#SBATCH --account=2021-3-15
# The name of the job in the queue
#SBATCH --job-name=steered_v2
#SBATCH -C Haswell

# Output file names for stdout and stderr
#SBATCH --error=slurm_out/steer.err
#SBATCH --output=slurm_out/steer.out

# Add your email below.
# Receive e-mails when your job fails
#SBATCH --mail-user=sergiopc@kth.se
#SBATCH --mail-type=FAIL

# Time should slightly greater than the time for one full iteration
# if you want to do one iteration per slurm job (recomended).
# If you can allocate big chunks of time in your cluster you can put the time
# of N-iterations and add this number of interations in the variable
# `max_iteration=$((($iteration+1)))` bellow.
#SBATCH --time=0:05:00

# Total number of nodes and MPI tasks
# This number of nodes and tasks has been found to work well for 60-80k atoms in beskow (@DelemotteLab).
# You can of course adapt it to your HPC environment following the guidelines of the main README.md

# Number of nodes and number of MPI tasks per node
#SBATCH --nodes=4
# In slurm jargon tasks is like MPI-ranks
#SBATCH --ntasks-per-node=32

# Choose version of gromacs
module unload gromacs
module load gromacs/2020.5

# Path to string-method-gmxapi
path_string_method_gmxapi=../../../../string-method-gmxapi/

# Path to anaconda3 instalation with string_method environment
my_conda_env=/cfs/klemming/nobackup/s/sergiopc/anaconda3
env >environment


######################  DON'T MODIFY ###############################
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
var=$("${my_conda_env}/bin/conda" 'shell.bash' 'hook' 2> /dev/null)
__conda_setup=$var
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "$my_conda_env/etc/profile.d/conda.sh" ]; then
        . "$my_conda_env/etc/profile.d/conda.sh"
    else
        export PATH="$my_conda_env/bin/$PATH"
    fi
fi
unset __conda_setup
conda activate string_method

cmd=" `which python` ${path_string_method_gmxapi}/stringmethod/main.py --config_file=config_steered.json --start_mode=steered"
echo $cmd
$cmd
