#!/bin/bash --login

# PBS job options (name, compute nodes, job time)
#PBS -N Scaling_MPI_Job
#PBS -l select=8
#PBS -l walltime=00:20:00

# Replace [budget code] below with your project code (e.g. t01)
#PBS -A e370

# Make sure any symbolic links are resolved to absolute path
export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)               

# Change to the directory that the job was submitted from
# (remember this should be on the /work filesystem)
cd $PBS_O_WORKDIR

# Set the number of threads to 1
#   This prevents any system libraries from automatically 
#   using threading.
export OMP_NUM_THREADS=1

# Launch the parallel job
#   Using 1536 MPI processes and 24 MPI processes per node
aprun -n 192 -N 24 ./felixrefine