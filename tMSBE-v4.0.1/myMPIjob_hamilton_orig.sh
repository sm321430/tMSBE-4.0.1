######################################
#!/bin/sh
#
#PBS -l select=1:ncpus=3:mpiprocs=3:ompthreads=1:mem=300GB
#PBS -l walltime=900:00:00
#PBS -j oe

# nonodes = select*mpiprocs
nonodes=3 # Number of MPI processes [NOTE: numnodes = nodes, for nodes found above]
numthds=1  # Number of cores assigned to each MPI process



# -- should not have to edit below this line ------------------------------------------------
# Switch to the working directory; by default PBS launches processes from your home directory.
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

### Output node information
echo Running on host `hostname`
echo Start Time is `date`
echo Directory is `pwd`
echo This jobs runs on the following processors:
echo `cat $PBS_NODEFILE`
# Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`

export MPI_DSM_VERBOSE=1
export OMP_NUM_THREADS=$numthds

#With Intel compiler versions 10.1.015 and later,
#you need to set KMP_AFFINITY to disabled
#to avoid the interference between dplace and
#Intel's thread affinity interface.
export KMP_AFFINITY=disabled

ulimit -S -s unlimited

#omplace -nt 10 -tm intel -vv ./a.out > output.txt
#mpirun -np $nonodes omplace -nt $numthds -vv ./a.out > output.txt
mpiexec_mpt -np $nonodes omplace -nt $numthds -c 672-674 -vv ./a.out > output.txt

sleep 10
qstat -n $PBS_JOB_ID > placement_info.out
echo End Time is `date`
########################################
#ps -C a.out -L -opsr,comm,time,pid,ppid,lwp > placement.out
