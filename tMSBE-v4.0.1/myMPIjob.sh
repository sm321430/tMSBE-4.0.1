######################################
#!/bin/sh
#
#PBS -l nodes=1:ppn=20
#PBS -l walltime=900:00:00
#PBS -j oe

#FOR WHOLE MOTHERBOARD USE
nonodes=10 # Number of MPI processes [NOTE: (numnodes = nodes/2) for nodes found above]
numthds=2 # number of OpenMP threads to use for each MPI node

# Calculate used threads for Hybrid OpenMP and MPI
let totThreads="$nonodes*$numthds"
let numCoresPerNode="$numthds"
let totCores="$numCoresPerNode*$nonodes"

echo "Number of nodes should be = " $totCores


# -- should not have to edit below this line ------------------------------------------------
# Switch to the working directory; by default PBS launches processes from your home directory.
echo "Working directory is " $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

### Output node information
echo "Running on host" $(hostname)
echo "Start Time is " $(date)
echo "Directory is " $(pwd)
echo "This jobs runs on the following processors: "
echo $(cat $PBS_NODEFILE)
# Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`

export MPI_DSM_VERBOSE=1
export OMP_NUM_THREADS=$totThreads

# Prevents OpenMP threads from spawning at random cores/sockets
export OMP_PROC_BIND=true

# Disable MPI processor affinity and allow for select placement
export MV2_ENABLE_AFFINITY=0

#With Intel compiler versions 10.1.015 and later,
#you need to set KMP_AFFINITY to disabled
#to avoid the interference between dplace and
#Intel's thread affinity interface.
export KMP_AFFINITY=disabled

ulimit -S -s unlimited

# BOOLE
#omplace -nt 10 -tm intel -vv ./a.out > output.txt
#mpirun -np $nonodes omplace -nt $numthds -vv ./a.out > output.txt

# WALTON/STOKES/BOYLE
if [[ "$numthds" -eq 1 ]]; then
# OPTION 1: SPACE NODES OVER CORE/SOCKET
mpirun -np $nonodes -map-by core ./a.out > output.txt
#mpirun -np $nonodes -map-by socket ./a.out > output.txt
else
# OPTION 2: RUN ALL NODES IN ORDER
mpirun -np $nonodes -map-by slot:pe=$numCoresPerNode ./a.out > output.txt
fi

sleep 10
qstat -n $PBS_JOB_ID > placement_info.out
echo "End Time is " $(date)
########################################
#ps -C a.out -L -opsr,comm,time,pid,ppid,lwp > placement.out
