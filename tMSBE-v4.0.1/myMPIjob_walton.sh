######################################

#FOR WHOLE MOTHERBOARD USE
#nonodes=1  # Number of MPI processes [NOTE: (numnodes = nodes/2) for nodes found above]
#numthds=44 # number of cores assigned to each MPI process

# FOR WHOLE CPU (2 per motherboard) USE
#nonodes=7 # Number of MPI processes [NOTE: numnodes = nodes, for nodes found above]
#numthds=6 # Number of OpenMP threads per MPI job
nonodes=1
numthds=3

# Calculate used threads for Hybrid OpenMP and MPI
let totThreads="$nonodes*$numthds"
let numCoresPerNode="$numthds"
let totCores="$numCoresPerNode*$nonodes"

# Disable MPI processor affinity and allow for select placement
#export MV2_ENABLE_AFFINITY=0

#With Intel compiler versions 10.1.015 and later,
#you need to set KMP_AFFINITY to disabled
#to avoid the interference between dplace and
#Intel's thread affinity interface.
#export KMP_AFFINITY=disabled

export MPI_DSM_VERBOSE=1
export OMP_NUM_THREADS=$totThreads
#export OMP_PROC_BIND=true
#export OMP_PLACES=sockets
#export GOMP_CPU_AFFINITY="0-10"
#export GOMP_CPU_AFFINITY="11-21"
#export GOMP_CPU_AFFINITY="22-32"
#export GOMP_CPU_AFFINITY="33-43"

# Disable MPI processor affinity and allow for select placement
export MV2_ENABLE_AFFINITY=0

#With Intel compiler versions 10.1.015 and later,
#you need to set KMP_AFFINITY to disabled
#to avoid the interference between dplace and
#Intel's thread affinity interface.
export KMP_AFFINITY=disabled

#numactl --hardware
#mpirun -np $nonodes omplace -nt $numthds -vv ./a.out > output.txt
#mpirun -np $nonodes ./a.out > output.txt

# Run on cpus on node 0, with preferred memory allocation on local node(s)
#mpirun -np $nonodes numactl --cpunodebind=1 --localalloc ./a.out > output.txt
#mpirun -np $nonodes --cpus-per-proc $numthds --report-bindings ./a.out > output.txt
#mpirun -np $nonodes taskset --cpu-list 0-10 ./a.out > output.txt
#mpirun -np $nonodes taskset --cpu-list 11-21 ./a.out > output.txt
#mpirun -np $nonodes taskset --cpu-list 22-32 ./a.out > output.txt
#mpirun -np $nonodes taskset --cpu-list 33-43 ./a.out > output.txt
#mpirun -np $nonodes numactl --localalloc taskset --cpu-list $1 ./a.out > output.txt
mpirun -np $nonodes numactl --localalloc taskset --cpu-list '41-43' ./a.out > output.txt
