nonodes=1 # Number of MPI processes [NOTE: numnodes = nodes, for nodes found above]
numthds=2  # Number of cores assigned to each MPI process

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

########################################
