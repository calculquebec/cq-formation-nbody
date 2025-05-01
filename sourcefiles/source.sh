if [ "$1" == "2020" ]; then
  module load StdEnv/2020
  module load gcc/11.3 openmpi/4.1.4
  module load python/3.11 scipy-stack/2023a mpi4py/3.1.4
else
  module load StdEnv/2023
  module load gcc/12.3 openmpi/4.1.5
  module load python/3.11 scipy-stack/2025a mpi4py/4.0
fi
