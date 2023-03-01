# Load GROMACS
export PATH="/work/donglab/software/gromacs-2021.5/build/bin/:$PATH"
module load gcc/10.1.0 cuda/11.2 cmake/3.18.1 openmpi/4.1.0-gcc10.1-cuda11.2

# Load AMBER
module load gcc/7.3.0 
module load intel/mkl-2020u2 
module load plumed/2.6.1-skylake 
module load oracle_java/jdk1.8.0_181 
module load amber/20-mpi 
source /shared/centos7/amber/amber20/skylake-mpi/amber.sh 
source /shared/centos7/amber/amber20/skylake-mpi/miniconda/bin/activate 
