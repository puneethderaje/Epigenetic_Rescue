#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --account=def-mmosmond
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1

module load CCEnv
module load StdEnv
module load nixpkgs/16.09
module load gcc/7.3.0
module load parallel/20160722

srun="srun --exclusive -N1 -n1 --cpus-per-task 40"

source ~/.virtualenvs/myenv/bin/activate

$srun python Stoch_PreChange_Diploid.py 

