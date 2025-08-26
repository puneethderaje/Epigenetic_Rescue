#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH --account=def-mmosmond
#SBATCH --cpus-per-task=80
#SBATCH --nodes=1


module load StdEnv
module load nixpkgs/16.09
module load gcc/7.3.0
module load parallel/20160722


source ~/.virtualenvs/newenv/bin/activate

srun="srun --exclusive -N1 -n1 "

parallel="parallel --delay 0.2 -j 5"

echo "First Run Started" 

$parallel $srun python ER_Haploid_PureEpi.py $1 {1} ::: 0.01 0.2 0.5 0.8 1.0




