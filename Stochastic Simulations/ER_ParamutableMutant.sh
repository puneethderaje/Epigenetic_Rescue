#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --account=def-mmosmond
#SBATCH --cpus-per-task=80
#SBATCH --nodes=5

module load CCEnv
module load StdEnv
module load nixpkgs/16.09
module load gcc/7.3.0
module load parallel/20160722


source ~/.virtualenvs/myenv/bin/activate

srun="srun --exclusive -N1 -n1 "

parallel="parallel --delay 0.2 -j 5 --joblog runtask_ER_Diploid.log --resume"

echo "First Run Started" 

$parallel $srun python ER_ParamutableMutant.py {1} ::: 0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99




