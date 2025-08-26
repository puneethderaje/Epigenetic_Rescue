#!/bin/bash
#SBATCH --time=07:00:00
#SBATCH --account=def-mmosmond
#SBATCH --cpus-per-task=80
#SBATCH --nodes=5

module load CCEnv
module load StdEnv
module load nixpkgs/16.09
module load gcc/7.3.0
module load parallel/20160722


source ~/.virtualenvs/myenv/bin/activate

python ER_ParamutableMutant.py 




