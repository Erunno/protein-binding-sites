#!/bin/bash

#SBATCH -p gpu-short
#SBATCH --gres=gpu:V100
#SBATCH -A nprg058s
#SBATCH --cpus-per-task=4
#SBATCH --time=0-02:00:00
#SBATCH --signal=B:USR1@30
#SBATCH --output=/home/brabecm4/diplomka/protein-binding-sites/data/logs/%A_%a.out

if [ -z "$1" ]; then
  echo "Usage: $0 <filename with command>"
  exit 1
fi

source /home/brabecm4/diplomka/protein-binding-sites/python_envs/pytorch_envs/bin/activate

echo "launching worker..." 

srun /home/brabecm4/diplomka/protein-binding-sites/run_scripts/run_commands.sh "$1"

echo "worker has ended"
