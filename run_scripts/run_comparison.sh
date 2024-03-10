#!/bin/bash

#SBATCH -p gpu-long
#SBATCH --gres=gpu:V100
#SBATCH -A nprg058s
#SBATCH --cpus-per-task=2
#SBATCH --time=0-20:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --signal=B:USR1@30
#SBATCH --output=/home/brabecm4/diplomka/protein-binding-sites/data/logs/%A_%a.out

if [ -z "$1" ] || [ -z "$2" ]; then
  echo "Usage: $0 <ligand> <radius>"
  exit 1
fi

echo "Starting the job"
start_time=$(date +%s.%N)

source /home/brabecm4/diplomka/protein-binding-sites/python_envs/pytorch_envs/bin/activate
python3 /home/brabecm4/diplomka/protein-binding-sites/netws/compare_estimators.py --ligand $1 --radius $2

end_time=$(date +%s.%N)
elapsed_time=$(echo "$end_time - $start_time" | bc)

echo "Job has finished in: $elapsed_time seconds"