#!/bin/bash

#SBATCH -p gpu-short
#SBATCH -A nprg058s
#SBATCH --cpus-per-task=2
#SBATCH --time=0-02:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --signal=B:USR1@30
#SBATCH --output=/home/brabecm4/diplomka/protein-binding-sites/data/logs/pdb_cache_renew_%A_%a.out

if [ -z "$1" ] || [ -z "$2" ]; then
  echo "Usage: $0 <start protein idx> <protein count>"
  exit 1
fi

source /home/brabecm4/diplomka/protein-binding-sites/python_envs/pytorch_envs/bin/activate

echo "launching worker..." 

python3 /home/brabecm4/diplomka/protein-binding-sites/data_prep/pdb_files_refresh_cache.py --running-from-worker --start $1 --count $2

echo "worker has ended"
