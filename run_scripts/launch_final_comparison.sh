#!/bin/bash

#SBATCH -p gpu-long
#SBATCH --gres=gpu:V100
#SBATCH -A nprg058s
#SBATCH --cpus-per-task=4
#SBATCH --time=0-20:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --signal=B:USR1@30
#SBATCH --output=/home/brabecm4/diplomka/protein-binding-sites/data/logs/comparison_%A_%a.out

if [ -z "$1" ] || [ -z "$2" ]; then
  echo "Usage: $0 <tag> <lignd>"
  exit 1
fi

source /home/brabecm4/diplomka/protein-binding-sites/python_envs/pytorch_envs/bin/activate

echo "launching worker..." 

# networks with compressing layer
python3 /home/brabecm4/diplomka/protein-binding-sites/netws/network.final_eval.composed.py --final-tag $1 --ligand $2

# double protrusion
# python3 /home/brabecm4/diplomka/protein-binding-sites/netws/network.final_eval.py --final-tag $1 --ligand $2 --second-radius $3

# normal networks
# python3 /home/brabecm4/diplomka/protein-binding-sites/netws/network.final_eval.py --final-tag $1 --ligand $2


echo "worker has ended"
