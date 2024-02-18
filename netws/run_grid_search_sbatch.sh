#!/bin/bash

#SBATCH -p gpu-short
#SBATCH --gres=gpu:V100
#SBATCH -A nprg058s
#SBATCH --cpus-per-task=64
#SBATCH --time 02:00:00
#SBATCH --output=/home/brabecm4/diplomka/protein-binding-sites/data/logs/grid_search.log

source /home/brabecm4/diplomka/protein-binding-sites/python_envs/pytorch_envs/bin/activate

python3 run_grid_search.py