import json
import os
import argparse
import os
import random
import signal
import subprocess
from itertools import product
from multiprocessing import Pool, cpu_count, Manager

seed = 42

# python /home/brabecm4/diplomka/protein-binding-sites/run_scripts/create_run_cases.py --out-file /home/brabecm4/diplomka/protein-binding-sites/data/run_data/run_cases_output

parser = argparse.ArgumentParser(description='Description of your script')
parser.add_argument('--out-file', help='Output File')

args = parser.parse_args()

out_file = args.out_file

def get_command(params):
    global seed

    command = [
        'python',
        '/home/brabecm4/diplomka/protein-binding-sites/netws/network.v2.py',
        # '/home/brabecm4/diplomka/protein-binding-sites/netws/network.composed.py',
        '--tag', 'protrusion_no_emb',
        '--verbose', 'True',
        # '--use-simple-model', 'True',
        '--hidden-layers', *map(str, params['hidden_layers']),
        '--seed', str(seed),
        '--ligand', params['ligand'],
        '--learning-rate', str(params['learning_rate']),
        '--epochs', str(params['epochs']),
        '--batch-size', str(params['batch_size']),
        '--embedder', str(params['embedder']),
        '--epoch-stats-interval', str(params['epoch_stats_interval']),

        # '--radius', str(params['radius']),
        
        # '--neighbors', str(params['neighbors']),
        # '--epochs-sec', str(params['epochs_sec']),
        # '--learning-rate-sec', str(params['learning_rate_sec']),
        # '--hidden-layers-sec', *map(str, params['hidden_layers_sec']),
    ]

    return ' '.join(command)

parameters_to_test = {
    # 'hidden_layers': [ [900, 500, 100, 30], [100, 50, 20], [400, 200, 50], [90, 20]],
    'hidden_layers': [
        # [1024, 512, 512, 256],
        # [1024, 512, 256, 32],
        # [1024, 512, 512, 256, 64],
        # [1024, 512, 256, 128, 128, 64],
        # [1024, 512, 256, 128],
        # [1024, 512, 128, 32],
        # [2048, 512, 512, 256, 128],
        # [2048, 512, 256, 128],
        # [2048, 256, 256, 64],
        # [2048, 256, 64, 32],
        # [2048, 256, 32],
        # [4096, 256, 256, 64],
        # [1024, 256, 64],
        # [4096, 2048, 1024, 512, 256, 128, 64, 32],
        # [4096, 2048, 512, 256, 32],
        # [1024, 512, 128, 32],
        # [2048, 256, 128, 32],
        # [2048, 256, 128],
        # [1024, 256, 128],
        # [2048, 1024, 512, 256],
        # [2048, 1024, 512, 256],
        # [8192, 2048, 1024, 512, 256],
        # [512, 512, 512, 512],
        # [256, 256, 128, 32],
        # [512, 512, 512, 32],
        # [512, 512, 32],
        # [512, 256],
        # [256, 64, 32],
        # [256, 256, 128, 128],
        # [256, 256, 32],
        # [512, 512, 128, 128],
        # [256, 32, 32],
        # [256, 128, 32, 32],
        # [512, 128, 128, 128],
        # [512, 128, 128, 32],

        # [256, 32, 4],
        # [128, 32, 4],
        [64, 32, 2],
        [16, 2],
        [8, 2],
        # [256, 64, 2],
        # [512, 64, 4],

        
        # these are the best of bests
        # [256, 128],
        # [512, 512, 32],
        # [256, 128, 128],
        # [256, 256, 256],
        # [256, 128, 128, 32],
        
#   [512, 512, 32, 32],
#   [256, 256, 256, 128],
#   [512, 256, 128, 128],
#   [256, 128],
#   [256, 128, 128, 32],
#   [512, 256, 128, 32],
#   [512, 256, 256, 128],
#   [256, 256, 256],
#   [256, 256, 128],
#   [128, 32],
#   [512, 512, 32],
#   [256, 128, 128],
#   [128, 128, 32],
#   [512, 512],
#   [128, 32, 32, 32],
#   [256, 32],
#   [512, 256, 256, 32],
#   [512, 512, 256, 128],
#   [512, 512, 256, 256],
#   [512, 256, 256],
#   [512, 32, 32, 32],
#   [128, 128, 128, 32],
#   [512],
#   [512, 512, 512, 128],
#   [128, 128, 32, 32],
#   [256, 256, 32, 32],
#   [512, 256, 256, 256],
#   [512, 256, 32, 32],
#   [256],
#   [512, 512, 256],
#   [512, 512, 128, 32],
#   [512, 32, 32],
#   [128, 128, 128],
#   [512, 512, 128],
#   [128],
#   [256, 256, 256, 32],
#   [512, 128, 32, 32],
#   [256, 128, 128, 128],
#   [256, 256, 256, 256],
#   [128, 128],
#   [512, 128],
#   [512, 512, 512, 256],
#   [512, 512, 512],
#   [128, 32, 32],
#   [256, 128, 32],
#   [512, 256, 32],
#   [256, 256],
#   [512, 512, 512, 512],
#   [256, 256, 128, 32],
#   [512, 512, 512, 32],
#   [256, 32, 32, 32],
#   [512, 256, 128],
#   [512, 512, 256, 32],
#   [512, 256],
#   [512, 128, 128],
#   [512, 128, 32],
#   [128, 128, 128, 128],
#   [512, 32],
#   [256, 256, 128, 128],
#   [256, 256, 32],
#   [512, 512, 128, 128],
#   [256, 32, 32],
#   [256, 128, 32, 32],
#   [512, 128, 128, 128],
#   [512, 128, 128, 32],
], 
    'ligand': ['ADP', 'AMP', 'ATP', 'CA', 'DNA', 'FE', 'GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN'],
    # 'ligand': ['GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN' ], # just some of the ligands
    'learning_rate': [0.01, 0.001 , 0.0001],
    # 'epochs': [20, 30],
    'epochs': [120],
    'batch_size': [1000],
    'epoch_stats_interval': [10],
    'embedder': ['ESM'],
    # 'embedder': ['ESM', 'T5', 'BERT'],

    # for protrusion
    # 'radius': [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],

    # for composed network
    # 'neighbors': [6, 10],
    # 'epochs_sec': [120],
    # 'learning_rate_sec': [0.01, 0.0001],
    # 'hidden_layers_sec': [
    #     [256, 128],
    #     [512, 512, 32],
    #     [256, 256, 256],
    #     [256, 128, 128, 32],
    # ],
}

param_combinations = [
    dict(zip(parameters_to_test.keys(), values)) 
    for values in product(*parameters_to_test.values())
]

commands = [get_command(param) for param in param_combinations]


lines = [f'NOT_RUN;{command}' for command in commands]
random.shuffle(lines)

content = '\n'.join(lines)

with open(out_file, 'w') as file:
    file.write(content)

print(f'Created {len(commands)} cases')
