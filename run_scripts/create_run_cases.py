import json
import os
import argparse
import os
import signal
import subprocess
from itertools import product
from multiprocessing import Pool, cpu_count, Manager

seed = 42

parser = argparse.ArgumentParser(description='Description of your script')
parser.add_argument('--out-file', help='Output File')

args = parser.parse_args()

out_file = args.out_file

def get_command(params):
    global seed

    command = [
        'python',
        '/home/brabecm4/diplomka/protein-binding-sites/netws/network.py',
        '--tag', 'basic_v2',
        '--verbose', 'True',
        '--use-simple-model', 'True',
        '--hidden-layers', *map(str, params['hidden_layers']),
        '--seed', str(seed),
        '--ligand', params['ligand'],
        '--learning-rate', str(params['learning_rate']),
        '--epochs', str(params['epochs']),
        '--batch-size', str(params['batch_size']),
        '--embedder', str(params['embedder']),
        '--epoch-stats-interval', str(params['epoch_stats_interval']),
        # if not present protrusion is not used
        # '--protrusion-data-file', r'/home/brabecm4/diplomka/protein-binding-sites/data/3d_proc/protrusions.big.json',
        '--pdb-mappings-fname', r'/home/brabecm4/diplomka/protein-binding-sites/data/3d_proc/mappings_to_pdbs.json',
    ]

    return ' '.join(command)

parameters_to_test = {
    # 'hidden_layers': [ [900, 500, 100, 30], [100, 50, 20], [400, 200, 50], [90, 20]],
    'hidden_layers': [
  [512, 512, 32, 32],
  [256, 256, 256, 128],
  [512, 256, 128, 128],
  [256, 128],
  [256, 128, 128, 32],
  [512, 256, 128, 32],
  [512, 256, 256, 128],
  [256, 256, 256],
  [256, 256, 128],
  [128, 32],
  [512, 512, 32],
  [256, 128, 128],
  [128, 128, 32],
  [512, 512],
  [128, 32, 32, 32],
  [256, 32],
  [512, 256, 256, 32],
  [512, 512, 256, 128],
  [512, 512, 256, 256],
  [512, 256, 256],
  [512, 32, 32, 32],
  [128, 128, 128, 32],
  [512],
  [512, 512, 512, 128],
  [128, 128, 32, 32],
  [256, 256, 32, 32],
  [512, 256, 256, 256],
  [512, 256, 32, 32],
  [256],
  [512, 512, 256],
  [512, 512, 128, 32],
  [512, 32, 32],
  [128, 128, 128],
  [512, 512, 128],
  [128],
  [256, 256, 256, 32],
  [512, 128, 32, 32],
  [256, 128, 128, 128],
  [256, 256, 256, 256],
  [128, 128],
  [512, 128],
  [512, 512, 512, 256],
  [512, 512, 512],
  [128, 32, 32],
  [256, 128, 32],
  [512, 256, 32],
  [256, 256],
  [512, 512, 512, 512],
  [256, 256, 128, 32],
  [512, 512, 512, 32],
  [256, 32, 32, 32],
  [512, 256, 128],
  [512, 512, 256, 32],
  [512, 256],
  [512, 128, 128],
  [512, 128, 32],
  [128, 128, 128, 128],
  [512, 32],
  [256, 256, 128, 128],
  [256, 256, 32],
  [512, 512, 128, 128],
  [256, 32, 32],
  [256, 128, 32, 32],
  [512, 128, 128, 128],
  [512, 128, 128, 32],
], 
    
    'ligand': ['ADP', 'AMP', 'ATP', 'CA', 'DNA', 'FE', 'GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN'],
    # 'ligand': ['GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN' ], # just some of the ligands
    'learning_rate': [0.01, 0.001, 0.0001],
    'epochs': [100],
    'batch_size': [1000],
    'epoch_stats_interval': [10],
    'embedder': ['T5', 'BERT', 'ESM'],
    # 'embedder': ['ESM', 'T5', 'BERT'],
}

param_combinations = [
    dict(zip(parameters_to_test.keys(), values)) 
    for values in product(*parameters_to_test.values())
]

commands = [get_command(param) for param in param_combinations]


lines = [f'NOT_RUN;{command}' for command in commands]
content = '\n'.join(lines)

with open(out_file, 'w') as file:
    file.write(content)

print('DONE')
