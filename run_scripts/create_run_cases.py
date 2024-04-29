import argparse
import random
from itertools import product

import numpy as np

seed = 42

# python /home/brabecm4/diplomka/protein-binding-sites/run_scripts/create_run_cases.py --out-file /home/brabecm4/diplomka/protein-binding-sites/data/run_data/run_cases_output

parser = argparse.ArgumentParser(description='Create run cases')
parser.add_argument('--out-file', help='Output File')

args = parser.parse_args()

out_file = args.out_file

def get_command(params):
    global seed

    command = [
        'python',
        # '/home/brabecm4/diplomka/protein-binding-sites/netws/network.v2.py',
        '/home/brabecm4/diplomka/protein-binding-sites/netws/network.composed.v2.py',
        '--tag', f'nei_compr_{params["neighbors"]}_small_v4_c',
        '--neighbors', str(params["neighbors"]),
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
    'neighbors': [5, 3],
    'hidden_layers': [
        [1024, 512, 512, 256],
        # ...
    ], 
    'ligand': ['ADP', 'AMP', 'ATP', 'CA', 'DNA', 'FE', 'GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN'],
    'learning_rate': [0.01, 0.001 , 0.0001],
    'epochs': [120],
    'batch_size': [1000],
    'epoch_stats_interval': [10],
    'embedder': ['ESM'],
    # 'radius': list(np.arange(1.0, 10.5, 0.5)),
    # 'embedder': ['ESM', 'T5', 'BERT'],
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
