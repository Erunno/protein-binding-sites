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
        '/home/brabecm4/diplomka/protein-binding-sites/netws/compare_estimators.py',
        '--ligand', params['ligand'],
    ]

    return ' '.join(command)

parameters_to_test = {
    'ligand': ['ADP', 'AMP', 'ATP', 'CA', 'DNA', 'FE', 'GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN'],
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
