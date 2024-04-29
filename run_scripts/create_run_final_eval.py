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

# python /home/brabecm4/diplomka/protein-binding-sites/run_scripts/create_run_final_eval.py --out-file /home/brabecm4/diplomka/protein-binding-sites/data/run_data/run_cases_output

parser = argparse.ArgumentParser(description='Description of your script')
parser.add_argument('--out-file', help='Output File')

args = parser.parse_args()

out_file = args.out_file

def get_command(params):
    global seed

    command = [
        'python',
        '/home/brabecm4/diplomka/protein-binding-sites/netws/network.final_eval.py',
        '--final-tag', params['tag'],
        '--ligand', params['ligand'],
    ]

    return ' '.join(command)

parameters_to_test = {
    # first batch:

    # 'tag': [
    #     'basic_v6',
    #     'prot_all_c',
    #     'protrusion_bypass_v4_c',
    #     'SASA_bypassed_v2_c',
    #     'neighboring_emb_5_v2_c',
    # ],
    
    # second batch:
    
    # 'tag': [
    #     'one_prot_fst_v3_c',
    #     'nei_emb_3_v2_c',
    #     'nei_emb_5_avrg_v2_c',
    #     'nei_emb_3_avrg_v2_c',
    # ],
    
    # third batch:
    
    'tag': [
        'nei_3_comprs_v2_c',
        'nei_5_comprs_v2_c',
        'SASA_fst_v2_c',
    ],

    'ligand': ['ADP', 'AMP', 'ATP', 'CA', 'DNA', 'FE', 'GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN'],
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
