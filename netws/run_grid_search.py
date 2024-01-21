import os
import signal
import subprocess
from itertools import product
from multiprocessing import Pool, cpu_count, Manager

seed = 42
cpus = 8

finished_runs = 0
total_runs = 0

def print_progress(shared):
    loading_bar_len = 50

    bar_progress = (shared['finished'] * loading_bar_len) // shared['total']

    print(f'\r[{"=" * bar_progress}{" " * (loading_bar_len - bar_progress)}] ({shared["finished"]}/{shared["total"]})    ', end='')

def run_program(params, shared):
    global seed

    command = [
        'python',
        './network.py',
        '--tag', 'protrusion_v1',
        '--hidden-layers', *map(str, params['hidden_layers']),
        '--seed', str(seed),
        '--ligand', params['ligand'],
        '--learning-rate', str(params['learning_rate']),
        '--epochs', str(params['epochs']),
        '--batch-size', str(params['batch_size']),
        '--embedder', str(params['embedder']),
        '--epoch-stats-interval', str(params['epoch_stats_interval']),
        '--protrusion-data-file', r'..\data\3d_proc\protrusions.big.json',
        '--pdb-mappings-fname', r'..\data\3d_proc\mappings_to_pdbs.json',
    ]
    
    subprocess.run(command)
    
    shared['finished'] += 1
    print_progress(shared)

parameters_to_test = {
    'hidden_layers': [ [900, 500, 100, 30], [100, 50, 20], [400, 200, 50]],
    #'ligand': ['ADP', 'AMP', 'ATP', 'CA', 'DNA', 'FE', 'GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN'],
    'ligand': ['GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN', 'AMP'], # just some of the ligands
    'learning_rate': [0.01, 0.001],
    'epochs': [20],
    'batch_size': [1000],
    'epoch_stats_interval': [10],
    'embedder': ['BERT'],
    # 'embedder': ['ESM', 'T5', 'BERT'],
}

param_combinations = [
    dict(zip(parameters_to_test.keys(), values)) 
    for values in product(*parameters_to_test.values())
]

def run_parallel():
    global cpus

    manager = Manager()
    shared = manager.dict({'finished': 0, 'total': len(param_combinations)})

    print('Starting the grid search')
    print(f'Running on {cpus} CPUs')
    print_progress(shared)

    with Pool(processes=cpus) as pool:
        pool.starmap(run_program, [(param, shared) for param in param_combinations])
    
    print('\n\nGrid search finished')

if __name__ == "__main__":
    run_parallel()
