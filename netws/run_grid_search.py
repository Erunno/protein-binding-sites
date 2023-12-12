import subprocess
from itertools import product
from multiprocessing import Pool, cpu_count

seed = 42

def run_program(params):
    global seed

    command = [
        'python',
        './network.py',
        '--verbose', 'False',
        '--hidden-layers', *map(str, params['hidden_layers']),
        '--seed', str(seed),
        '--ligand', params['ligand'],
        '--learning-rate', str(params['learning_rate']),
        '--epochs', str(params['epochs']),
        '--batch-size', str(params['batch_size']),
        '--embedder', str(params['embedder']),
        '--epoch-stats-interval', str(params['epoch_stats_interval'])
    ]
    
    subprocess.run(command)


parameters_to_test = {
    'hidden_layers': [[1000, 800, 200, 100], [400, 200, 50], [50, 10, 5]],
    # 'ligand': ['ADP'], #, 'AMP', 'ATP', 'CA', 'DNA', 'FE', 'GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN'],
    'ligand': ['ADP', 'AMP', 'ATP', 'CA', 'DNA', 'FE', 'GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN'],
    'learning_rate': [0.01, 0.001],
    'epochs': [50],
    'batch_size': [1000],
    'epoch_stats_interval': [10],
    'embedder': ['BERT']
}

param_combinations = [
    dict(zip(parameters_to_test.keys(), values)) 
    for values in product(*parameters_to_test.values())
]

def run_parallel():
    with Pool(processes=cpu_count()) as pool:
        pool.map(run_program, param_combinations)

if __name__ == "__main__":
    run_parallel()