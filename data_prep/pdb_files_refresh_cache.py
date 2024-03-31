import argparse
import json
import os
import re
import time

import numpy as np
from file_cache import use_cache
import pdb_files_db
import traceback
import datasets_db
from datetime import timedelta

# python3 /home/brabecm4/diplomka/protein-binding-sites/data_prep/pdb_files_refresh_cache.py

# constant parameters

RESET = '\033[0m'
RED = '\033[91m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
SPACE = '                       '

# init database

parser = argparse.ArgumentParser()
parser.add_argument('--check-presence', action='store_true', help='Check if all methods are cached', required=False)

parser.add_argument('--remove-old', action='store_true', help='Remove old cache', required=False)
parser.add_argument('--running-from-worker', action='store_true', help='Append new lines to progress bar.', required=False)
parser.add_argument('--prots', type=str, nargs='+', help='List of proteins to refresh')

parser.add_argument('--start', type=int, help='Start protein index')
parser.add_argument('--count', type=int, help='Count of proteins')

args = parser.parse_args()

protrusion_radii = list(np.arange(1.0, 10.5, 0.5))
protrusion_function_names = pdb_files_db.ProtrusionFunctionsCollection.names_to_functions.keys()

_from, _to = 0, 100000
end_line = ''

if args.running_from_worker:
    end_line = '\n'
    RESET, RED, GREEN, YELLOW, BLUE, SPACE = '', '', '', '', '', ''

if args.start is not None:
    _from = args.start
    _to = args.start + args.count
    print ('From: ', _from, ', To: ', _to)

cache_storage = '/home/brabecm4/diplomka/protein-binding-sites/data/cache_data/production'

print(f'{YELLOW}Initializing...{RESET}')
pdb_db = pdb_files_db.PdbFilesDb()
chains_db = datasets_db.SeqDatasetDb()

if args.remove_old:
    print(f'{YELLOW}Removing old cache...{RESET}')
    for file_name in os.listdir(cache_storage):
        os.remove(os.path.join(cache_storage, file_name))

print(f'{YELLOW}Loading data...{RESET}')
all_chains = chains_db.get_all_chain_records()
all_chains = sorted(all_chains, key=lambda x: x.full_id())[_from:_to]

if args.prots is not None: 
    all_chains = [ch for ch in all_chains if ch.full_id() in args.prots]

start_time = time.time()
def estimate_to_end(done, total):
    current_time = time.time() - start_time
    progress = done / total

    if progress == 0:
        return '???'

    estimated_total_time = current_time / progress
    remaining_time = int(estimated_total_time - current_time)

    seconds = remaining_time % 60
    minutes = (remaining_time // 60) % 60
    hours = remaining_time // (60 * 60)

    return f'{hours:2}h {minutes:2}m {seconds:2}s'

print(f'{YELLOW}Running methods...{RESET}')

def run_functions_to_cache(chain_structure: pdb_files_db.Chain3dStructure):

    #################
    # residue count #
    #################

    residue_count = chain_structure.get_residue_count()

    ####################
    # nearest residues #
    ####################

    for res_i in range(residue_count):
        nearest = chain_structure.get_nearest_residue_indexes(res_i)

    ############
    # sequence #
    ############
        
    sequence = chain_structure.compute_sequence()

    ##############
    # protrusion #
    ##############

    for radius in protrusion_radii:

        # the version with default algorithm
        protrusion_vector = chain_structure.get_protrusion_vector(
            radius=radius
        )

        for algorithm in protrusion_function_names:
            protrusion_vector = chain_structure.get_protrusion_vector(
                radius=radius,
                protrusion_algorithm=algorithm
            )

    ###############
    # SASA vector #
    ###############
            
    SASA_vector = chain_structure.get_SASA_vector()

for i, chain in enumerate(all_chains):
    print(f'\r{BLUE}[{chain.full_id():5}] computing chain {YELLOW}{i + 1:4}{BLUE} / {YELLOW}{len(all_chains)}{BLUE} ... remaining: {estimate_to_end(i, len(all_chains))}' + \
            f' {RESET}', end=end_line, flush=True)
    
    chain_structure = pdb_db.get_chain_structure(chain.protein_id(), chain.chain_id())

    if not args.check_presence:
        # for checking the presence of data in the cache all that is needed it to
        # not load 3d file of the chain. If the method call misses, the program will fail
        # as the chain below does not have data to compute with
        chain_structure.load(preferred_sequence=chain.sequence())

    with use_cache(chain_structure, data_folder=cache_storage) as cached_chain_structure:
        run_functions_to_cache(cached_chain_structure)

    chain_structure.free_memory()


print(f'\n{YELLOW}DONE{RESET}')