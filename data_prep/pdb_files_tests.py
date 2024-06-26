import argparse
import json
import os
import sys
import time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import config.config as config

from file_cache import use_cache
import pdb_files_db
import traceback
import datasets_db
from datetime import timedelta

##################################
# set path to cache store folder #
##################################
path_to_test_cache = '/home/brabecm4/diplomka/protein-binding-sites/data/cache_data/for_tests'

# python3 /home/brabecm4/diplomka/protein-binding-sites/data_prep/pdb_files_tests.py

parser = argparse.ArgumentParser(description='Test pdb db')
parser.add_argument('--tests', type=str, nargs='+', help='List of tests to be run')
args = parser.parse_args()

# constant parameters

structure_storage_folder = config.pdbs_folder
path_to_embeddings_folder = config.embeddings_folder
path_to_sequences_folder = config.yu_sequences_folder

origin_file_key = 'structure_file'
sequence_key = 'sequence_from_3d_structure'

sequence_errors_threshold = 0.02
unmarked_alpha_carbon_count_threshold = 6

RESET = '\033[0m'
RED = '\033[91m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
SPACE = '                       '

# init database

print(f'{YELLOW}Initializing...')
pdb_db = pdb_files_db.PdbFilesDb(storage_folder=structure_storage_folder)
chains_db = datasets_db.SeqDatasetDb(
    sequences_folder=path_to_sequences_folder, 
    embeddings_folder=path_to_embeddings_folder)

def sequence_error_does_not_exceed_threshold():
    all_chains = chains_db.get_all_chain_records()
    
    wrong_count = 0

    for i, chain in enumerate(all_chains):
        print(f'\r{BLUE}[ ] {sequence_error_does_not_exceed_threshold.__name__} ... testing chain {YELLOW}{i}{BLUE} - wrong chains: {RED}{wrong_count}{BLUE} {RESET}', end='', flush=True)
        reference_sequence = chain.sequence()
        
        chain_structure = pdb_db.get_chain_structure(chain.protein_id(), chain.chain_id())
        chain_structure.load(preferred_sequence=reference_sequence)

        sequence = chain_structure.compute_sequence()

        chain_structure.free_memory()

        if reference_sequence != sequence:
            wrong_count += 1
    
    total_count = len(all_chains)
    err_percent = wrong_count / total_count

    if err_percent > sequence_errors_threshold:
        return False, f'{wrong_count} sequences differs -- error rate {err_percent:.4f} is grater than allowed threshold {sequence_errors_threshold}'

    return True, f'({RED}{wrong_count}{GREEN} wrong out of {YELLOW}{total_count}{GREEN}) ... {wrong_count / total_count * 100:.2f} % wrong'

def nearest_residues_does_not_raise_an_error():
    all_chains = chains_db.get_all_chain_records()
    start_time = time.time()

    for i, chain in enumerate(all_chains):
        print(f'\r{BLUE}[ ] {nearest_residues_does_not_raise_an_error.__name__} ... testing chain {YELLOW}{i}{BLUE} {RESET}', end='', flush=True)
        
        chain_structure = pdb_db.get_chain_structure(chain.protein_id(), chain.chain_id())
        chain_structure.load(preferred_sequence=chain.sequence())

        residue_count = chain_structure.get_residue_count()
        
        for res_i in range(residue_count):
            n_nearest = 5
            nearest = chain_structure.get_nearest_residue_indexes(
                residue_index=res_i)[:n_nearest]

            if len(nearest) != n_nearest:
                return False, f'Unexpected length of nearest residue list - expected {n_nearest}, got {len(nearest)}'            

        chain_structure.free_memory()

    end_time = time.time()
    elapsed_time_secs = end_time - start_time

    return True, f'elapsed time: {timedelta(seconds=elapsed_time_secs)}'

def nearest_residues_are_cashable():
    # remove all data in cache
    for file_name in os.listdir(path_to_test_cache):
        os.remove(os.path.join(path_to_test_cache, file_name))

    all_chains = chains_db.get_all_chain_records()[:10]

    times = [0, 0, 0]
    
    def to_str(time):
        return f'{timedelta(seconds=time)}'

    for i, chain in enumerate(all_chains):
        print(f'\r{BLUE}[ ] {nearest_residues_are_cashable.__name__} ... testing chain {YELLOW}{i}{BLUE} ... ' + \
              f'not cached: {to_str(times[0])}, cached: {to_str(times[1])}, cold cache: {to_str(times[2])} {RESET}', end='', flush=True)
        
        chain_structure = pdb_db.get_chain_structure(chain.protein_id(), chain.chain_id())
        chain_structure.load(preferred_sequence=chain.sequence())

        residue_count = chain_structure.get_residue_count()

        def run_nearest_residues_test(run_index, tested_chain_structure):
            start_time = time.time()

            for res_i in range(residue_count):
                n_nearest = 5

                nearest = tested_chain_structure.get_nearest_residue_indexes(res_i)[:n_nearest]

                if len(nearest) != n_nearest:
                    return False, f'Unexpected length of nearest residue list - expected {n_nearest}, got {len(nearest)}'            

            end_time = time.time()
            elapsed_time_secs = end_time - start_time
            
            times[run_index] += elapsed_time_secs

        
        with use_cache(chain_structure, data_folder=path_to_test_cache) as cached_chain_structure:
            for run_index in range(2):
                run_nearest_residues_test(run_index, cached_chain_structure)

        # cold cache test
        with use_cache(chain_structure, data_folder=path_to_test_cache) as cached_chain_structure:
            run_nearest_residues_test(2, cached_chain_structure)
                
        chain_structure.free_memory()


    return True, f'elapsed time: not cached: {to_str(times[0])}, cached: {to_str(times[1])}, cold cache: {to_str(times[2])}'

def check_residues_without_marked_alpha_carbon():
    all_chains = chains_db.get_all_chain_records()

    chains_with_unmarked_alpha_carbon_count = 0
    valid_chains_with_unmarked_count = 0

    for i, chain in enumerate(all_chains):
        print(f'\r{BLUE}[ ] {check_residues_without_marked_alpha_carbon.__name__} ... testing chain {YELLOW}{i}{BLUE}, ' + \
              f'unmarked {YELLOW}{chains_with_unmarked_alpha_carbon_count}{BLUE} of which {YELLOW}{valid_chains_with_unmarked_count}{BLUE} are valid chains{RESET}', end='', flush=True)
        
        chain_structure = pdb_db.get_chain_structure(chain.protein_id(), chain.chain_id())
        chain_structure.load(preferred_sequence=chain.sequence())
        
        if not chain_structure.all_residues_has_alpha_carbon():
            chains_with_unmarked_alpha_carbon_count += 1

            if chain_structure.compute_sequence() == chain.sequence():
                valid_chains_with_unmarked_count += 1

        chain_structure.free_memory()

    if chains_with_unmarked_alpha_carbon_count > unmarked_alpha_carbon_count_threshold:
        return False, f'Unmarked alpha carbon count has exceeded the threshold of {unmarked_alpha_carbon_count_threshold} as it was {chains_with_unmarked_alpha_carbon_count}'

    if valid_chains_with_unmarked_count > 0:
        return False, f'Unmarked alpha carbon found in valid chain ... total unmarked {chains_with_unmarked_alpha_carbon_count} of which {valid_chains_with_unmarked_count} are valid chains'

    return True, f'{chains_with_unmarked_alpha_carbon_count} chains has unmarked CA of which {valid_chains_with_unmarked_count} are valid chains (...which is within the threshold of {unmarked_alpha_carbon_count_threshold})'

def check_protrusion_does_not_throw_an_error():
    all_chains = chains_db.get_all_chain_records()
    start_time = time.time()

    for i, chain in enumerate(all_chains):
        print(f'\r{BLUE}[ ] {nearest_residues_does_not_raise_an_error.__name__} ... testing chain {YELLOW}{i}{BLUE} {RESET}', end='', flush=True)
        
        chain_structure = pdb_db.get_chain_structure(chain.protein_id(), chain.chain_id())
        chain_structure.load(preferred_sequence=chain.sequence())

        residue_count = chain_structure.get_residue_count()

        protrusion = chain_structure.get_protrusion_vector(radius=10)

        if protrusion is None or residue_count != len(protrusion):
            return False, f'Unexpected protrusion vector of length {len(protrusion)}, expected {residue_count}'

        chain_structure.free_memory()

    end_time = time.time()
    elapsed_time_secs = end_time - start_time

    return True, f'elapsed time: {timedelta(seconds=elapsed_time_secs)}'    

def check_SASA_vector_does_not_throw_an_error():
    all_chains = chains_db.get_all_chain_records()
    start_time = time.time()

    for i, chain in enumerate(all_chains):
        print(f'\r{BLUE}[ ] {nearest_residues_does_not_raise_an_error.__name__} ... testing chain {YELLOW}{i}{BLUE} {RESET}', end='', flush=True)
        
        chain_structure = pdb_db.get_chain_structure(chain.protein_id(), chain.chain_id())
        chain_structure.load(preferred_sequence=chain.sequence())

        residue_count = chain_structure.get_residue_count()

        SASA_vector = chain_structure.get_SASA_vector()

        if SASA_vector is None or residue_count != len(SASA_vector):
            return False, f'Unexpected SASA vector of length {len(SASA_vector)}, expected {residue_count}'

        chain_structure.free_memory()

    end_time = time.time()
    elapsed_time_secs = end_time - start_time

    return True, f'elapsed time: {timedelta(seconds=elapsed_time_secs)}'    


# run test cases

def _run(test_case):
    print(f'{BLUE}  ...running {test_case.__name__}{RESET}', end='', flush=True)
    try:
        success, message = test_case()
    except Exception as e:
        success = False
        message = f'Exception occurred: "{e}"\n\n{traceback.format_exc()}'

    if success:
        print(f'\r{GREEN}[✔] {test_case.__name__} {"- " + message if message is not None else ""}{RESET}{SPACE}')
    else:
        print(f'\r{RED}[✘] {test_case.__name__}{RESET}{SPACE}')
        print(f'  Error: {message}')

def _skip(test_case):
    print(f'{YELLOW}[ ] {test_case.__name__} {RESET}', flush=True)

print(f'\n{BLUE}Running tests...{RESET}\n')

all_symbols = globals()
tests = [
    func 
    for func in all_symbols.values() 
    if callable(func) and 
       func.__module__ == __name__ and
       not func.__name__.startswith('_')]

test_to_be_run = args.tests or [t.__name__ for t in tests] 

for test_case in tests:
    if test_case.__name__ in test_to_be_run:
        _run(test_case)
    else:
        _skip(test_case)

print()