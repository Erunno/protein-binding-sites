import re

import os
import sys
from typing import List
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))

import numpy as np
import datasets_db
import pdb_files_db
import traceback

# constant parameters
path_to_protrusion_file = '/home/brabecm4/diplomka/protein-binding-sites/data/3d_proc/protrusion.max-neighbors.big.json'
path_to_embeddings_folder = '/home/brabecm4/diplomka/protein-binding-sites/data/embedded_sequences'
path_to_sequences_folder = '/home/brabecm4/diplomka/protein-binding-sites/data/orig/yu_sequences'
tested_embedder = 'ESM'
all_protrusion_radii = list(np.arange(1.0, 10.5, 0.5))

# constants
ALL_LIGANS = ['ADP', 'AMP', 'ATP', 'CA', 'DNA', 'FE', 'GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN']

RESET = '\033[0m'
RED = '\033[91m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
SPACE = '                            '

# init database

print(f'{YELLOW}Initializing...')
db = datasets_db.SeqDatasetDb(
    sequences_folder=path_to_sequences_folder, 
    embeddings_folder=path_to_embeddings_folder)

pdb_db = pdb_files_db.PdbFilesDb()
db.set_pdb_db(pdb_db) 

print(f'Loading protrusion...{RESET}')
db.load_protrusion_data_file(path_to_protrusion_file)

def can_load_all_datasets():
    all_chains = db.get_all_chain_records()

    if all_chains is None or len(all_chains) == 0:
        return False, 'No chains found'

    return True, None
    
def can_access_protein_id(): 
    all_chains = db.get_all_chain_records()

    for chain in all_chains:
        prot_id = chain.protein_id()

        if prot_id is None or prot_id.strip() == '':
            return False, f'Protein id was "{prot_id}"'
    
    return True, None

def can_access_chain_id(): 
    all_chains = db.get_all_chain_records()

    for chain in all_chains:
        ch_id = chain.chain_id()

        if ch_id is None or ch_id.strip() == '':
            return False, f'Chain id was "{ch_id}"'

    return True, None

def can_access_full_id(): 
    all_chains = db.get_all_chain_records()

    for chain in all_chains:
        full_id = chain.full_id()

        if full_id is None or full_id.strip() == '':
            return False, f'Chain id was "{full_id}"'

    return True, None

def protein_id_has_correct_format(): 
    all_chains = db.get_all_chain_records()

    for chain in all_chains:
        prot_id = chain.protein_id()

        if not re.match('^[a-z0-9]{4}$', prot_id):
            return False, f'Protein id had incorrect format: "{prot_id}"'

    return True, None


def chain_id_has_correct_format(): 
    all_chains = db.get_all_chain_records()

    for chain in all_chains:
        ch_id = chain.chain_id()

        if not re.match('^[A-Z0-9]$', ch_id):
            return False, f'Chain id had incorrect format: "{ch_id}"'
        
    return True, None

def full_id_has_correct_format(): 
    all_chains = db.get_all_chain_records()

    for chain in all_chains:
        full_id = chain.full_id()

        if not re.match('^[a-z0-9]{4}[A-Z0-9]$', full_id):
            return False, f'Full id had incorrect format: "{full_id}"'
        
    return True, None

def binding_sights_are_merged_correctly_from_all_records(): 
    for ligand in datasets_db.all_ligands:
        ds = db.get_dataset_for(ligand)

        chains_per_binding_sight = ds.testing_per_binding_sight() + ds.training_per_binding_sight()

        aggregated_chains = {}
        for chain in chains_per_binding_sight:
            if chain.full_id() in aggregated_chains:
                aggregated_chains[chain.full_id()].append(chain)
            else:
                aggregated_chains[chain.full_id()] = [chain]

        all_chains = ds.testing() + ds.training()
        
        for chain in all_chains:
            aggr_chains: List[datasets_db.ChainRecord] = aggregated_chains[chain.full_id()]
            binding_sights = np.array([0] * len(chain.sequence()))

            for ch in aggr_chains:
                binding_sights = binding_sights | ch.original_binding_sights()

            binding_sights_from_ds = chain.binding_sights()

            if np.any(binding_sights != binding_sights_from_ds):
                return False, f'Unexpected binding sights -- lig: {ligand}, prot: {chain.full_id()}'
    
    return True, None

def there_are_some_chains_with_protrusion_records():
    all_chains = db.get_all_chain_records()
    all_chains = datasets_db.Helpers.filter_chains_with_protrusion(all_chains)

    if len(all_chains) == 0:
        return False, 'No chains with protrusion record'

    return True, None 

def protrusion_radii_are_correct(): 
    all_chains = db.get_all_chain_records()
    all_chains = datasets_db.Helpers.filter_chains_with_protrusion(all_chains)

    for chain in all_chains:
        radii = all_protrusion_radii

        if radii is None or len(radii) == 0:
            return False, f'Unexpected radii list: "{radii}"'      

    return True, None

def all_ligand_datasets_not_empty():
    for ligand in ALL_LIGANS:
        ds = db.get_dataset_for(ligand)

        testing = ds.testing()
        training = ds.training()

        if testing is None or len(testing) == 0:
            return False, f'Unexpected testing data for ligand {ligand}, got "{testing}"'

        if training is None or len(training) == 0:
            return False, f'Unexpected testing data for ligand {ligand}, got "{training}"'

    return True, None

def protrusion_does_not_throw_exception(): 
    all_chains = db.get_all_chain_records()
    all_chains = datasets_db.Helpers.filter_chains_with_protrusion(all_chains)

    for i, chain in enumerate(all_chains):
        print(f'\r {BLUE} {protrusion_does_not_throw_exception.__name__} ... {YELLOW}{i}{BLUE} checked out of {YELLOW}{len(all_chains)}{BLUE}{RESET}', end='', flush=True)

        radii = all_protrusion_radii

        for r in radii: 
            protrusion = chain.protrusion_vector_for(r)

            if protrusion is None or len(protrusion) == 0:
                return False, f'Unexpected protrusion values: "{protrusion}"'        

    return True, None

# def all_radii_function_on_db_returns_correct_value():
#     ligand_ds = db.get_dataset_for('AMP')
#     all_chains = datasets_db.Helpers.filter_chains_with_protrusion(ligand_ds.all()) 

#     expected_all_radii = sorted(all_chains[0].get_all_protrusion_radii())
#     actual_all_radii = ligand_ds.get_all_radii()

#     are_same = all(x == y for x, y in zip(expected_all_radii, actual_all_radii))

#     if not are_same:
#         return False, f'Get all radii failed - expected {expected_all_radii}, got {actual_all_radii}'

#     return True, None

def most_of_the_records_matches_protrusion_sequence():
    max_invalid_percent = 0.02

    all_chains = db.get_all_chain_records()
    
    valid, invalid = datasets_db.Helpers.split_chains_to_valid_and_invalid_protrusion(all_chains)

    if len(valid) + len(invalid) != len(all_chains):
        return False, 'valid + invalid != all'

    invalid_percent = len(invalid) / len(all_chains)

    if invalid_percent > max_invalid_percent:
        return False, f'Max number of invalid chains exceeded - was {invalid_percent * 100} %'
    
    return True, None

def data_accessor_for_binding_sights_works():
    all_chains = db.get_all_chain_records()
    binding_sights_accessor = datasets_db.DataAccessors.biding_sights_vect()

    for chain in all_chains:
        expected = chain.binding_sights()
        actual = binding_sights_accessor(chain)

        if not (expected == actual).all():
            return False, f'Binding sights accessor returned different value'
        
    return True, None

def data_accessor_for_protrusion_works():
    all_chains = db.get_all_chain_records()
    all_chains = datasets_db.Helpers.filter_chains_with_protrusion(all_chains)

    binding_sights_accessor = datasets_db.DataAccessors.protrusion(1.0, 2.0, 4.0)

    for chain in all_chains:
        prot1 = chain.protrusion_vector_for(1.0)
        prot2 = chain.protrusion_vector_for(2.0)
        prot4 = chain.protrusion_vector_for(4.0)

        expected = np.array([[prot1[i], prot2[i], prot4[i]] for i in range(len(prot1))])
        actual = binding_sights_accessor(chain)

        if not (expected == actual).all():
            return False, f'Protrusion sight accessor returned different value'

    return True, None

def concat_chain_data_works():
    the_chain_1 = 1 
    the_chain_2 = 10 

    data = [the_chain_1, the_chain_2]

    def mock_accessor_1(chain):
        return np.array([1, 2, 3]) * chain
    
    def mock_accessor_2(chain):
        return np.array([[11, 12], [21, 22], [31, 32]]) * chain
    
    def mock_accessor_3(chain):
        return np.array([100, 200, 300]) * chain

    actual = datasets_db.Helpers.concat_chain_data(
        mock_accessor_1, mock_accessor_2, mock_accessor_3,
        chains=data
    )

    expected = np.array([
        [1, 11, 12, 100],
        [2, 21, 22, 200],
        [3, 31, 32, 300],
        [10, 110, 120, 1000],
        [20, 210, 220, 2000],
        [30, 310, 320, 3000],
    ])

    if not (actual == expected).all():
        return False, f'Unexpected return value from {datasets_db.Helpers.concat_chain_data.__name__}'

    return True, None    

def can_load_all_embeddings():
    all_chains = db.get_all_chain_records()

    for i, chain in enumerate(all_chains):
        embeddings = chain.embeddings(tested_embedder)
        expected_len = len(chain.sequence())

        print(f'\r  {BLUE} {can_load_all_embeddings.__name__} {YELLOW}{i}{BLUE} checked out of {YELLOW}{len(all_chains)}{BLUE}{RESET}', end='', flush=True)

        if embeddings is None or len(embeddings) != expected_len:
            return False, f'Unexpected embeddings: len={len(embeddings)} should be {expected_len} \n"{embeddings}"'
        
        chain.free_embeddings_from_RAM()

    return True, None

def data_accessor_for_embeddings_works():
    all_chains = db.get_all_chain_records()
    
    embedding_accessor = datasets_db.DataAccessors.embeddings(tested_embedder)

    i = 0
    for chain in all_chains:
        expected_embeddings = chain.embeddings(tested_embedder)
        actual = embedding_accessor(chain)

        i += 1
        if i % 20 == 0:
            print(f'\r  {BLUE}checking embeddings - {YELLOW}{i}{BLUE} checked out of {YELLOW}{len(all_chains)}{BLUE}{RESET}', end='', flush=True)

        if not (expected_embeddings == actual).all():
            return False, f'Embeddings accessor returned different value'
        
        chain.free_embeddings_from_RAM()

    return True, None

def embeddings_format_is_correct():
    all_chains = db.get_all_chain_records()
    expected_embedding_vector_size = len(all_chains[0].embeddings(tested_embedder)[0]) 

    i = 0
    for chain in all_chains:
        embeddings = chain.embeddings(tested_embedder)
        bindings_flags = chain.binding_sights()

        i += 1
        if i % 20 == 0:
            print(f'\r  {BLUE}checking embeddings - {YELLOW}{i}{BLUE} checked out of {YELLOW}{len(all_chains)}{BLUE}{RESET}', end='', flush=True)

        if len(bindings_flags) != len(embeddings):
            return False, 'Binding sights vector size and the embeddings vector size does not match'

        for embedding in embeddings:
            if len(embedding) != expected_embedding_vector_size:
                return False, f'Unexpected embedding vector size - expected {expected_embedding_vector_size} got {len(embedding)}'
            
        chain.free_embeddings_from_RAM()

    return True, None

# run test cases

def _run(test_case):
    print(f'{BLUE}  ...running {test_case.__name__}{RESET}', end='', flush=True)
    try:
        success, err_message = test_case()
    except Exception as e:
        success = False
        err_message = f'Exception occurred: "{e}"\n\n{traceback.format_exc()}'

    if success:
        print(f'\r{GREEN}[✔] {test_case.__name__}{RESET}{SPACE}')
    else:
        print(f'\r{RED}[✘] {test_case.__name__}{RESET}{SPACE}')
        print(f'  Error: {err_message}')
        
print(f'\n{BLUE}Running tests...{RESET}\n')

all_symbols = globals()
tests = [
    func 
    for func in all_symbols.values() 
    if callable(func) and 
       func.__module__ == __name__ and
       not func.__name__.startswith('_')]

for test_case in tests:
    _run(test_case)

print()