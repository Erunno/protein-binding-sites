import json
import re

import numpy as np
import datasets_db
import traceback

# constant parameters
path_to_protrusion_file = '/home/brabecm4/diplomka/protein-binding-sites/data/3d_proc/protrusion.max-neighbors.big.json'
path_to_embeddings_folder = '/home/brabecm4/diplomka/protein-binding-sites/data/embedded_sequences'
path_to_sequences_folder = '/home/brabecm4/diplomka/protein-binding-sites/data/orig/yu_sequences'
tested_embedder = 'ESM'

# constants
ALL_LIGANS = ['ADP', 'AMP', 'ATP', 'CA', 'DNA', 'FE', 'GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN']

RESET = '\033[0m'
RED = '\033[91m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
SPACE = '                       '

# init database

print(f'{YELLOW}Initializing...')
db = datasets_db.SeqDatasetDb(
    sequences_folder=path_to_sequences_folder, 
    embeddings_folder=path_to_embeddings_folder)

print(f'Loading protrusion...{RESET}')
db.load_protrusion_data_file(path_to_protrusion_file)

all_chains = db.get_all_chain_records()

res = { 
    ch.full_id(): 
    { 
        'sequence_from_3d_structure': ch.sequence(),
        'structure_file': ch.get_3d_structure_filename() 
    }
    for ch in all_chains 
}

with open('pdb_files_tests.reference.json', 'w') as json_file:
    json.dump(res, json_file, indent=4)

print ('done')