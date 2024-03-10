import json
import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import config.config as config
from typing import Any, Callable, Dict, List, Union
import numpy as np
from numpy import ndarray
from Levenshtein import distance as lev_distance
import datasets_db as datasets
import shutil

ligand = 'MG'
protrusion_fname = '/home/brabecm4/diplomka/protein-binding-sites/data/3d_proc/protrusion.max-neighbors.big.json'
binding_sights_fname = f'/home/brabecm4/diplomka/protein-binding-sites/data/p2rank_ds/binding_sights_{ligand}.txt'
pdbs_sources_fname = f'/home/brabecm4/diplomka/protein-binding-sites/data/p2rank_ds/source_pdbs_for_{ligand}.txt'
relevant_pdbs_folder = f'/home/brabecm4/diplomka/protein-binding-sites/data/p2rank_ds/pdbs_for_{ligand}/'
datasets_folder = '/home/brabecm4/diplomka/protein-binding-sites/data/p2rank_ds/'


db = datasets.SeqDatasetDb()
db.load_protrusion_data_file(protrusion_fname)
ds = db.get_dataset_for(ligand)

test, train = ds.testing(), ds.training()

test  = datasets.Helpers.filter_chains_with_valid_protrusion(test)
train = datasets.Helpers.filter_chains_with_valid_protrusion(train)

output_str = ""

for chain in test + train:
    id = chain.full_id()
    sequence = chain.sequence()
    binding_sights_string = ''.join([f'{val:.0f}' for val in chain.binding_sights()])

    output_str += f'>{id}\n{sequence}\n{binding_sights_string}\n'

with open(binding_sights_fname, 'w') as file:
    file.write(output_str)

file_mappings = {}

for chain in test + train:
    new_pdb_fname = chain.full_id()
    path_to_pdb = chain.get_3d_structure_filename()

    if path_to_pdb.endswith('pdb'):
        new_pdb_fname += '.pdb'
    elif path_to_pdb.endswith('cif'):
        new_pdb_fname += '.cif'
    elif path_to_pdb.endswith('ent'):
        new_pdb_fname += '.ent'
    else:
        raise Exception(f'invalid file name: {path_to_pdb}')

    file_mappings[chain.full_id()] = new_pdb_fname

    shutil.copy(path_to_pdb, relevant_pdbs_folder + new_pdb_fname)

my_datasets = [(train, 'training'), (test, 'testing')]

for ds, ds_label in my_datasets:

    ds_fname = f'{ligand}_{ds_label}.ds'

    output_str = ''

    for chain in ds:
        output_str += f'pdbs_for_{ligand}/{file_mappings[chain.full_id()]}\n'

    with open(datasets_folder + ds_fname, 'w') as file:
        file.write(output_str)
