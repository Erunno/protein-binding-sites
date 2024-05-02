import json
import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import config.config as config
from typing import Any, Callable, Dict, List, Union
import numpy as np
from numpy import ndarray
from Levenshtein import distance as lev_distance
import data_prep.datasets_db as datasets
import data_prep.pdb_files_db as pdb_files_db
import shutil

# python /home/brabecm4/diplomka/protein-binding-sites/p2rank_prep/p2rank_dataset_preparation.py

# datasets_folder = '/home/brabecm4/diplomka/protein-binding-sites/data/p2rank_ds/test_train_data/'
datasets_folder = '/home/brabecm4/diplomka/protein-binding-sites/p2rank/p2rank-llmex-data'

db = datasets.SeqDatasetDb()
pdb_db = pdb_files_db.PdbFilesDb()
db.set_pdb_db(pdb_db)

def run_for_ligand(ligand):
    working_folder = os.path.join(datasets_folder, '_' + ligand) 
    try:
        os.mkdir(working_folder)
    except:
        pass 

    ds = db.get_dataset_for(ligand)
    test, train = ds.testing(), ds.training()

    test  = datasets.Helpers.filter_chains_with_valid_3D_file(test)
    train = datasets.Helpers.filter_chains_with_valid_3D_file(train)

    output_str = ""

    ids = {}

    for chain in test + train:
        id = chain.full_id()
        
        if id in ids: 
            print(f'ERROR: duplicity found, chain:{chain.full_id()}, lig: {ligand}')

        ids[id] = chain
        
        binding_sights_string = ''.join([f'{val:.0f}' for val in chain.binding_sights()])

        output_str += f'>{id}\n{chain.sequence()}\n{binding_sights_string}\n'

    bindings_fname = f'binding_sights_{ligand}.txt'

    with open(os.path.join(working_folder, bindings_fname), 'w') as file:
        file.write(output_str)

    for type, chains in [['test', test], ['train', train]]:

        content = f'PARAM.RESIDUE_LABELING_FORMAT=sprint\n' + \
                  f'PARAM.RESIDUE_LABELING_FILE={bindings_fname}\n\n' + \
                  f'HEADER:  protein  chains\n\n'
        
        for chain in chains:
            content += f'{chain.get_chain_structure().structure_file} {chain.chain_id()}\n'

        with open(os.path.join(working_folder, f'{ligand}_{type}.ds'), 'w') as file:
            file.write(content)


for ligand in datasets.all_ligands:
    run_for_ligand(ligand)
