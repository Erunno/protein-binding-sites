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

# python /home/brabecm4/diplomka/protein-binding-sites/p2rank_prep/SAS_points_ds_prep.py

ds_folder = '/home/brabecm4/diplomka/protein-binding-sites/data/p2rank_ds/ds_per_chain/'

db = datasets.SeqDatasetDb()
pdb_db = pdb_files_db.PdbFilesDb()
db.set_pdb_db(pdb_db)


all_chains = db.get_all_chain_records()

##################
# datasets (.ds) #
##################

for chain in all_chains:
    fname = chain.get_chain_structure().structure_file

    content = f'PARAM.RESIDUE_LABELING_FORMAT=sprint\n' + \
              f'PARAM.RESIDUE_LABELING_FILE=/home/brabecm4/diplomka/protein-binding-sites/data/p2rank_ds/ds_per_chain/__binds.txt\n\n' + \
              f'HEADER:  protein  chains\n\n' + \
              f'{fname} {chain.chain_id()}\n'

    with open(ds_folder + chain.full_id() + '.ds', "w") as file:
        file.write(content)

##################
# mock binding   #
##################

mock_bind=''    
for chain in all_chains:
    mock_bind += f'>{chain.full_id()}\n{chain.sequence()}\n{"1" + "0" * (len(chain.sequence()) - 1)}\n\n'

with open(ds_folder + '__binds.txt', "w") as file:
    file.write(mock_bind)
