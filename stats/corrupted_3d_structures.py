import argparse
import os
import sys
from typing import List
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import data_prep.datasets_db as datasets
import data_prep.pdb_files_db as pdb_files_db

RESET = '\033[0m'
RED = '\033[91m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
BLUE = '\033[94m'

# python3 /home/brabecm4/diplomka/protein-binding-sites/stats/corrupted_3d_structures.py

pdb_db = pdb_files_db.PdbFilesDb()
db = datasets.SeqDatasetDb()
db.set_pdb_db(pdb_db)

all_ligands = datasets.SeqDatasetDb.all_ligands()

print()

def print_stats_of(label, valid: List[datasets.ChainRecord], invalid: List[datasets.ChainRecord]):
    print(f'  {label + ":":10}  {YELLOW}{len(valid + invalid):4}{BLUE}  ', end='')
    print(f'({GREEN}{len(valid):4}{BLUE}/{RED}{len(invalid):3}{BLUE})  ', end='')
    print(f' {RED if len(invalid) != 0 else GREEN}{len(invalid) / len(valid + invalid) * 100:.2f} %{BLUE}', end='')

    if len(invalid) == 0:
        print()
        return

    lev_distances = [ch.get_lev_distance_of_3d_sequence_and_sequence() for ch in invalid]
    valid_lev_distances = [dist for dist in lev_distances if dist is not None]

    missing = len(lev_distances) - len(valid_lev_distances)
    corrupted = len(valid_lev_distances)
    avr_lev_dist = sum(valid_lev_distances) / len(valid_lev_distances)

    print(f' ... missing:{RED}{missing:2}{BLUE}, corrupted: {RED}{corrupted:2}{BLUE} (avrg lev: {RED}{avr_lev_dist:.2f}{BLUE})') 


for ligand in all_ligands:
    ds = db.get_dataset_for(ligand)
    
    testing = ds.testing()
    training = ds.training()

    all_chains = ds.all()

    test_valid, test_invalid = datasets.Helpers.split_chains_to_valid_and_invalid_3D_file(testing)
    train_valid, train_invalid = datasets.Helpers.split_chains_to_valid_and_invalid_3D_file(training)
    all_valid, all_invalid = datasets.Helpers.split_chains_to_valid_and_invalid_3D_file(all_chains)

    print(f'{YELLOW}{ligand}{BLUE}')

    print_stats_of('training', train_valid, train_invalid)
    print_stats_of('testing', test_valid, test_invalid)
    print_stats_of(' -> total', all_valid, all_invalid)

    print()

print(f'--- {YELLOW} All data {BLUE} ---\n')

all_chains = db.get_all_chain_records()
valid, invalid = datasets.Helpers.split_chains_to_valid_and_invalid_3D_file(all_chains)

print_stats_of('All chains', valid, invalid)

print(RESET)