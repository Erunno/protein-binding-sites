import argparse
import os
import sys
from typing import List
import matplotlib.pyplot as plt
import numpy as np
from numpy import ndarray
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import data_prep.datasets_db as datasets
import data_prep.pdb_files_db as pdb_structures
import config.config as config


# python3 /home/brabecm4/diplomka/protein-binding-sites/stats/yu_ds_summary.py 

pdb_db = pdb_structures.PdbFilesDb()
db = datasets.SeqDatasetDb()
db.set_pdb_db(pdb_db) 

def get_stats_for(ligand):
    ds = db.get_dataset_for(ligand)

    train, test = ds.training(), ds.testing()
    train = datasets.Helpers.filter_chains_with_valid_protrusion(train)
    test = datasets.Helpers.filter_chains_with_valid_protrusion(test)

    train_chain_count = len(train)
    test_chain_count = len(test)
    total_chain_count = train_chain_count + test_chain_count

    train_residue_count = sum(len(chain.sequence()) for chain in train)
    test_residue_count = sum(len(chain.sequence()) for chain in test)
    total_residue_count = train_residue_count + test_residue_count

    train_binding = sum(sum(chain.binding_sights()) for chain in train)
    train_non_binding = train_residue_count - train_binding

    test_binding = sum(sum(chain.binding_sights()) for chain in test)
    test_non_binding = test_residue_count - test_binding

    total_binding = sum(sum(chain.binding_sights()) for chain in train + test)
    total_non_binding = total_residue_count - total_binding
    
    # print('ligand: ', ligand ,
    #       ', total chains:', total_chain_count, 
    #       ', train chains:', train_chain_count, f'({train_chain_count/total_chain_count * 100:.2f}%)',
    #       ', test chains:', test_chain_count, f'({test_chain_count/total_chain_count * 100:.2f}%)',

    #       ', total residues:', total_residue_count,
    #       ', train residues:', train_residue_count, f'({train_residue_count/total_residue_count * 100:.2f}%)',
    #       ', test residues:', test_residue_count, f'({test_residue_count/total_residue_count * 100:.2f}%)',

    #       ', total binding:', total_binding, f'({total_binding/total_residue_count * 100:.2f}%)',
    #       ', total non-binding:', total_non_binding, f'({total_non_binding/total_residue_count * 100:.2f}%)'
    #       )

    # columns: ligand, type, chains, residues, binding, non-binding
    l = '{'
    r = '}'
    print(f'''
        {ligand} & Training & {train_chain_count:,} & \\textcolor{l}prct{r}{l}{train_chain_count/total_chain_count * 100:.2f}{r} & {train_residue_count:,} & \\textcolor{l}prct{r}{l}{train_residue_count/total_residue_count * 100:.2f}{r} & {train_binding:,} & \\textcolor{l}prct{r}{l}{train_binding/train_residue_count*100:.2f}{r} & {train_non_binding:,} & \\textcolor{l}prct{r}{l}{train_non_binding/train_residue_count*100:.2f}{r} \\\\
                 & Testing  & {test_chain_count:,} & \\textcolor{l}prct{r}{l}{test_chain_count/total_chain_count * 100:.2f}{r} & {test_residue_count:,} & \\textcolor{l}prct{r}{l}{test_residue_count/total_residue_count * 100:.2f}{r} & {test_binding:,} & \\textcolor{l}prct{r}{l}{test_binding/test_residue_count*100:.2f}{r} & {test_non_binding:,} & \\textcolor{l}prct{r}{l}{test_non_binding/test_residue_count*100:.2f}{r} \\\\
                 & Total    & {total_chain_count:,} & - & {total_residue_count:,} & - & {total_binding:,} & \\textcolor{l}prct{r}{l}{total_binding/total_residue_count * 100:.2f}{r} & {total_non_binding:,} & \\textcolor{l}prct{r}{l}{total_non_binding/total_residue_count * 100:.2f}{r} \\\\
        \cmidrule(lr)''', '{1-10}')


all_ligands = ['ADP', 'AMP', 'ATP', 'CA', 'DNA', 'FE', 'GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN']
for ligand in all_ligands:
    get_stats_for(ligand)

# get_stats_for("AMP")