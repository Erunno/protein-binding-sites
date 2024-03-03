import argparse
import os
import sys
from typing import List
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import data_prep.datasets_db as datasets

parser = argparse.ArgumentParser(description='Show corrupted protrusion records stats')
parser.add_argument('--protrusion-file', help='Specify input protrusion', required=True)

args = parser.parse_args()
protrusion_fname = args.protrusion_file

db = datasets.SeqDatasetDb()
db.load_protrusion_data_file(protrusion_fname)

all_chains = db.get_all_chain_records()
all_chains = datasets.Helpers.filter_chains_with_valid_protrusion(all_chains)

def split_vector(flag_vector, data_vector):
    zero, one = [], []
    for i in range(flag_vector):
        if i == 0:
            zero.append(data_vector[i])
        elif i == 1:     
            one.append(data_vector[i])     
        else:
            raise Exception(f'invalid flag vector value {data_vector[i]}')

def to_data(chains):
    


