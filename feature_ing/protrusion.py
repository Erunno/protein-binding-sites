import os
import sys

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import config.config as config
import argparse
import json
import os
from Bio.PDB import PDBParser
from Bio.PDB import MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch
from pathlib import Path
import numpy as np
import pprint
from Levenshtein import distance as lev_distance
from datetime import datetime, timedelta

parser = argparse.ArgumentParser(description='Process input and output directories')
parser.add_argument('--output-file', help='Specify the output file', required=False)
parser.add_argument('--input-dirs', nargs='+', help='Specify the input directories', required=True)
parser.add_argument('--radius', nargs='+', type=float, help='List of radii')

args = parser.parse_args()

output_file = args.output_file
input_dirs = args.input_dirs
radii = args.radius

UP_char = "\033[A"

def get_prot_id(pdb_fname: str):
    if (pdb_fname.endswith('.pdb')):
        return pdb_fname[-9:-5]
    
    if (pdb_fname.endswith('.ent')):
        return pdb_fname[-8:-4]
    
    if (pdb_fname.endswith('.cif')):
        return pdb_fname[-8:-4]

def build_structure(pdb_file: str):
    if (pdb_file.endswith('cif')):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    id = get_prot_id(pdb_file)
    structure = parser.get_structure(id, pdb_file)
    
    structure_atoms = list(structure.get_atoms())
    ns = NeighborSearch(structure_atoms)

    return structure, ns

def search_chain(chain, ns, radius):
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
        'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
        'TYR': 'Y', 'VAL': 'V'
    }

    def map_name(name):
        if name in three_to_one:
            return three_to_one[name]
        return f'#({name})'

    protrusions = []
    sequence = ""
    
    for residue in chain:
        # Exclude hetero residues or non-standard residues with insertion code ' '
        if residue.get_id()[0] == " ":
            sequence += map_name(residue.get_resname())

            max_neighbors_count = 0

            for atom in residue:
                neighbors = ns.search(atom.get_coord(), radius, 'R')
                num_neighbors = len(neighbors)

                if num_neighbors > max_neighbors_count:
                    max_neighbors_count = num_neighbors
        
            protrusions.append(max_neighbors_count)
    
    return protrusions, sequence

def get_neighborhood_counts(pdb_file, radii):
    structure, ns = build_structure(pdb_file)
    protein_id = structure.id

    results_per_chain = []

    for model in structure: 
        for chain in model:

            results_per_radius = []
            
            for radius in radii:
                protrusion, sequence = search_chain(chain, ns, radius)
                results_per_radius.append({
                    'radius': radius,
                    'protrusion': protrusion
                })

            results_per_chain.append({
                'chain': chain.id.upper(),
                'results': results_per_radius,
                'sequence': sequence,
                'origin_file': os.path.abspath(pdb_file)
            })

    return protein_id, results_per_chain

start_time = datetime.now() 
def print_progress(done, total, bar_len, result_object = None):
    estimate = '??? remaining'
    now = datetime.now()
    elapsed_time = now - start_time
    
    if (elapsed_time.total_seconds() > 1):
        progress_rate = done / elapsed_time.total_seconds()
        remaining_time_seconds = (total - done) / progress_rate
        remaining_time = str(timedelta(seconds=int(remaining_time_seconds)))
        estimate = f'eta {remaining_time} remaining'

    result_size = ''
    if (result_object is not None):
        size = sys.getsizeof(result_object)
        result_size = f'(result has {size} bytes)'

    loading_bar_progress = bar_len * done // total
    bar_string = f'[{loading_bar_progress * "="}{(bar_len - loading_bar_progress) * " "}] ({done}/{total}) ... {estimate} {result_size}'

    print('\r', bar_string, end='', flush=True)

def list_files_in_directory(directory):
    path = Path(directory)
    files = [file.name for file in path.iterdir() if file.is_file()]
    return files

def get_all_pdbs(dirs_path):
    result = []

    for dir in dirs_path:
        files = [os.path.join(dir, file) 
                    for file in list_files_in_directory(dir) 
                    if file.endswith('.pdb') or \
                       file.endswith('.cif') or \
                       file.endswith('.ent')
                ]
    
        result = result + files

    return result

def merge_results(res1, res2, protein_id):
    chains1 = [res['chain'] for res in res1]
    chains2 = [res['chain'] for res in res2]

    intersection = set(chains1) & set(chains2)

    if (len(intersection) != 0):
        print(f"\nERR: duplicity found - protein {protein_id} - dup chains {intersection}")

    return res1 + res2

all_pdbs = get_all_pdbs(input_dirs)

done = 0
all_count = len(all_pdbs)
loading_bar_len = 40

result = {}

print ('Starting the protrusion', flush=True)
print_progress(done, all_count, loading_bar_len)

for file_path in all_pdbs:

    protein_id, protrusion_results = get_neighborhood_counts(file_path, radii)

    if protein_id in result:
        protrusion_results = merge_results(protrusion_results, result[protein_id], protein_id)            

    result[protein_id] = protrusion_results

    done += 1
    print_progress(done, all_count, loading_bar_len, result)

result_string = pprint.pformat(result, compact=True, width=100).replace("'",'"')

print ('\n\nDONE')

if output_file is None:
    print(f'\n{result_string}')
else:
    with open(output_file, 'w') as file:
        file.write(result_string)
