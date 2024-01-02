import argparse
import json
import os
from Bio.PDB import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch
from pathlib import Path
import numpy as np
import pprint

file_path = r"C:\Users\mbrabec\Desktop\MFF\diplomka\data\biolip_structures\receptor\1a0bA.pdb"

parser = argparse.ArgumentParser(description='Process input and output directories')
parser.add_argument('--output-file', help='Specify the output file', required=False)
parser.add_argument('--input-dir', help='Specify the input directory', required=True)
parser.add_argument('--radius', nargs='+', type=float, help='List of radii')

args = parser.parse_args()

output_file = args.output_file
input_dir = args.input_dir
radii = args.radius

def build_structure(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("my_structure", pdb_file)

    structure_atoms = list(structure.get_atoms())
    ns = NeighborSearch(structure_atoms)

    return structure, ns

def search_structure(structure, ns, radius):
    neighbors_counts = []
    
    for model in structure:
        for chain in model:
            for residue in chain:

                # Exclude hetero residues or non-standard residues with insertion code ' '
                if residue.get_id()[0] == " ":
                
                    center = residue.center_of_mass()
                    
                    neighbors = ns.search(center=center, radius=radius, level='R') 
                    count = len(neighbors)

                    neighbors_counts.append(count)
    
    return neighbors_counts

def get_neighborhood_counts(pdb_file, radii):
    structure, ns = build_structure(pdb_file)

    results_per_radius = []

    for radius in radii:
        results_per_radius.append({
            'radius': radius,
            'neighbors': search_structure(structure, ns, radius),
        }) 

    return results_per_radius

def print_progress(done, total, bar_len):
    loading_bar_progress = bar_len * done // total 
    print(f'\r [{loading_bar_progress * "="}{(bar_len - loading_bar_progress) * " "}] ({done}/{total})', end='')

def list_files_in_directory(directory):
    path = Path(directory)
    files = [file.name for file in path.iterdir() if file.is_file()]
    return files

all_pdbs = [f for f in list_files_in_directory(input_dir) if f.endswith('.pdb')]

done = 0
all_count = len(all_pdbs)
loading_bar_len = 40

result = {}

for file in all_pdbs:
    full_pdb_path = os.path.join(input_dir, file)

    protrusion = get_neighborhood_counts(full_pdb_path, radii)

    protein_id = file.split('.')[0]
    result[protein_id] = protrusion 

    done += 1
    print_progress(done, all_count, loading_bar_len)

result_string = pprint.pformat(result, compact=True).replace("'",'"')

if output_file is None:
    print(f'\n{result_string}')
else:
    with open(output_file, 'w') as file:
        file.write(result_string)
