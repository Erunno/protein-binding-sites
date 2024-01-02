
import argparse
import os
from Bio.PDB import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch
from pathlib import Path

parser = argparse.ArgumentParser(description='Process input and output directories')
parser.add_argument('--output-file', help='Specify the output file', required=False)
parser.add_argument('--input-dir', help='Specify the input directory', required=True)

args = parser.parse_args()

output_file = args.output_file
input_dir = args.input_dir


def pdb_to_fasta(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("my_structure", pdb_file)

    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
        'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
        'TYR': 'Y', 'VAL': 'V'
    }

    sequence = ""
    count = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                count += 1
                if residue.get_id()[0] == " ":
                    # Exclude hetero residues or non-standard residues with insertion code ' '
                    long_name = residue.get_resname()
                    short_name = three_to_one[long_name]

                    sequence += short_name

    return sequence

def list_files_in_directory(directory):
    path = Path(directory)
    files = [file.name for file in path.iterdir() if file.is_file()]
    return files

def print_progress(done, total, bar_len):
    loading_bar_progress = bar_len * done // total 
    print(f'\r [{loading_bar_progress * "="}{(bar_len - loading_bar_progress) * " "}]', end='')


fasta_string = ''

all_pdbs = [f for f in list_files_in_directory(input_dir) if f.endswith('.pdb')]

done = 0
all_count = len(all_pdbs)
loading_bar_len = 40

for file in all_pdbs:

    print_progress(done, all_count, loading_bar_len)

    full_pdb_path = os.path.join(input_dir, file)

    protein_name = file.split('.')[0]
    protein_record = f'>{protein_name}\n{pdb_to_fasta(full_pdb_path)}\n'

    fasta_string += protein_record
    done += 1

print_progress(done, all_count, loading_bar_len)

if output_file is None:
    print('\n' + fasta_string)
else:
    with open(output_file, 'w') as file:
        file.write(fasta_string)
