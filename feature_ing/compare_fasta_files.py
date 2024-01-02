import argparse
import math
import pprint
import sys
from Levenshtein import distance as lev_distance

UNSPECIFIED_CHAIN = '<unspecified chain>'

parser = argparse.ArgumentParser(description='Process input and output directories')
parser.add_argument('file1', help='File to compare')
parser.add_argument('file2', help='File to compare')

parser.add_argument('--mappings-output-file', help='File to compare')

args = parser.parse_args()

f1 = args.file1
f2 = args.file2

mappings_out_file = args.mappings_output_file

def parse_fasta_file(file_path):
    sequences = {}
    with open(file_path, 'r') as file:
        lines = [ line.strip() for line in file.readlines() if line.strip() != '']
        
    current_id = None
    current_chain = None

    for line in lines:
        if line.startswith('>'):
            full_id = line[1:]

            current_id = line[1:][:4]
            current_chain = full_id[4] if len(full_id) == 5 else UNSPECIFIED_CHAIN

        else:
            if current_id not in sequences:
                sequences[current_id] = {}

            sequences[current_id][current_chain] = line

    return sequences

def get_missing_keys(base_set, database):
    return [key for key in base_set if key not in database ]

def get_present_keys(base_set, database):
    return [key for key in base_set if key in database ]

def get_best_matches(protein_id, base_set, database):
    base_record = base_set[protein_id]
    compare_record = database[protein_id]

    best_matches = {}

    for base_prot_chain in base_record:
        
        best_edit_distance = sys.maxsize
        best_match_chain = None

        for compare_prot_chain in compare_record:
            base_seq = base_record[base_prot_chain]
            compare_seq = compare_record[compare_prot_chain]

            edit_distance = lev_distance(base_seq, compare_seq)

            if edit_distance < best_edit_distance:
                best_edit_distance = edit_distance
                best_match_chain = compare_prot_chain

        best_matches[base_prot_chain] = best_match_chain

    return best_matches


base_set = parse_fasta_file(f1)
database = parse_fasta_file(f2)

missing_keys = get_missing_keys(base_set, database)
present_keys = get_present_keys(base_set, database)

print (f'present: {len(present_keys)}, missing: {len(missing_keys)}')

mappings = {}

for protein_id in present_keys:
    best_matches = get_best_matches(protein_id, base_set, database)
    
    protein_chains = base_set[protein_id]

    for chain in protein_chains:
        best_matching_chain = best_matches[chain]
        chain_text = f' <chain {chain}>' if chain != UNSPECIFIED_CHAIN else ''

        base_seq = protein_chains[chain]
        compare_seq = database[protein_id][best_matching_chain]   

        protein_chain_id = protein_id + (chain if chain != UNSPECIFIED_CHAIN else '')

        mappings[protein_chain_id] = {
            'distance': lev_distance(base_seq, compare_seq),
            'chain_id': best_matching_chain
        }

matching = [key for key in mappings if mappings[key]['distance'] == 0]
different = [key for key in mappings if mappings[key]['distance'] != 0]

print (f'matching sequences: {len(matching)}, non matching sequences: {len(different)}\n')

if len(different) != 0:
    print ('Non matching:')
    
    for protein_id in different:
        best_match = mappings[protein_id]

        print (f'{protein_id}: {best_match["distance"]} <chain {best_match["chain_id"]}>')


# print('\nMissing proteins:')

# for protein_id in missing_keys:
#     print (protein_id)

result = {
    'mappings': mappings,
    'missing': missing_keys
}

result_string = pprint.pformat(result, compact=True).replace("'",'"')

if mappings_out_file is None:
    print(f'\n{result_string}')
else:
    with open(mappings_out_file, 'w') as file:
        file.write(result_string)
