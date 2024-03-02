import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import config.config as config
import datasets_db


yu_ds_folder = config.yu_sequences_folder
out_folder = config.yu_as_fasta_files_folder

proteins_in_one_fasta_file_count = 500

def write_to_fasta_file(fname, fasta_records):
    with open(fname, 'w') as f:
        for record in fasta_records:
            f.write(record)

def split_list(list, items_per_list):
    list_len = len(list)
    chunk_count = (list_len + items_per_list - 1) // items_per_list
    
    return [ 
        list[chunk_idx * items_per_list : (chunk_idx + 1) * items_per_list]
        for chunk_idx in range(chunk_count)
    ]

db = datasets_db.SeqDatasetDb()
all_chain_records = db.get_all_chain_records()

fasta_records = [chain.to_fasta_record() for chain in all_chain_records]

fasta_recs_splitted = split_list(fasta_records, proteins_in_one_fasta_file_count)

for idx, fasta_rec_chunk in enumerate(fasta_recs_splitted):
    print(f'writing chunk {idx}')
    out_fname = os.path.join(out_folder, f'yu_{idx}.fasta')
    write_to_fasta_file(out_fname, fasta_rec_chunk)

print('DONE')
