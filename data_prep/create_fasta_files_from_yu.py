import os

yu_ds_folder = r'C:\Users\mbrabec\Desktop\MFF\diplomka\data\yu_sequences'
out_folder = r'.\data\all_proteins_of_YU'

proteins_in_one_fasta_file_count = 500

protein_id_column_idx = 0
protein_sequence_idx = 5
column_separator = ';'

def get_sequences_from_file(fname): 
    with open(fname, 'r') as file:
        protein_records = [
            line_record.split(column_separator)
            for line_record in file.read().replace('\r', '').split('\n')
            if line_record != ''
        ]

    return { record[protein_id_column_idx].lower(): record[protein_sequence_idx] for record in protein_records }

def get_all_files_paths_of_yu(yu_ds_folder):
    return [
        os.path.join(root, file_name)
        for root, _, files in os.walk(yu_ds_folder)
        for file_name in files
        if file_name.endswith('.txt')
    ]

def get_fasta_records_from_sequences(sequences):
    return [
        f'''>{id}
{sequences[id]}
'''
        for id in sequences
    ]

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

yu_file_paths = get_all_files_paths_of_yu(yu_ds_folder)

all_sequences = {}
for yu_file_path in yu_file_paths:
    sequences_from_one_file = get_sequences_from_file(yu_file_path)
    all_sequences.update(sequences_from_one_file)

fasta_records = get_fasta_records_from_sequences(all_sequences)

fasta_recs_splitted = split_list(fasta_records, proteins_in_one_fasta_file_count)

for idx, fasta_rec_chunk in enumerate(fasta_recs_splitted):
    print(f'writing chunk {idx}')
    out_fname = os.path.join(out_folder, f'yu_{idx}.fasta')
    write_to_fasta_file(out_fname, fasta_rec_chunk)

print('DONE')
