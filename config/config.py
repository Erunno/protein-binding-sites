import os

data_top_folder = '/home/brabecm4/diplomka/protein-binding-sites/data'

yu_sequences_folder = f'{data_top_folder}/orig/yu_sequences'
binding_sights_file = f'{data_top_folder}/binding_sights/binding_sights_by_ligand.json'
yu_as_fasta_files_folder = f'{data_top_folder}/all_proteins_of_YU'
proteins_by_datasets_file = f'{data_top_folder}/yu_datasets/proteins_by_datasets.json'
embeddings_folder = f'{data_top_folder}/embedded_sequences'
networks_results_folder = f'{data_top_folder}/netw_results/netw_runs'
model_comparisons_folder = f'{data_top_folder}/netw_results/comparisons'
protrusion_data_file = f'{data_top_folder}/3d_proc/mappings_to_pdbs.json'
graphs_folder = f'{data_top_folder}/graphs'
cache_folder = f'{data_top_folder}/cache_data/production'
pdbs_folder = f'{data_top_folder}/orig/biolip_structures'
