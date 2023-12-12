import json
import os
import table_printer as printer

results_folder = r'C:\Users\mbrabec\Desktop\MFF\diplomka\neural_netw_emb\data\netw_results'  

ligands = ['ADP', 'AMP', 'ATP', 'CA', 'DNA', 'FE', 'GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN']
embedders = ['BERT', 'ESM', 'T5'] 

def load_results_objects_from(folder):
    results = []
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        if filename.endswith('.json'):
            with open(file_path, 'r') as file:
                result = json.load(file)
                result['file_name'] = filename
                results.append(result)

    return results

def get_the_highest_mcc_stat_within_the_result(result):
    return max(result['all_stats'], key=lambda x: x['mcc'])

def get_run_with_best_mcc(results):
    if (len(results) == 0):
        return None

    best = max(results, key=lambda res: get_the_highest_mcc_stat_within_the_result(res)['mcc'])
    return best

def get_best_mcc_for_embedder(results, embedder):
    results_of_embedder = [res for res in results if res['embedder'] == embedder]
    best_run = get_run_with_best_mcc(results_of_embedder)

    if best_run is None: 
        return '---', ''

    # here more info about the best run can be extracted    
    best_mcc = get_the_highest_mcc_stat_within_the_result(best_run)['mcc']
    best_mcc_string = "{:.3f}".format(best_mcc)
    
    return best_mcc_string, best_run['file_name']

def get_line_for_ligand(results, ligand):
    global embedders

    results_for_ligand = [res for res in results if res['ligand'] == ligand]
    best_mcc_s_with_files = [get_best_mcc_for_embedder(results_for_ligand, e) for e in embedders]

    columns = [
        ligand, 
        embedders, 
        [mcc for mcc, fname in best_mcc_s_with_files], 
        [fname for mcc, fname in best_mcc_s_with_files]
    ]

    return columns
    


results = load_results_objects_from(results_folder)

table = [ get_line_for_ligand(results, ligand) for ligand in ligands ]
header = ['ligand', 'embedder', 'mcc', 'file name']

printer.print_table([header] + table)
