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

def get_best_result_for_embedder(results, embedder):
    results_of_embedder = [res for res in results if res['embedder'] == embedder]
    best_run = get_run_with_best_mcc(results_of_embedder)

    if best_run is None: 
        return {
            'mcc': '---', 
            'fname': '---', 
            'layers': '---',
            'epoch': '---'
        }

    # here more info about the best run can be extracted
    best_stat = get_the_highest_mcc_stat_within_the_result(best_run)
    best_mcc = best_stat['mcc']
    best_mcc_string = "{:.3f}".format(best_mcc)
    
    layers = '-'.join([str(layer) for layer in best_run['hidden_layers']])
    epoch = str(best_stat['epochs'])
    learning_rate = str(best_run['learning_rate'])

    return {
        'mcc': best_mcc_string, 
        'fname': best_run['file_name'], 
        'layers': layers,
        'epoch': epoch,
        'learning_rate': learning_rate,
    }

def get_line_for_ligand(results, ligand):
    global embedders

    results_for_ligand = [res for res in results if res['ligand'] == ligand]
    results = [get_best_result_for_embedder(results_for_ligand, e) for e in embedders]

    columns = [
        ligand, 
        embedders, 
        
        [res['mcc'] for res in results], 
        [res['layers'] for res in results], 
        [res['epoch'] for res in results],
        [res['learning_rate'] for res in results],

        [ res['fname'] if len(res['fname']) < 10 else f'{res["fname"][:7]}...{res["fname"][-10:]}' 
         for res in results]
    ]

    return columns
    


results = load_results_objects_from(results_folder)

table = [ get_line_for_ligand(results, ligand) for ligand in ligands ]
header = ['ligand', 'embedder', 'mcc', 'layers', 'epoch', 'learning rate', 'file name']

printer.print_table([header] + table)
