import json
import os
import table_printer as printer
import argparse


ligands = ['ADP', 'AMP', 'ATP', 'CA', 'DNA', 'FE', 'GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN']
embedders = [] # to be defined later in the script 

parser = argparse.ArgumentParser(description='Description of your script')
parser.add_argument('--tag', nargs='+', help='List of tags')
parser.add_argument('--results-folder', help='Folder with results')
args = parser.parse_args()

tags = []

if args.results_folder:
    results_folder = args.results_folder  
else:
    print("Err: no '--results-folder' option defined")

if args.tag:
    tags = args.tag
else:
    print("WARNING: No tags provided (using --tag <tag1> <tag2> ...). Taking all results into consideration...")

def load_results_objects_from(folder):
    global tags

    results = []

    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        if filename.endswith('.json'):
            with open(file_path, 'r') as file:
                result = json.load(file)
                result['file_name'] = filename
                results.append(result)

    if (len(tags) == 0): 
        return results
    
    return [res for res in results if res['result_tag'] in tags]

def get_all_embedders(results):
    return list(dict.fromkeys([res['embedder'] for res in results]))

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
            'epoch': '---',
            'learning_rate': '---',
            'tag': '---'
        }

    # here more info about the best run can be extracted
    best_stat = get_the_highest_mcc_stat_within_the_result(best_run)
    best_mcc = best_stat['mcc']
    best_mcc_string = "{:.3f}".format(best_mcc)
    
    layers = '-'.join([str(layer) for layer in best_run['hidden_layers']])
    epoch = str(best_stat['epochs'])
    learning_rate = str(best_run['learning_rate'])
    tag = best_run['result_tag']

    return {
        'mcc': best_mcc_string, 
        'fname': best_run['file_name'], 
        'layers': layers,
        'epoch': epoch,
        'learning_rate': learning_rate,
        'tag': tag,
    }

def get_line_for_ligand(results, ligand):
    global embedders

    results_for_ligand = [res for res in results if res['ligand'] == ligand]
    results = [get_best_result_for_embedder(results_for_ligand, e) for e in embedders]

    columns = [
        ligand, 
        embedders, 
        
        [res['mcc'] for res in results], 
        [res['tag'] for res in results], 
        [res['layers'] for res in results], 
        [res['epoch'] for res in results],
        [res['learning_rate'] for res in results],

        [ res['fname'] if len(res['fname']) < 10 else f'{res["fname"][:7]}...{res["fname"][-10:]}' 
         for res in results]
    ]

    return columns
    


results = load_results_objects_from(results_folder)
embedders = get_all_embedders(results)

table = [ get_line_for_ligand(results, ligand) for ligand in ligands ]
header = ['ligand', 'embedder', 'mcc', 'tag', 'layers', 'epoch', 'learning rate', 'file name']


printer.print_table([header] + table)
