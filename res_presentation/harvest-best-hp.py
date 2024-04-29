import json
import os
import table_printer as printer
import argparse
import results_loader 

# python3 /home/brabecm4/diplomka/protein-binding-sites/res_presentation/harvest-best-hp.py --save-to /home/brabecm4/diplomka/protein-binding-sites/data/final_eval/best_HPs.json --results-folder /home/brabecm4/diplomka/protein-binding-sites/data/netw_results/netw_runs

parser = argparse.ArgumentParser(description='Harvest best hyperparameters')
parser.add_argument('--results-folder', help='Folder with results')
parser.add_argument('--save-to', help='File to save the hyperparameters')

args = parser.parse_args()

if args.results_folder:
    results_folder = args.results_folder  
else:
    print("Err: no '--results-folder' option defined")


    
def group_everything_by_tag_ligand(results):
    grouped = {}
    def add_to_results(tag, ligand, item):
        if tag not in grouped:
            grouped[tag] = {}
        
        if ligand not in grouped[tag]:
             grouped[tag][ligand] = [] 

        
        grouped[tag][ligand].append(item)

    for result in results:
        add_to_results(result['result_tag'], result['ligand'], result)

    return grouped

def get_all(key, results):
    return list(dict.fromkeys([res[key] for res in results]))

def get_all_embedders(results):
    return get_all('embedder', results)

def get_all_ligands(results):
    return get_all('ligand', results)

def get_all_tags(results):
    return get_all('result_tag', results)

def get_the_highest_mcc_stat_within_the_result(result):
    return max(result['all_stats'] + [result['final_stats']], key=lambda x: x['mcc'])

def get_run_with_best_mcc(results):
    if (len(results) == 0):
        return None

    best = max(results, key=lambda res: get_the_highest_mcc_stat_within_the_result(res)['mcc'])
    return best

def find_best_hyperparameters(grouped_results):
    best_HP = {}

    for tag in grouped_results:
        tag_HPs = {}

        for ligand in grouped_results[tag]:
            runs = grouped_results[tag][ligand]

            best_run = get_run_with_best_mcc(runs)
            best_epoch_stat = get_the_highest_mcc_stat_within_the_result(best_run)

            tag_HPs[ligand] = {
                'mcc': best_epoch_stat['mcc'],
                'test_mcc': best_epoch_stat['test_data_stats']['mcc'] \
                    if 'test_data_stats' in best_epoch_stat else None,

                'layers': best_run['hidden_layers'],
                'epochs': best_epoch_stat['epochs'],
                'learning_rate': best_run['learning_rate'],
                'radius': best_run['radius'] if 'radius' in best_run else None,
                'threshold': best_epoch_stat['threshold'] if 'threshold' in best_epoch_stat else None,
                'tag': tag,
                'ligand': ligand,
                'seed': best_run['seed'],
                'batch_size': best_run['batch_size'],
                'file_name': best_run['file_name'],
            }

        best_HP[tag] = tag_HPs

    return best_HP

def save_to_file(fname, object):
    with open(fname, 'w') as f:
        json.dump(object, f, indent=4, sort_keys=True)    

def print_all_tags(grouped_runs):
    print ('Tags found:')

    for tag in grouped_runs:
        print (tag)

##########################
#         main           # 
##########################

runs = results_loader.load_from(args.results_folder)
grouped_runs = group_everything_by_tag_ligand(runs)

best_hps = find_best_hyperparameters(grouped_runs)
save_to_file(args.save_to, best_hps)

print_all_tags(grouped_runs)
