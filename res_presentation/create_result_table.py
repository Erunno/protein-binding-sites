import json
import os
import sys
import table_printer as printer
import argparse
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import config.config as config

# python /home/brabecm4/diplomka/protein-binding-sites/res_presentation/create_result_table.py --embedder-aliases ProtT5:-T5 ProtBert:-BERT --compare --tag basic_v6 one_prot_fst_v3_c

parser = argparse.ArgumentParser(description='Description of your script')
parser.add_argument('--tag', nargs='+', help='List of tags')
parser.add_argument('--embedder-aliases', nargs='+', help='List of embedders aliases e.g. ProtT5:-T5')
parser.add_argument('--compare', action='store_true', help='Compare results based on tags.')

args = parser.parse_args()

compare_results = True if args.compare else False
display_result_table = not compare_results
embedder_aliases = [ alias.split(':-') for alias in (args.embedder_aliases if args.embedder_aliases else []) ]
embedder_aliases = { alias[0]: alias[1] for alias in embedder_aliases }

results_folder = config.networks_results_folder  

tags = []
if args.tag:
    tags = args.tag
else:
    print("WARNING: No tags provided (using --tag <tag1> <tag2> ...). ", end='')

    if compare_results:
        print("Comparing all results...")
    else:
        print("Taking all results into consideration...")

def load_results_objects_from(folder, embedder_aliases = [], tags = []):
    results = []

    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        if filename.endswith('.json'):
            with open(file_path, 'r') as file:
                try:
                    result = json.load(file)
                except:
                    print (f'ERROR: cannot read file: {file_path}')
                    exit(1)

                result['file_name'] = filename

                embedder = result['embedder']
                if embedder in embedder_aliases:
                    result['embedder'] = embedder_aliases[embedder]
        
                results.append(result)

    if (len(tags) == 0): 
        return results
    
    return [res for res in results if res['result_tag'] in tags]

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

def get_best_result_for_embedder(results, embedder):
    results_of_embedder = [res for res in results if res['embedder'] == embedder]
    best_run = get_run_with_best_mcc(results_of_embedder)

    if best_run is None: 
        return {
            'mcc val': '---', 
            'mcc test': '---', 
            'fname': '---', 
            'layers': '---',
            'epoch': '---',
            'learning_rate': '---',
            'tag': '---',
            'radius': '---',
            'threshold': '---'
        }

    # here more info about the best run can be extracted
    best_stat = get_the_highest_mcc_stat_within_the_result(best_run)
    best_mcc = best_stat['mcc']
    best_mcc_string = "{:.3f}".format(best_mcc)
    test_mcc_string = "{:.3f}".format(best_stat['test_data_stats']['mcc']) \
        if 'test_data_stats' in best_stat else '---'

    layers = '-'.join([str(layer) for layer in best_run['hidden_layers']])
    epoch = str(best_stat['epochs'])
    learning_rate = str(best_run['learning_rate'])
    tag = best_run['result_tag']

    radius = best_run['radius'] if 'radius' in best_run else '---'
    threshold = best_stat['threshold'] if 'threshold' in best_stat else '---'

    return {
        'mcc val': best_mcc_string, 
        'mcc test': test_mcc_string, 
        'fname': best_run['file_name'], 
        'layers': layers,
        'epoch': epoch,
        'learning_rate': learning_rate,
        'tag': tag,
        'radius': radius,
        'threshold': threshold
    }

def get_line_for_ligand(result_objects, ligand, embedders):
    results_for_ligand = [res for res in result_objects if res['ligand'] == ligand]
    results = [get_best_result_for_embedder(results_for_ligand, e) for e in embedders]

    columns = [
        ligand, 
        embedders, 
        
        [res['mcc val'] for res in results], 
        [res['mcc test'] for res in results], 
        [res['tag'] for res in results], 
        [res['layers'] for res in results], 
        [res['epoch'] for res in results],
        [res['learning_rate'] for res in results],
        [res['radius'] for res in results],
        [res['threshold'] for res in results],

        # [ res['fname'] if len(res['fname']) < 10 else f'{res["fname"][:7]}...{res["fname"][-10:]}' 
        [ res['fname'] for res in results]
    ]

    return columns

def get_comparison_line_for_ligand(results_by_tags, ligand, tags, embedders):
    def get_tag_mcc_values(tag):
        result_objects = results_by_tags[tag]
        results_for_ligand = [res for res in result_objects if res['ligand'] == ligand]
        results = [get_best_result_for_embedder(results_for_ligand, e) for e in embedders]

        return [ res['mcc test'] for res in results ]
    
    return [
        ligand, embedders, *[ get_tag_mcc_values(tag) for tag in tags ]
    ]

##########################
#         main           # 
##########################

if display_result_table:
    results = load_results_objects_from(results_folder, tags=tags, embedder_aliases=embedder_aliases)

    embedders = get_all_embedders(results)
    ligands = get_all_ligands(results)

    table = [ get_line_for_ligand(results, ligand, embedders) for ligand in ligands]
    header = ['ligand', 'embedder', 'mcc val', 'mcc test', 'tag', 'layers', 'epoch', 'learning rate', 'radius', 'threshold', 'file name']

if compare_results:
    all_results_in_folder = load_results_objects_from(results_folder, embedder_aliases) 

    if len(tags) == 0:
        tags = get_all_tags(all_results_in_folder)
    
    results_by_tag = {
        tag: [res for res in all_results_in_folder if res['result_tag'] == tag]
        for tag in tags
    }

    merged_results = []
    for tag in results_by_tag:
        merged_results.extend(results_by_tag[tag])

    embedders = get_all_embedders(merged_results)
    ligands = get_all_ligands(merged_results)

    table = [
        get_comparison_line_for_ligand(results_by_tag, ligand, tags, embedders) 
        for ligand in ligands
    ]

    header = ['ligand', 'embedder' ] + tags


printer.print_table([header] + table)
