import json
import os
import random
import sys

from matplotlib.lines import Line2D
import table_printer as printer
import matplotlib.patches as mpatches
import argparse
import results_loader 
import numpy as np
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import config.config as config

# python /home/brabecm4/diplomka/protein-binding-sites/res_presentation/analyze_performance.py --results-folder /home/brabecm4/diplomka/protein-binding-sites/data/netw_results/netw_runs --save-to-folder /home/brabecm4/diplomka/protein-binding-sites/data/graphs/losses --tag basic_v6

RESET = '\033[0m'
RED = '\033[91m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
MAGENTA = '\033[95m'

title_mappings = {
    'basic_v6': 'Baseline',
    'neighboring_emb_5_v2_c': 'Direct Input - 4 neighbors',
    'nei_emb_5_avrg_v2_c': 'Averaged Embeddings - 4 neighbors',
    'prot_all_c': 'Many Protrusion Values in First Layer',
    'protrusion_bypass_v4_c': 'One Protrusion Value Bypassing',
    'SASA_bypassed_v2_c': 'SASA Value Bypassing',
    'one_prot_fst_v3_c': 'One Protrusion Value in First Layer',
    'nei_emb_3_avrg_v2_c': 'Averaged Embeddings - 2 neighbors',
    'nei_emb_3_v2_c': 'Direct Input - 2 neighbors',
    'SASA_fst_v2_c': 'SASA Value in First Layer',
    'nei_5_comprs_v3_cc': 'Compressing Layer of Networks - 4 neighbors',
    'nei_3_comprs_v3_cc': 'Compressing Layer of Networks - 2 neighbors',
}

parser = argparse.ArgumentParser(description='Harvest best hyperparameters')
parser.add_argument('--results-folder', help='Folder with results')
parser.add_argument('--save-to-folder', help='Folder to save the plots')
parser.add_argument('--tag', type=str, nargs='+', help='Tag of the desired result, special value ALL (takes all known)', required=True)

args = parser.parse_args()

displayed_losses_limit = 500
displayed_mcc_s_limit = 1000


with open(config.best_HPs_file, 'r') as file:
    all_best_hps = json.load(file)
    
def plot_losses(data):
    for losses in data[:displayed_losses_limit]:
        plt.plot(losses, color='lightgray', alpha=1 / 256)

    average_losses = np.mean(data, axis=0)
    median_losses = np.median(data, axis=0)

    max = np.max(average_losses)
    min = np.min(median_losses)

    plt.plot(average_losses, color='orange', linewidth=2, label='Average Loss')
    plt.plot(median_losses, color='teal', linewidth=2, label='Median Loss')

    plt.ylim(0, max + min)

    title = title_mappings[tag] if tag in title_mappings else f'Plot for: {tag}'
    plt.title(title)
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()
    plt.grid(True)
    plt.legend(loc='center right')

    plot_fname = f'losses_of_{tag}.png'
    plt.savefig(os.path.join(args.save_to_folder, plot_fname),
                dpi=100, bbox_inches='tight')
    
    plt.clf()

def get_best_filenames(data):
    ligands = list(set([rec['ligand'] for rec in data]))
    tag = data[0]['result_tag']

    if tag not in all_best_hps:
        return None
    
    fnames = []
    for ligand in ligands:
        if ligand in all_best_hps[tag]: 
            fnames.append(all_best_hps[tag][ligand]['file_name'])

    return fnames if len(fnames) != 0 else None

def plot_mcc_s(data): 
    best_runs_fnames = get_best_filenames(data)
    
    if best_runs_fnames:
        stats_of_best = [rec['all_stats'] + [rec['final_stats']] for rec in data if rec['file_name'] in best_runs_fnames]

        best_test = np.array([np.array([0] + [stat['test_data_stats']['mcc'] for stat in stats]) for stats in stats_of_best])
        best_validation = np.array([np.array([0] + [stat['mcc'] for stat in stats]) for stats in stats_of_best])

    ligands = list(set([rec['ligand'] for rec in data]))
    stats_per_ligand = { lig: [] for lig in ligands }
    for rec in data: 
        stats_per_ligand[rec['ligand']].append(rec['all_stats'] + [rec['final_stats']])
    
    stats_of_all = [rec['all_stats'] + [rec['final_stats']] for rec in data]

    stats_of_all = [rec['all_stats'] + [rec['final_stats']] for rec in data]
    mcc_s_test_all = [np.array([0] + [stat['test_data_stats']['mcc'] for stat in stats]) for stats in stats_of_all]
    mcc_s_test = np.array(mcc_s_test_all)
    
    mcc_s_of_best_test = [np.array([0] + [stat['test_data_stats']['mcc'] for stat in stats]) for stats in stats_of_all]
    mcc_s_of_best_validation = [np.array([0] + [stat['mcc'] for stat in stats]) 
                                for stats in stats_of_all 
                                ]

    all_mcc_s_color = 'blue'
    for mcc_s in mcc_s_test_all[:displayed_mcc_s_limit]:
        alpha = 1 / 50
        plt.plot(mcc_s, color=all_mcc_s_color, alpha=alpha)

    # if best_runs_fnames:
    #     for mcc_s in best_test[:displayed_mcc_s_limit]:
    #         plt.plot(mcc_s, color='orange', alpha=1 / 50)

    #     for mcc_s in best_validation[:displayed_mcc_s_limit]:
    #         plt.plot(mcc_s, color='blue', alpha=1 / 50)

    mcc_s_per_lig_color = 'red'
    for ligand in ligands:
        mcc_s_test_ligand = [np.array([0] + [stat['test_data_stats']['mcc'] for stat in stats]) for stats in stats_per_ligand[ligand]]

        aggr = np.median(mcc_s_test_ligand, axis=0)
        alpha = 0.8

        plt.plot(aggr, color=mcc_s_per_lig_color, alpha=alpha)
    
    
    plt.plot([100, 100], color=mcc_s_per_lig_color, label='MCC median per ligand')
    plt.plot([100, 100], color=all_mcc_s_color, label='MCC of a single model')

    

    test_data_mcc_aggregate = np.median(mcc_s_test, axis=0)
    # plt.plot(test_data_mcc_aggregate, color='red', linewidth=2, label='Median On Test Data')

    # plt.plot(best_test, color='orange', linewidth=2, label='Test Data MCC')
    # plt.plot(best_validation, color='blue', linewidth=2, label='Validation Data MCC')

    max = 1
    min = 0

    plt.ylim(0, max + min)

    title = title_mappings[tag] if tag in title_mappings else f'Plot for: {tag}'
    plt.title(title)
    plt.xlabel('Epoch')
    plt.ylabel('MCC')
    plt.legend()
    plt.grid(True)
    plt.legend(loc='upper right')
    labels = np.arange(len(test_data_mcc_aggregate)) * 10
    plt.xticks(ticks=np.arange(len(test_data_mcc_aggregate)), labels=labels)


    plot_fname = f'mcc_s_of_{tag}.png'
    plt.savefig(os.path.join(args.save_to_folder, plot_fname),
                dpi=100, bbox_inches='tight')
    
    plt.clf()


def plot_for_tag(results, tag):
    print (f'{BLUE}plotting for tag: {YELLOW}{tag:25}{RESET} ', end='')

    results = [res for res in results if res['result_tag'] == tag]
    random.shuffle(results)

    if len(results) == 0:
        print('NO DATA')
        return

    losses_data = [np.array(res['losses']) for res in results]

    plot_losses(losses_data[:displayed_losses_limit])
    print (f'{GREEN}losses  {RESET}', end='')

    plot_mcc_s(results)
    print (f'{GREEN}MCCs{RESET}')
    

# main

results = results_loader.load_from(args.results_folder)

requested_tags = args.tag if args.tag[0] != 'ALL' else title_mappings.keys() 

print ('analyzing following tags: ', requested_tags)

for tag in requested_tags:
    plot_for_tag(results, tag)