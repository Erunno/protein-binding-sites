import argparse
import os
import sys
from typing import List
import matplotlib.pyplot as plt
import numpy as np
from numpy import ndarray
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import data_prep.datasets_db as datasets
import data_prep.pdb_files_db as pdb_structures
import config.config as config

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 12

parser = argparse.ArgumentParser(description='Show corrupted protrusion records stats')
parser.add_argument('--out-graph-path', help='Specify output file', required=True)

# python3 /home/brabecm4/diplomka/protein-binding-sites/stats/protrusion_distribution_per_ligand.py --out-graph-path /home/brabecm4/diplomka/protein-binding-sites/data/graphs/binding_vs_non_binding_per_ligand.ATP_FE.png

args = parser.parse_args()

def load_histograms():
    pdb_db = pdb_structures.PdbFilesDb()
    db = datasets.SeqDatasetDb()
    db.set_pdb_db(pdb_db)    

    all_radii = [5.5, 8.5]
    all_ligands = ["FE", "ATP"] # has to be 2 ligands

    hists_per_ligand = {}

    for ligand in all_ligands:
        all_chains = db.get_dataset_for(ligand).all()
        all_chains = datasets.Helpers.filter_chains_with_valid_protrusion(all_chains)

        binding_sights = datasets.Helpers.concat_chain_data(
            datasets.DataAccessors.biding_sights_vect(),
            chains=all_chains
        )

        def get_histogram(vect):
            histogram = {}

            for val in vect:
                histogram[val] = histogram[val] + 1 if val in histogram else 1

            return histogram

        histograms = []

        for radius in all_radii:
            print('working on radius: ', radius)
            entire_protrusion = datasets.Helpers.concat_chain_data(
                datasets.DataAccessors.protrusion(radius),
                chains=all_chains
            )

            binding_sights_protrusion, non_binding_sights_protrusion = \
                entire_protrusion[binding_sights == 1], entire_protrusion[binding_sights == 0] 

            histograms.append({
                'binding': get_histogram(binding_sights_protrusion),
                'non_binding': get_histogram(non_binding_sights_protrusion),
            })

        hists_per_ligand[ligand] = histograms

    used_histograms = []
    used_radii = []
    used_ligands = []

    for i, radius in enumerate(all_radii):
        used_histograms = used_histograms + [
            hists_per_ligand[all_ligands[0]][i], 
            hists_per_ligand[all_ligands[1]][i] 
        ]

        used_radii = used_radii + [radius, radius]
        used_ligands = used_ligands + [all_ligands[0], all_ligands[1]]
        

    return used_histograms, used_radii, used_ligands

def create_graph(histograms, all_radii, all_ligands):
    bar_width = 0.3

    def normalize(histogram):
        total_items = sum(histogram.values())
        for key in histogram:
            histogram[key] = histogram[key] / total_items 

    for radius_hist in histograms:
        normalize(radius_hist['binding'])
        normalize(radius_hist['non_binding'])

    num_rows = len(histograms)
    fig, axes = plt.subplots((num_rows + 1) // 2, 2, figsize=(30, 5 * num_rows))

    for i, entry in enumerate(histograms):
        x_values = sorted(set(entry['binding'].keys()).union(entry['non_binding'].keys()))

        binding_values = [entry['binding'].get(x, 0) for x in x_values]
        non_binding_values = [entry['non_binding'].get(x, 0) for x in x_values]

        x_values_binding = np.array(x_values) - bar_width / 2
        x_values_non_binding = np.array(x_values) + bar_width / 2

        axes[i // 2, i % 2].bar(x_values_binding, binding_values, width=bar_width, color='blue', label='binding')
        axes[i // 2, i % 2].bar(x_values_non_binding, non_binding_values, width=bar_width, color='orange', label='non-binding')

        axes[i // 2, i % 2].set_title(f"Ligand {all_ligands[i]}, Radius {all_radii[i]} Å", fontsize=44)

    for ax in axes.flat:
        ax.legend(fontsize=30)

    plot_fname = args.out_graph_path
    plt.savefig(os.path.join(config.graphs_folder, plot_fname),
                dpi=100, bbox_inches='tight')

histograms, all_radii, all_ligands = load_histograms()
create_graph(histograms, all_radii, all_ligands)
