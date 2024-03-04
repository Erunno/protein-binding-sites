import argparse
import os
import sys
from typing import List
import matplotlib.pyplot as plt
import numpy as np
from numpy import ndarray
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import data_prep.datasets_db as datasets
import config.config as config

parser = argparse.ArgumentParser(description='Show corrupted protrusion records stats')
parser.add_argument('--protrusion-file', help='Specify input protrusion', required=True)

# python3 protrusion_distribution.py --protrusion-file /home/brabecm4/diplomka/protein-binding-sites/data/3d_proc/protrusion.max-neighbors.big.json

args = parser.parse_args()
protrusion_fname = args.protrusion_file

def get_file_name_without_extension(file_path):
    base_name = os.path.basename(file_path)
    name, extension = os.path.splitext(base_name)
    return name

def load_histograms():
    db = datasets.SeqDatasetDb()
    db.load_protrusion_data_file(protrusion_fname)

    all_chains = db.get_all_chain_records()
    all_chains = datasets.Helpers.filter_chains_with_valid_protrusion(all_chains)

    all_radii = datasets.Helpers.get_all_radii(all_chains)

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

    return histograms, all_radii

def create_graph(histograms, all_radii):
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
        axes[i // 2, i % 2].bar(x_values_non_binding, non_binding_values, width=bar_width, color='orange', label='non binding')

        axes[i // 2, i % 2].set_title(f"Radius {all_radii[i]}")

    for ax in axes.flat:
        ax.legend()

    plot_fname = f'binding_vs_non_binding_for-{get_file_name_without_extension(protrusion_fname)}.png'
    plt.savefig(os.path.join(config.graphs_folder, plot_fname),
                dpi=300, bbox_inches='tight')

histograms, all_radii = load_histograms()
create_graph(histograms, all_radii)
