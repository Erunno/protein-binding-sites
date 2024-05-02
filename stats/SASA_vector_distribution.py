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

parser = argparse.ArgumentParser(description='Show corrupted protrusion records stats')
parser.add_argument('--out-graph-path', help='Specify output file', required=True)

# python3 /home/brabecm4/diplomka/protein-binding-sites/stats/SASA_vector_distribution.py --out-graph-path /home/brabecm4/diplomka/protein-binding-sites/data/graphs/binding_vs_non_binding.SASA.png

args = parser.parse_args()

def load_histogram():
    pdb_db = pdb_structures.PdbFilesDb()
    db = datasets.SeqDatasetDb()
    db.set_pdb_db(pdb_db)    

    all_chains = db.get_all_chain_records_with_merged_binding_sites()
    all_chains = datasets.Helpers.filter_chains_with_valid_3D_file(all_chains)

    binding_sights = datasets.Helpers.concat_chain_data(
        datasets.DataAccessors.biding_sights_vect(),
        chains=all_chains
    )

    SASA_vectors = datasets.Helpers.concat_chain_data(
        datasets.DataAccessors.SASA_vector(),
        chains=all_chains
    )

    binding_sights_SASA, non_binding_sights_SASA = \
        SASA_vectors[binding_sights == 1], SASA_vectors[binding_sights == 0] 

    return {
        'binding': binding_sights_SASA,
        'non_binding': non_binding_sights_SASA,
    }

def create_graph(histogram):
    bins=50
    alpha=1
    hist_type='step'
    
    plt.hist(histogram['non_binding'], bins=bins, alpha=alpha, label='non-binding', color='orange', density=True, histtype=hist_type)
    plt.hist(histogram['binding'], bins=bins, alpha=alpha, label='binding', color='blue', density=True, histtype=hist_type)

    # Add labels and title
    plt.xlabel('Accessible Surface', fontsize=17)
    plt.ylabel('Frequency', fontsize=17)
    plt.title('SASA Distributions', fontsize=20)

    # Add legend
    plt.legend(fontsize=13)

    plt.gca().set_xticks([])  # Hide x-axis numbers
    plt.gca().set_yticks([])  # Hide y-axis numbers

    plot_fname = args.out_graph_path
    plt.savefig(os.path.join(config.graphs_folder, plot_fname),
                dpi=200, bbox_inches='tight')

histogram = load_histogram()
create_graph(histogram)
