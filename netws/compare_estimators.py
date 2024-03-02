
import os
import sys

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import config.config as config
import json
import random
import string
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import data_loader as dl
import model_builder as ModelBuilder
import arguments_parser as ArgumentParser
from seed_network import seed_all 
from evaluator import get_statistics
from datetime import datetime
from estimators.basic import BasicNetwork
from mlxtend.evaluate import paired_ttest_5x2cv
from sklearn.metrics import matthews_corrcoef


binding_sights_db_filename = config.binding_sights_file
dataset_db_filename = config.proteins_by_datasets_file
embeddings_top_folder = config.embeddings_folder
results_folder = config.networks_results_folder
protrusion_data_fname = config.protrusion_data_file

# tag = args.tag 
# ligand = args.ligand
ligand = "ZN"
# batch_size = args.batch_size
# n_epochs = args.epochs
# seed = args.seed
# learning_rate = args.learning_rate
# hidden_layers = args.hidden_layers
# stats_interval = args.epoch_stats_interval
# verbose = True if args.verbose else False
# use_simple_model = args.use_simple_model
# embedder = args.embedder
embedder = "T5"
# protrusion_data_fname = args.protrusion_data_file
# pdb_mappings_fname = args.pdb_mappings_fname
# protrusion_radii = args.used_protrusion_radii

seed_all(42)

data_loader = None 
radii_count = 0

print("loading data ...", flush=True)

data_loader = dl.DataLoader(
    binding_sights_db_fname=binding_sights_db_filename,
    dataset_by_ligands_db_fname=dataset_db_filename,
    embeddings_folder=os.path.join(embeddings_top_folder, embedder),
    verbose=True
)

print("parsing data ...", flush=True)

X_train, y_train, X_test, y_test = data_loader.get_data_set_for(ligand)

X_train = torch.tensor(X_train, dtype=torch.float32)
y_train = torch.tensor(y_train, dtype=torch.float32)

X_test = torch.tensor(X_test, dtype=torch.float32)
y_test = torch.tensor(y_test, dtype=torch.float32)

X = torch.cat((X_train, X_test), axis=0)
y = torch.cat((y_train, y_test), axis=0)

print("defining estimators ...", flush=True)

input_size = len(X_train[0])

model_1 = BasicNetwork(
    input_size=input_size,
    epochs=10,
    batch_size=1024,
    hidden_sizes=[512, 256, 64],
    learning_rate=0.001
)
model_1.set_name("model 1")
model_1.set_verbose()

model_2 = BasicNetwork(
    input_size=input_size,
    epochs=6,
    batch_size=1024,
    hidden_sizes=[64, 32],
    learning_rate=0.01
)
model_2.set_name("model 2")
model_2.set_verbose()

def mcc_scorer(model, X, y):
    pred_y = model.predict(X)
    mcc = matthews_corrcoef(y, pred_y)


    def print_ones_and_zeros(label, y):
        zeros_count = np.count_nonzero(y == 0)
        ones_count = np.count_nonzero(y == 1)

        print(f'   {label} ({ones_count}/{zeros_count}) (ones/zeros)')

    print(f'  mcc of "{model.get_name()}" = {mcc}')
    print_ones_and_zeros('reference:', y)
    print_ones_and_zeros('real:     ', pred_y)

    return mcc 

print("evaluating...")

# t, p = paired_ttest_5x2cv(
#     estimator1=model_1, estimator2=model_2, 
#     X=X, y=y.int(), scoring='accuracy', random_seed=1)
t, p = paired_ttest_5x2cv(
    estimator1=model_1, estimator2=model_2, 
    X=X, y=y.int(), scoring=mcc_scorer, random_seed=1)

print('P-value: %.3f, t-Statistic: %.3f' % (p, t))

if p <= 0.05:
    print('Difference between mean performance is probably real')
else:
    print('Algorithms probably have the same performance')

