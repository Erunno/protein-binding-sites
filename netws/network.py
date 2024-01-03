binding_sights_db_filename = r'C:\Users\mbrabec\Desktop\MFF\diplomka\neural_netw_emb\data\binding_sights\binding_sights_by_ligand.json'
dataset_db_filename = r'C:\Users\mbrabec\Desktop\MFF\diplomka\neural_netw_emb\data\yu_datasets\proteins_by_datasets.json'
embeddings_top_folder = r'C:\Users\mbrabec\Desktop\MFF\diplomka\neural_netw_emb\data\embedded_sequences'
results_folder = r'C:\Users\mbrabec\Desktop\MFF\diplomka\neural_netw_emb\data\netw_results'

import json
import os
import random
import string
import torch
import torch.nn as nn
import torch.optim as optim
import data_loader as dl
import model_builder as ModelBuilder
import arguments_parser as ArgumentParser
from seed_network import seed_all 
from evaluator import get_statistics
from datetime import datetime

args = ArgumentParser.get_args()

tag = args.tag 

ligand = args.ligand
batch_size = args.batch_size
n_epochs = args.epochs
seed = args.seed
learning_rate = args.learning_rate
hidden_layers = args.hidden_layers
stats_interval = args.epoch_stats_interval
verbose = args.verbose
embedder = args.embedder

protrusion_data_fname = args.protrusion_data_file
pdb_mappings_fname = args.pdb_mappings_fname
protrusion_radii = args.used_protrusion_radii

protrusion_used = protrusion_data_fname is not None  

def generate_random_string(length):
    letters = string.ascii_letters + string.digits  # You can include other characters as needed
    return ''.join(random.choice(letters) for _ in range(length))

random_suffix = generate_random_string(5)

result_file_name = f'{ligand}_hl{"-".join([str(n) for n in hidden_layers])}.{datetime.now().strftime("%d-%m-%y-%Hh%Mm%Ss.%f")}.{random_suffix}.json'

seed_all(seed)

data_loader = None 
if protrusion_used:
    data_loader = dl.ProtrusionDataLoader(
        binding_sights_db_fname=binding_sights_db_filename,
        dataset_by_ligands_db_fname=dataset_db_filename,
        embeddings_folder=os.path.join(embeddings_top_folder,embedder),
        protrusion_data_fname=protrusion_data_fname,
        pdb_mappings_fname=pdb_mappings_fname,
        radii=protrusion_radii,
        verbose=verbose
    )
else:
    data_loader = dl.DataLoader(
        binding_sights_db_fname=binding_sights_db_filename,
        dataset_by_ligands_db_fname=dataset_db_filename,
        embeddings_folder=os.path.join(embeddings_top_folder,embedder),
        verbose=verbose
    )

X_train, y_train, X_test, y_test = data_loader.get_data_set_for(ligand)

X_train = torch.tensor(X_train, dtype=torch.float32)
y_train = torch.tensor(y_train, dtype=torch.float32)
 
X_test = torch.tensor(X_test, dtype=torch.float32)
y_test = torch.tensor(y_test, dtype=torch.float32)
 
# define the model
input_size = len(X_train[0])
output_size = 1

model = ModelBuilder.create_model(
    input_size=input_size,
    output_size=output_size,
    hidden_sizes=hidden_layers
)

if (verbose):
    print(str(model))
 
# train the model
loss_fn   = nn.BCELoss()  # binary cross entropy
optimizer = optim.Adam(model.parameters(), lr=learning_rate)

all_stats = []

for epoch in range(n_epochs):
    for i in range(0, len(X_train), batch_size):
        X_batch = X_train[i:i+batch_size]
        
        y_pred = model(X_batch)[:, 0]
        y_batch = y_train[i:i+batch_size]
        
        loss = loss_fn(y_pred, y_batch)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    if (verbose):
        print(f'Epoch: {epoch}, loss: {loss :.5f}')

    if (epoch % stats_interval == 0):
        stats = get_statistics(model, X_test, y_test, epoch, print_res=verbose)
        all_stats.append(stats)


final_stats = get_statistics(model, X_test, y_test, n_epochs, print_res=verbose)

all_stats.append(final_stats)

training_report = {
    'result_tag': tag,
    'ligand': ligand,
    'embedder': embedder,
    'model_to_string': str(model),
    'batch_size': batch_size,
    'total_epochs': n_epochs,
    'seed': seed,
    'learning_rate': learning_rate,
    'hidden_layers': hidden_layers,
    'all_stats': all_stats,
    'final_stats': final_stats
}

if protrusion_used:
    training_report['protrusion_radii'] = data_loader.radii

result_path = os.path.join(results_folder, result_file_name)

with open(result_path, 'w') as file:
    json.dump(training_report, file, indent=2)