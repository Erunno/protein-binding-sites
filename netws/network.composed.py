import argparse
from datetime import datetime
import json
import os
import random
import sys

import torch
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import config.config as config
from estimators.basic import BasicNetwork
from estimators.bypass import BypassedInputsNetwork
import data_prep.datasets_db as dataset
import data_prep.pdb_files_db as pdb_data
from evaluator import get_statistics
import string
from seed_network import seed_all 
import torch.nn as nn
from sklearn.preprocessing import StandardScaler

allowed_ligands = dataset.SeqDatasetDb.all_ligands()
# allowed_embedder = ['BERT', 'ESM', 'T5']
allowed_embedder = ['ESM']

parser = argparse.ArgumentParser(description='Ligand binding sites neural network')

# srun -p gpu-short --gres=gpu:V100 -A nprg058s --cpus-per-task=4 --mem-per-cpu=32G python /home/brabecm4/diplomka/protein-binding-sites/netws/network.composed.py --verbose True --seed 42 --learning-rate 0.01 --epochs 10 --batch-size 1000 --embedder ESM --epoch-stats-interval 10 --tag test --hidden-layers 256 256 256 32  --hidden-layers-sec 1024 32 5 --epochs-sec 10 --learning-rate-sec 0.01 --neighbors 6 --ligand MN
 
parser.add_argument('--tag', type=str, help='Tag the final result', required=True)
parser.add_argument('--ligand', type=str, choices=allowed_ligands, help='Name of the ligand')
parser.add_argument('--embedder', type=str, choices=allowed_embedder, help='Embedder to be used.')
parser.add_argument('--hidden-layers', type=int, nargs='+', help='List of hidden layer sizes')
parser.add_argument('--hidden-layers-sec', type=int, nargs='+', help='List of hidden layer sizes')
parser.add_argument('--batch-size', type=int, help='Size of one batch')
parser.add_argument('--epochs', type=int, help='Number of epochs')
parser.add_argument('--epochs-sec', type=int, help='Number of epochs')
parser.add_argument('--epoch-stats-interval', type=int, help='Interval of gathering statistics.')
parser.add_argument('--seed', type=int, help='Seed of random.')
parser.add_argument('--neighbors', type=int, help='Seed of random.')
parser.add_argument('--learning-rate', type=float, help='Learning rate')
parser.add_argument('--learning-rate-sec', type=float, help='Learning rate')
parser.add_argument('--verbose', type=bool, default=False, help='Print intermediate results.')

# parser.add_argument('--use-simple-model', type=bool, default=False, help='Use simple sequential model.')
# parser.add_argument('--protrusion-data-file', type=str, help='Path to protrusion data')
# parser.add_argument('--pdb-mappings-fname', type=str, help='Path to mappings to pdb files')
# parser.add_argument('--used-protrusion-radii', type=float, nargs='+', help='List of used radii for protrusion. (only those found in "--protrusion-data-file" can be specified)')


args = parser.parse_args()
print('args: ', args)

args.verbose = True if args.verbose else False

seed_all(args.seed)

def generate_random_string(length):
    letters = string.ascii_letters + string.digits
    return ''.join(random.choice(letters) for _ in range(length))

random_suffix = generate_random_string(5)
result_file_name = f'{args.ligand}_hl{"-".join([str(n) for n in args.hidden_layers])}.{datetime.now().strftime("%d-%m-%y-%Hh%Mm%Ss.%f")}.{random_suffix}.json'


print ('loading data ...', flush=True)

pdb_db = pdb_data.PdbFilesDb()

db = dataset.SeqDatasetDb()
db.set_pdb_db(pdb_db)

ds = db.get_dataset_for(args.ligand)

print ('defining training and testing data ...', flush=True)

X_train, y_train, X_test, y_test = ds.get_train_test_data(
    # [dataset.DataAccessors.average_neighborhood_embeddings(args.embedder, 5)],
    # [dataset.DataAccessors.neighborhood_embeddings(args.embedder, 5)],
    [dataset.DataAccessors.embeddings(args.embedder)],

     filters=[dataset.Helpers.filter_chains_with_valid_protrusion]
)

print ('initializing network ...', flush=True)

model = BasicNetwork(batch_size=args.batch_size, 
                     input_size=len(X_train[0]),
                     hidden_sizes=args.hidden_layers,
                     epochs=args.epochs,
                     learning_rate=args.learning_rate)

all_stats_in_first = []
first_losses = []

def first_epoch_callback(epoch, loss, model):
    if (args.verbose):
        print(f'Epoch: {epoch}, loss: {loss :.5f}', flush=True)

    first_losses.append(float(loss))

    if (epoch % args.epoch_stats_interval == 0):
        stats = get_statistics(model, 
                               X_test, y_test, epoch, print_res=args.verbose)
        all_stats_in_first.append(stats)

model.register_epoch_callback(first_epoch_callback)

print ('fitting ... the first model', flush=True)
model.fit(X_train, y_train)


print ('\nfinal of first network: ')
final_stats_of_first = get_statistics(
    model, 
    X_test, y_test, args.epochs, print_res=args.verbose)


print ('\ncreating custom embeddings...')

sequential_model = model.underling_model()
layers = list(sequential_model.children())
trained_embedder = nn.Sequential(*layers[:-1])

def transform_embeddings(embeddings):
    X = StandardScaler().fit_transform(embeddings)
    X = torch.tensor(X, dtype=torch.float32).to('cuda')
    with torch.no_grad():
        return trained_embedder(X.to('cuda'))[:, 0].cpu()

X_train, y_train, X_test, y_test = ds.get_train_test_data(
    [dataset.DataAccessors.neighborhood_with_custom_embeddings(
        args.embedder, transform_embeddings, args.neighbors)],

     filters=[dataset.Helpers.filter_chains_with_valid_protrusion]
)

print ('fitting ... second model')

second_layer_model = BasicNetwork(
    batch_size=args.batch_size, 
    input_size=len(X_train[0]),
    hidden_sizes=args.hidden_layers_sec,
    epochs=args.epochs_sec,
    learning_rate=args.learning_rate_sec)

all_stats_in_second = []
second_losses = []

def second_epoch_callback(epoch, loss, model):
    if (args.verbose):
        print(f'Epoch: {epoch}, loss: {loss :.5f}', flush=True)

    second_losses.append(float(loss))

    if (epoch % args.epoch_stats_interval == 0):
        stats = get_statistics(model, 
                               X_test, y_test, epoch, print_res=args.verbose)
        all_stats_in_second.append(stats)

second_layer_model.register_epoch_callback(second_epoch_callback)

print ('fitting ... the second model', flush=True)
second_layer_model.fit(X_train, y_train)

final_stats = get_statistics(
    second_layer_model, 
    X_test, y_test, args.epochs_sec, print_res=args.verbose)

training_report = {
    'result_tag': args.tag,
    'ligand': args.ligand,
    'embedder': args.embedder,
    'model_to_string': f'first {model.underling_model()}; second {second_layer_model.underling_model()}',
    'batch_size': args.batch_size,
    'total_epochs': f'first: {args.epochs}, second: {args.epochs_sec}',
    'seed': args.seed,
    'learning_rate': f'first: {args.learning_rate}, second: {args.learning_rate}',
    'hidden_layers': args.hidden_layers + ['compose'] + args.hidden_layers_sec,
    'all_stats': all_stats_in_second,
    'all_stats_in_first': all_stats_in_first,
    'final_stats': final_stats,
    'losses': {
        'first': first_losses,
        'second': second_losses
    },
    'neighbors': args.neighbors
}

result_path = os.path.join(config.networks_results_folder, result_file_name)

with open(result_path, 'w') as file:
    json.dump(training_report, file, indent=2)

print ('report: ', training_report)

if (args.verbose):
    print('successfully finished', flush=True)