import argparse
from datetime import datetime
import json
import os
import random
import sys
import torch

from sklearn.model_selection import train_test_split

from estimators.get_compressor import get_compressor_function
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

tag_of_the_pretrained_model = 'basic_v6'

# srun -p gpu-short --gres=gpu:V100 -A nprg058s --cpus-per-task=4 --mem-per-cpu=32G python /home/brabecm4/diplomka/protein-binding-sites/netws/network.composed.v2.py --verbose True --seed 42 --learning-rate 0.01 --epochs 10 --batch-size 1000 --embedder ESM --epoch-stats-interval 10 --tag test --hidden-layers 256 256 256 32 --ligand FE --neighbors 3

allowed_ligands = dataset.SeqDatasetDb.all_ligands()
# allowed_embedder = ['BERT', 'ESM', 'T5']
allowed_embedder = ['ESM']

parser = argparse.ArgumentParser(description='Ligand binding sites neural network')

parser.add_argument('--tag', type=str, help='Tag the final result', required=True)
parser.add_argument('--ligand', type=str, choices=allowed_ligands, help='Name of the ligand')
parser.add_argument('--embedder', type=str, choices=allowed_embedder, help='Embedder to be used.')
parser.add_argument('--hidden-layers', type=int, nargs='+', help='List of hidden layer sizes')
parser.add_argument('--batch-size', type=int, help='Size of one batch')
parser.add_argument('--epochs', type=int, help='Number of epochs')
parser.add_argument('--epoch-stats-interval', type=int, help='Interval of gathering statistics.')
parser.add_argument('--seed', type=int, help='Seed of random.')
parser.add_argument('--neighbors', type=int, help='Seed of random.')
parser.add_argument('--learning-rate', type=float, help='Learning rate')
parser.add_argument('--verbose', type=bool, default=False, help='Print intermediate results.')

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

print ('training compressors...')
print ('loading data for compressor...')

X_train_validate, y_train_validate, _, _ = ds.get_train_test_data(
    [
        dataset.DataAccessors.embeddings(args.embedder),
    ],
    filters=[dataset.Helpers.filter_chains_with_valid_3D_file]
)

X_train, _, y_train, _ = train_test_split(
    X_train_validate, y_train_validate, train_size=0.8, random_state=args.seed)

with open(config.best_HPs_file) as f:
	best_params_for_compressor = json.load(f)[tag_of_the_pretrained_model][args.ligand]

print ('defining compressor with params: ', best_params_for_compressor)

compressor = get_compressor_function(
    X_train, y_train, 
    best_params_for_compressor)

print('compressor finished...')

print ('defining training and testing data for the network ...', flush=True)

X_train_validate, y_train_validate, X_test, y_test = ds.get_train_test_data(
    [
        dataset.DataAccessors.neighborhood_with_custom_embeddings(
            args.embedder, compressor, args.neighbors),
    ],

     filters=[dataset.Helpers.filter_chains_with_valid_3D_file]
)

X_train, X_validate, y_train, y_validate = train_test_split(
    X_train_validate, y_train_validate, train_size=0.8, random_state=args.seed)

print()

print ('initializing network ...', flush=True)

model = BasicNetwork(batch_size=args.batch_size, 
                     input_size=len(X_train[0]),
                     hidden_sizes=args.hidden_layers,
                     epochs=args.epochs,
                     learning_rate=args.learning_rate)

all_stats = []
losses = []

def get_stats(model, epoch):
    if (args.verbose):
        print ('validation data stats:')

    stats = get_statistics(model, 
                            X_validate, y_validate, 
                            epoch, print_res=args.verbose)
    
    if (args.verbose):
        print ('test data stats:')
    
    stats['test_data_stats'] = get_statistics(model, 
                               X_test, y_test, 
                               epoch, print_res=args.verbose)
    
    return stats

def epoch_callback(epoch, loss, model):
    if (args.verbose):
        print(f'Epoch: {epoch}, loss: {loss :.5f}', flush=True)

    losses.append(float(loss))

    if (epoch % args.epoch_stats_interval == 0):
        model.optimize_threshold_for_mcc(X_train, y_train)

        all_stats.append(get_stats(model, epoch))

model.register_epoch_callback(epoch_callback)

print (f'fitting ... (data shape: {X_train.shape})', flush=True)
model.fit(X_train, y_train)

print ('predicting ...', flush=True)
model.predict(X_test)

final_stats = get_stats(model, args.epochs)

training_report = {
    'result_tag': args.tag,
    'ligand': args.ligand,
    'embedder': args.embedder,
    'model_to_string': str(model.underling_model()),
    'batch_size': args.batch_size,
    'total_epochs': args.epochs,
    'seed': args.seed,
    'learning_rate': args.learning_rate,
    'hidden_layers': args.hidden_layers,
    'all_stats': all_stats,
    'final_stats': final_stats,
    'losses': losses,
    'params_of_compressor': best_params_for_compressor,
    'neighbors': args.neighbors
}

result_path = os.path.join(config.networks_results_folder, result_file_name)

with open(result_path, 'w') as file:
    json.dump(training_report, file, indent=2)

print ('report: ', training_report)

if (args.verbose):
    print('successfully finished', flush=True)