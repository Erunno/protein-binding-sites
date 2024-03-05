import argparse
from datetime import datetime
import json
import os
import random
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import config.config as config
from estimators.basic import BasicNetwork
import data_prep.datasets_db as dataset
from evaluator import get_statistics
import string
from seed_network import seed_all 

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
parser.add_argument('--learning-rate', type=float, help='Learning rate')
parser.add_argument('--verbose', type=bool, default=False, help='Print intermediate results.')
# parser.add_argument('--use-simple-model', type=bool, default=False, help='Use simple sequential model.')
# parser.add_argument('--protrusion-data-file', type=str, help='Path to protrusion data')
# parser.add_argument('--pdb-mappings-fname', type=str, help='Path to mappings to pdb files')
# parser.add_argument('--used-protrusion-radii', type=float, nargs='+', help='List of used radii for protrusion. (only those found in "--protrusion-data-file" can be specified)')

args = parser.parse_args()
args.verbose = True if args.verbose else False

seed_all(args.seed)

def generate_random_string(length):
    letters = string.ascii_letters + string.digits
    return ''.join(random.choice(letters) for _ in range(length))

random_suffix = generate_random_string(5)
result_file_name = f'{args.ligand}_hl{"-".join([str(n) for n in args.hidden_layers])}.{datetime.now().strftime("%d-%m-%y-%Hh%Mm%Ss.%f")}.{random_suffix}.json'


print ('loading data ...', flush=True)

db = dataset.SeqDatasetDb()
ds = db.get_dataset_for(args.ligand)

print ('defining training and testing data ...', flush=True)

X_train, y_train, X_test, y_test = ds.get_train_test_data(
    [dataset.DataAccessors.embeddings(args.embedder)]
)

print ('initializing network ...', flush=True)

model = BasicNetwork(batch_size=args.batch_size, 
                     input_size=len(X_train[0]),
                     hidden_sizes=args.hidden_layers,
                     epochs=args.epochs,
                     learning_rate=args.learning_rate)

all_stats = []

def epoch_callback(epoch, loss, model):
    if (args.verbose):
        print(f'Epoch: {epoch}, loss: {loss :.5f}', flush=True)

    if (epoch % args.epoch_stats_interval == 0):
        stats = get_statistics(model, 
                               X_test, y_test, epoch, print_res=args.verbose)
        all_stats.append(stats)

model.register_epoch_callback(epoch_callback)

print ('fitting ...', flush=True)
model.fit(X_train, y_train)

print ('predicting ...', flush=True)
model.predict(X_test)

final_stats = get_statistics(
    model, 
    X_test, y_test, args.epochs, print_res=args.verbose)

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
    'final_stats': final_stats
}

result_path = os.path.join(config.networks_results_folder, result_file_name)

with open(result_path, 'w') as file:
    json.dump(training_report, file, indent=2)

print ('report: ', training_report)

if (args.verbose):
    print('successfully finished', flush=True)