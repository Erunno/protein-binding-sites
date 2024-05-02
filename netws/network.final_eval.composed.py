import argparse
import json
import os
import random
import sys

import numpy as np
from sklearn.metrics import matthews_corrcoef
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import config.config as config
from estimators.basic import BasicNetwork
import data_prep.datasets_db as dataset
import data_prep.pdb_files_db as pdb_data
from seed_network import seed_all 
from estimators.get_compressor import get_compressor_function

tag_of_the_pretrained_model = 'basic_v6'

# srun -p gpu-short --gres=gpu:V100 -A nprg058s --cpus-per-task=4 --mem-per-cpu=32G python /home/brabecm4/diplomka/protein-binding-sites/netws/network.final_eval.composed.py --final-tag nei_3_comprs_v2_c --ligand AMP

print ('running final eval run...')

allowed_ligands = dataset.SeqDatasetDb.all_ligands()
embedder = 'ESM'

parser = argparse.ArgumentParser(description='Ligand binding sites neural network')

parser.add_argument('--final-tag', type=str, help='This tag\'s best hyperparameters are used', required=True)
parser.add_argument('--ligand', type=str, choices=allowed_ligands, help='Name of the ligand')

args = parser.parse_args()
print('args: ', args)

with open(config.best_HPs_file, 'r') as file:
    hyper_params = json.load(file)[args.final_tag][args.ligand]

seed_all(hyper_params['seed'])

print('tested hyperparameters: \n', hyper_params)

print ('loading data ...', flush=True)

pdb_db = pdb_data.PdbFilesDb()
db = dataset.SeqDatasetDb()
db.set_pdb_db(pdb_db)
ds = db.get_dataset_for(args.ligand)

print ('loading data ...', flush=True)

with open(config.best_HPs_file) as f:
    best_params_for_compressor = json.load(f)[tag_of_the_pretrained_model][args.ligand]

neighbors = {
    'nei_5_comprs_v3_cc': 5, 'nei_3_comprs_v3_cc': 3,
    'nei_compr_5_small_v4_c': 5, 'nei_compr_3_small_v4_c': 3}[args.final_tag]

X_train_validate, y_train_validate, X_test, y_test = ds.get_train_test_data(
    [dataset.DataAccessors.neighborhood_embeddings(embedder, neighbors)],
    filters=[dataset.Helpers.filter_chains_with_valid_3D_file]
)

def get_model(input_size):
    new_model = BasicNetwork(
        batch_size=hyper_params['batch_size'], 
        input_size=input_size,
        hidden_sizes=hyper_params['layers'],
        epochs=hyper_params['epochs'],
        learning_rate=hyper_params['learning_rate'])
    new_model.set_verbose()
    return new_model

split_random = random.Random()
split_random.seed(42)

def split_data(X, y): 
    combined = list(zip(X, y))
    split_random.shuffle(combined)
    X_s, y_s = zip(*combined)
    
    half = len(X_s) // 2

    return \
        np.array(X_s[:half]), np.array(X_s[half:]), \
        np.array(y_s[:half]), np.array(y_s[half:])

def get_split_mask(count):
    items = range(count)
    split_random.shuffle(items)

    flags = [0] * count
    half = count // 2

    for i in items[half:]:
        flags[i] = 1

    return np.array(flags)

def slice_matrix_column_wise(matrix):
    num_cols = matrix.shape[1] // neighbors
    matrix_parts = []
        
    for i in range(neighbors):
        part = matrix[:, num_cols * i:num_cols * (i + 1)]
        matrix_parts.append(part)

    return matrix_parts


def get_results(y_true, y_pred):
    y_true, y_pred = np.array(y_true), np.array(y_pred)
    y_pred_class = (y_pred >= hyper_params['threshold']).astype(int)
    
    tp = np.sum((y_pred_class == 1) & (y_true == 1)).item()
    fp = np.sum((y_pred_class == 1) & (y_true == 0)).item()
    tn = np.sum((y_pred_class == 0) & (y_true == 0)).item()
    fn = np.sum((y_pred_class == 0) & (y_true == 1)).item()
    confusion_matrix = {'tp': tp, 'fp': fp, 'tn': tn, 'fn': fn}

    mcc = matthews_corrcoef(y_true, y_pred_class)
    
    return {
        'mcc': mcc,
        'confusion_matrix': confusion_matrix,
    }

print('running 10 2-cross validations...')

_2cv_results = []

all_X = np.concatenate((X_train_validate, X_test), axis=0)
all_y = np.concatenate((y_train_validate, y_test), axis=0)

for iter in range(10):

    print('\n... running iteration ', iter , '\n')

    X1, X2, y1, y2 = split_data(all_X, all_y)

    for X_train_2cv, X_test_2cv, y_train_2cv, y_test_2cv in [[X1, X2, y1, y2], [X2, X1, y2, y1]]:
       
        print ('orig', X_train_2cv.shape, '\n')
        print ('training compressor ...', flush=True)

        train_parts = slice_matrix_column_wise(X_train_2cv)
        X_train_2cv = None

        compressor = get_compressor_function(
            train_parts[0], y_train_2cv,
            best_params_for_compressor)
        
        print ('compressed', compressor(train_parts[0]).shape, '\n')

        print ('compressing vectors...')
        test_parts = slice_matrix_column_wise(X_test_2cv)
        X_test_2cv = np.hstack([compressor(part) for part in test_parts])

        X_train_2cv = np.hstack([compressor(part) for part in train_parts])

        print ('new', X_train_2cv.shape, '\n')

        print ('initializing second network ...', flush=True)
        model = get_model(input_size=len(X_train_2cv[0]))

        print ('fitting second network ...', flush=True)
        model.fit(X_train_2cv, y_train_2cv)

        print ('predicting ...', flush=True)
        pred_y = model.predict(X_test_2cv)

        fold_result = get_results(y_test_2cv, pred_y)
        _2cv_results.append(fold_result)

        print (f'   [mcc = {fold_result["mcc"]:.3f}]')

print ('running network for on entire train data...')

print ('compressing ...')
train_parts = slice_matrix_column_wise(X_train_validate)
compressor = get_compressor_function(
            train_parts[0], y_train_validate,
            best_params_for_compressor)

X_train_validate = np.hstack([compressor(part) for part in train_parts])

test_parts = slice_matrix_column_wise(X_test)
X_test = np.hstack([compressor(part) for part in test_parts])

print ('initializing second network ...', flush=True)
model = get_model(input_size=len(X_train_validate[0]))

print ('fitting')
model.fit(X_train_validate, y_train_validate)

print ('predicting')
y_pred = model.predict(X_test)

final_result = get_results(y_test, y_pred)

report = {
    '2cv_results': _2cv_results,
    'final_result': final_result,
    'hyper_params': hyper_params
}

result_file_name = f'{hyper_params["tag"]}_{hyper_params["ligand"]}.report.json'
result_path = os.path.join(config.final_networks_results_folder, result_file_name)

with open(result_path, 'w') as file:
    json.dump(report, file, indent=2)

print ('report: ', report)

print('successfully finished', flush=True)