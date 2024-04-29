import argparse
from datetime import datetime
import json
import os
import random
import sys

import numpy as np
from sklearn.metrics import matthews_corrcoef
from sklearn.model_selection import train_test_split
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import config.config as config
from estimators.basic import BasicNetwork
from estimators.bypass import BypassedInputsNetwork
import data_prep.datasets_db as dataset
import data_prep.pdb_files_db as pdb_data
from evaluator import get_statistics
import string
from seed_network import seed_all 
from estimators.get_compressor import get_compressor_function


# srun -p gpu-short --gres=gpu:V100 -A nprg058s --cpus-per-task=4 --mem-per-cpu=32G python /home/brabecm4/diplomka/protein-binding-sites/netws/network.final_eval.py --final-tag basic_v6 --ligand AMP

allowed_ligands = dataset.SeqDatasetDb.all_ligands()
embedder = 'ESM'

parser = argparse.ArgumentParser(description='Ligand binding sites neural network')

parser.add_argument('--final-tag', type=str, help='Tag the final result', required=True)
parser.add_argument('--ligand', type=str, choices=allowed_ligands, help='Name of the ligand')

parser.add_argument('--second-radius', type=float,  help='Second protrsuion radius')

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

print ('defining model and data ...', flush=True)

tags_with_bypassed_input = ['protrusion_bypass_v4_c', 'SASA_bypassed_v2_c']

accessors = { 
    'basic_v6': [
        dataset.DataAccessors.embeddings(embedder)
    ],
    
    'prot_all_c': [
        dataset.DataAccessors.embeddings(embedder),
        dataset.DataAccessors.protrusion(1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0),
    ],
    
    'protrusion_bypass_v4_c': [
        dataset.DataAccessors.embeddings(embedder),
        dataset.DataAccessors.protrusion(hyper_params['radius']),
    ],
    'SASA_bypassed_v2_c': [
        dataset.DataAccessors.embeddings(embedder),
        dataset.DataAccessors.SASA_vector(),
    ],

    'neighboring_emb_5_v2_c': [
        dataset.DataAccessors.neighborhood_embeddings(embedder, 5),
    ],

    'one_prot_fst_v3_c': [
        dataset.DataAccessors.embeddings(embedder),
        dataset.DataAccessors.protrusion(hyper_params['radius']),
    ],

    'SASA_fst_v2_c': [
        dataset.DataAccessors.embeddings(embedder),
        dataset.DataAccessors.SASA_vector(),
    ],

    'nei_emb_3_v2_c': [
        dataset.DataAccessors.neighborhood_embeddings(embedder, 3),
    ],

    'nei_emb_5_avrg_v2_c': [
        dataset.DataAccessors.average_neighborhood_embeddings(embedder, 5),
    ],

    'nei_emb_3_avrg_v2_c': [
        dataset.DataAccessors.average_neighborhood_embeddings(embedder, 3),
    ],

    'two_prot_v2_c': [
        dataset.DataAccessors.embeddings(embedder),
        dataset.DataAccessors.protrusion(8.5),   
        dataset.DataAccessors.protrusion(args.second_radius),   
    ]
}[args.final_tag]


X_train_validate, y_train_validate, X_test, y_test = ds.get_train_test_data(
    accessors,
    filters=[dataset.Helpers.filter_chains_with_valid_protrusion]
)

def get_model():
    if args.final_tag not in tags_with_bypassed_input: 
        model = BasicNetwork(
            batch_size=hyper_params['batch_size'], 
            input_size=len(X_train_validate[0]),
            hidden_sizes=hyper_params['layers'],
            epochs=hyper_params['epochs'],
            learning_rate=hyper_params['learning_rate'])

    if args.final_tag in tags_with_bypassed_input: 
        model = BypassedInputsNetwork(
            batch_size=hyper_params['batch_size'], 
            input_size=len(X_train_validate[0]),
            hidden_sizes=hyper_params['layers'],
            epochs=hyper_params['epochs'],
            learning_rate=hyper_params['learning_rate'],
            bypassed_inputs=1)
            
    model.set_verbose()
    return model

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

all_X = np.concatenate((X_train_validate, X_test), axis=0)
all_y = np.concatenate((y_train_validate, y_test), axis=0)

_2cv_results = []

for iter in range(10):
    print('\n... running iteration ', iter , '\n')
    X1, X2, y1, y2 = split_data(all_X, all_y)

    for X_train_2cv, X_test_2cv, y_train_2cv, y_test_2cv in [[X1, X2, y1, y2], [X2, X1, y2, y1]]:
        print ('initializing network ...', flush=True)
        model = get_model()

        print ('fitting ...', flush=True)
        model.fit(X_train_2cv, y_train_2cv)

        print ('predicting ...', flush=True)
        pred_y = model.predict(X_test_2cv)

        fold_result = get_results(y_test_2cv, pred_y)
        _2cv_results.append(fold_result)

        print (f'   [mcc = {fold_result["mcc"]:.3f}]')

print ('running network for on entire train data...')

print ('initializing network ...', flush=True)
model = get_model()

print ('fitting')
model.fit(X_train_validate, y_train_validate)

print ('predicting')
y_pred = model.predict(X_test)

final_result = get_results(y_test, y_pred)

if args.second_radius:
    hyper_params['tag'] = hyper_params['tag'] + f'__{args.second_radius * 10:2.0f}'

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