import argparse
import math
import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import config.config as config
import numpy as np
import random
import string
import torch
from seed_network import seed_all 
from datetime import datetime
from estimators.basic import BasicNetwork
from estimators.bypass import BypassedInputsNetwork
from estimators.partial import PartialInputNetwork
from mlxtend.evaluate import paired_ttest_5x2cv
from sklearn.metrics import matthews_corrcoef
import data_prep.datasets_db as datasets
import pprint

parser = argparse.ArgumentParser(description='Ligand binding sites model comparing')
parser.add_argument('--ligand', type=str, help='Name of the ligand')
parser.add_argument('--radius', type=float, help='Protrusion radius')
args = parser.parse_args()

protrusion_fname = '/home/brabecm4/diplomka/protein-binding-sites/data/3d_proc/protrusion.max-neighbors.big.json'

tag = "v2_radii_tests"
ligand = args.ligand
batch_size = 1024
epochs = 80
learning_rate = 0.001
hidden_layers = [512, 256, 128]
embedder = "ESM"
radius=args.radius

seed_all(42)

data_loader = None 
radii_count = 0

print("loading data ...", flush=True)

db = datasets.SeqDatasetDb()
db.load_protrusion_data_file(protrusion_fname)

ds = db.get_dataset_for(ligand)

X_train, y_train, X_test, y_test = ds.get_train_test_data(
    [datasets.DataAccessors.embeddings(embedder), datasets.DataAccessors.protrusion(radius)],
    [datasets.Helpers.filter_chains_with_valid_protrusion]
)

X_train = torch.tensor(X_train, dtype=torch.float32)
y_train = torch.tensor(y_train, dtype=torch.float32)

X_test = torch.tensor(X_test, dtype=torch.float32)
y_test = torch.tensor(y_test, dtype=torch.float32)

X = torch.cat((X_train, X_test), axis=0)
y = torch.cat((y_train, y_test), axis=0)

print("defining estimators ...", flush=True)

input_size = len(X_train[0])

model_1 = PartialInputNetwork(
    input_size=input_size - 1,
    model = BasicNetwork(
        input_size=input_size - 1,
        epochs=epochs,
        batch_size=batch_size,
        hidden_sizes=hidden_layers,
        learning_rate=learning_rate
    ) 
) 
model_1.set_name("embeddings_only")
model_1.set_verbose()


model_2 = BypassedInputsNetwork(
    input_size=input_size,
    epochs=epochs,
    batch_size=batch_size,
    hidden_sizes=hidden_layers,
    learning_rate=learning_rate,
    bypassed_inputs=1
)
model_2.set_name("protrusion_in_last_layer")
model_2.set_verbose()

scores = {}
scores[model_1.get_name()] = []
scores[model_2.get_name()] = []

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

    result = scores[model.get_name()]
    result.append({
        'mcc': mcc,
        'counts': {
            'reference': {
                'zeros': np.count_nonzero(y == 0), 
                'ones': np.count_nonzero(y == 1) 
            },
            'model': {
                'zeros': np.count_nonzero(pred_y == 0), 
                'ones': np.count_nonzero(pred_y == 1) 
            },
        }
    })

    return mcc 

print("evaluating...")

t, p = paired_ttest_5x2cv(
    estimator1=model_1, estimator2=model_2, 
    X=X, y=y.int(), scoring=mcc_scorer, random_seed=1)

print('P-value: %.3f, t-Statistic: %.3f' % (p, t))

if p <= 0.05:
    print('Difference between mean performance is probably real')
else:
    print('Algorithms probably have the same performance')

report = {
    'tag': tag,
    'models': [model_1.get_name(), model_2.get_name()],
    'scores': scores,
    'ligand': ligand,
    'p-value': p if not math.isnan(p) else 'nan',
    't-statistic': t if not math.isnan(t) else 'nan',
    'batch_size': batch_size,
    'epochs': epochs,
    'learning_rate': learning_rate,
    'hidden_layers': hidden_layers,
    'embedder': embedder,
    'radius': radius,
}

str_report = pprint.pformat(report, compact=True, width=100).replace("'",'"')

def generate_random_string(length):
    letters = string.ascii_letters + string.digits
    return ''.join(random.choice(letters) for _ in range(length))

random_suffix = generate_random_string(5)
result_file_name = f'{ligand}_hl{"-".join([str(n) for n in hidden_layers])}.{datetime.now().strftime("%d-%m-%y-%Hh%Mm%Ss.%f")}.{random_suffix}.json'

result_path = os.path.join(config.model_comparisons_folder, result_file_name)
with open(result_path, 'w') as file:
    file.write(str_report)

print ('report: ', str_report)

print('successfully finished', flush=True)