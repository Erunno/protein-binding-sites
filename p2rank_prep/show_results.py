import os
import re

# python /home/brabecm4/diplomka/protein-binding-sites/p2rank_prep/show_results.py

all_ligands = ['ADP', 'AMP', 'ATP', 'CA', 'DNA', 'FE', 'GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN']

results_folder = '/home/brabecm4/diplomka/protein-binding-sites/data/netw_results/netw_runs/'

for ligand in all_ligands:
    print(f'{ligand:5} ', end='')

    res_log_fname = f'/home/brabecm4/diplomka/protein-binding-sites/p2rank/p2rank-llmex-data/OUTPUT/_{ligand}/traineval_{ligand}_train_clean_{ligand}_test_clean/run.log'

    try:
        with open(res_log_fname, 'r') as file:
            content = file.read().split('\n')
    except:
        print('... no log file')
        continue

    res = {}
    for metric in ['MCC', 'ACC']:
        relevant_lines = [line for line in content if line.startswith(metric)]

        if len(relevant_lines) == 0:
            print('... no MCC values in log file')
            continue

        vals = []

        for line in relevant_lines:
            match = re.search(rf'{metric}:\s*([0-9.]+)', line)
            value = float(match.group(1))

            vals.append(value)

        vals = sorted(set(vals), reverse=True)
        
        res[metric] = vals[0]

        print (' ', metric, '->', vals, end='')

    print()

    res_object = '''
{
  "ligand": "{ligand}",
  "embedder": "p2rank",
  "result_tag": "p2rank",
  "model_to_string": "?",
  "batch_size": "?",
  "total_epochs": "?",
  "seed": "?",
  "learning_rate": "?",
  "hidden_layers": [
    "?"
  ],
  "all_stats": [
    {
      "epochs": "?",
      "acc": {acc},
      "mcc": {mcc},
      "test_data_stats": {
        "epochs": "?",
        "acc": {acc},
        "mcc": {mcc}
      }
    }
  ],
  "final_stats": {
    "epochs": "?",
    "acc": {acc},
    "mcc": {mcc}
  }
}
    '''

    res_object = res_object \
                    .replace('{ligand}', str(ligand)) \
                    .replace('{acc}', str(res['ACC'])) \
                    .replace('{mcc}', str(res['MCC']))
    
    with open(os.path.join(results_folder, f'_p2rank_{ligand}.json'), 'w') as f:
        f.write(res_object)
