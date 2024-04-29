import json
import os
import statistics

from scipy import stats
import table_printer as printer
import argparse

# python /home/brabecm4/diplomka/protein-binding-sites/res_presentation/latex_table_printer.py --baseline-model-tag nei_emb_3_avrg_v2_c --results-folder /home/brabecm4/diplomka/protein-binding-sites/data/netw_results/final_runs --tags prot_all_c SASA_bypassed_v2_c
# python /home/brabecm4/diplomka/protein-binding-sites/res_presentation/latex_table_printer.py --baseline-model-tag basic_v6 --results-folder /home/brabecm4/diplomka/protein-binding-sites/data/netw_results/final_runs --tags prot_all_c SASA_bypassed_v2_c

all_ligands = ['ADP', 'AMP', 'ATP', 'CA', 'DNA', 'FE', 'GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN']
samples = 10

tag_mappings = {
    'basic_v6': '',
    'neighboring_emb_5_v2_c': '',
    'nei_emb_5_avrg_v2_c': '',
    'prot_all_c': '',
    'protrusion_bypass_v4_c': '',
    'SASA_bypassed_v2_c': '',
}

parser = argparse.ArgumentParser(description='Description of your script')
parser.add_argument('--baseline-model-tag', help='Tag of the model to compare to')
parser.add_argument('--tags', type=str, nargs='+', help='List of tags to compare')
parser.add_argument('--results-folder', help='File to save the hyperparameters')
args = parser.parse_args()

base_tag = args.baseline_model_tag

results = { ligand: {} for ligand in  all_ligands}
tags = set()

def get_mcc_s_of(ligand, tag, comparing_tag):
    all_mcc_s = [rec['mcc'] for rec in results[ligand][tag]['2cv_results']]
    all_mcc_s_comparing = [rec['mcc'] for rec in results[ligand][comparing_tag]['2cv_results']]
    return [mcc for mcc, comp_mcc in zip(all_mcc_s, all_mcc_s_comparing) if mcc > 0.01 and comp_mcc > 0.01][:samples]

def get_p_value(ligand, tag):
    if (tag not in results[ligand]) or (base_tag not in results[ligand]):
        return None

    base_dist = get_mcc_s_of(ligand, base_tag, comparing_tag=tag)
    tag_dist = get_mcc_s_of(ligand, tag, comparing_tag=base_tag)

    t_stat, p_value = stats.ttest_ind(base_dist, tag_dist)
    return p_value

def get_avrg(ligand, tag, comparing_tag):
    mcc_s = get_mcc_s_of(ligand, tag, comparing_tag)
    return sum(mcc_s) / len(mcc_s) 
    
def get_range(ligand, tag, comparing_tag):
    mcc_s = get_mcc_s_of(ligand, tag, comparing_tag)
    return f'[{min(*mcc_s):.3f}-{max(*mcc_s):.3f}]'

def get_stdv(ligand, tag, comparing_tag):
    mcc_s = get_mcc_s_of(ligand, tag, comparing_tag)
    return statistics.stdev(mcc_s)

def flatten(nested_list):
    return [item for sublist in nested_list for item in sublist]

def col(col, text):
    if isinstance(text, float):
        text = f'{text:.3f}'
    return '\\textcolor{' + col + '}{' + text + '}'

def bb(text):
    if isinstance(text, float):
        text = f'{text:.3f}'
    return '\\textbf{' + text + '}'

def ul(text):
    if isinstance(text, float):
        text = f'{text:.3f}'
    return '\\underline{' + text + '}'

def better_col(text, model_is_better):
    if model_is_better:
        text = ul(bb(text))

    return col('higher', text)

def worst_col(text):
    return col('lower', text)

def print_latex_header(tags):
    tags_count = len(tags)
    header = """
\\begin{table}[H]
    \\centering
    \\renewcommand{\\arraystretch}{\\tableLineHeightFactor}
    \\resizebox{1.0\\textwidth}{!}{
        \\begin{tabular}{@{}l""" + ('ccc' * tags_count) + """@{}}
        \\toprule
        \\textbf{Ligand} & \\multicolumn{2}{c}{\\textbf{Baseline}} & """ + \
\
        '&'.join([' \\multicolumn{3}{c}{\\textbf{' + tag.replace('_', '-') + '}} ' for i, tag in enumerate(tags[1:])]) + \
\
""" \\\\ 
        \\cmidrule(lr){2-3}""" + \
\
         ' '.join([' \\cmidrule(lr){' + f'{(i * 3) + 4}-{(i * 3) + 6}' + '}' for i in range(len(tags) - 1)]) + \
\
"""
        & MCC & Average MCC """ + (' & MCC & Average MCC & p-value ' * (tags_count - 1)) + """ \\\\
        \\midrule
"""
    print(header)

def print_latex_bottom():

    bottom = """
        \\bottomrule
    \\end{tabular}
    }
    \\caption{}
    \\label{tab:table_label}
\\end{table}
"""


    print(bottom)


def print_the_contents(tags):
    items_per_model = 3

    for ligand in all_ligands:
        lig_results = results[ligand]

        insignificant_line = 0 == len(['' for tag in tags if tag in results[ligand] and get_p_value(ligand, tag) < 0.05]) 

        b_col = 'muted' if insignificant_line else 'black'

        avrg = get_avrg(ligand, tags[0], comparing_tag=tags[0])
        stddev = get_stdv(ligand, tags[0], tags[0])
        avrg_entire_s = f'{avrg:.3f} (' + '\\pm ' + f'{stddev:.3f})'

        print(f'{col(b_col, ligand)} & {col(b_col, lig_results[tags[0]]["final_result"]["mcc"])} ' + \
              f'& {col(b_col, avrg_entire_s)} ', end='')

        for tag in tags[1:]:
            if tag not in results[ligand]:
                print(' & \\textcolor{muted}{?} ' * items_per_model, end='')
                continue


            record = results[ligand][tag]
            mcc = record['final_result']['mcc']
            base_mcc = results[ligand][base_tag]['final_result']['mcc']
            
            p_val = get_p_value(ligand, tag)
            
            avrg = get_avrg(ligand, tag, comparing_tag=base_tag)
            base_avrg = get_avrg(ligand, base_tag, comparing_tag=tag)

            stddev = get_stdv(ligand, tag, comparing_tag=base_tag)
            stddev_s = '(\\pm ' + f'{stddev:.3f})'
            
            insignificant_part = p_val > 0.05

            model_is_better = avrg > base_avrg and mcc > base_mcc


            if insignificant_part:
                mcc_s = col('muted', mcc)
                avrg_s = col('muted', avrg)
                p_val_s = col('muted', p_val)
                stddev_s = col('muted', stddev_s)
            else:
                mcc_s = better_col(mcc, model_is_better) if mcc > base_mcc else worst_col(mcc)
                avrg_s = better_col(avrg, model_is_better) if avrg > base_avrg else worst_col(avrg)
                p_val_s = f'{p_val:.3f}'

            print(f'& {mcc_s} & {avrg_s} {stddev_s} & ' + \
                  f'{p_val_s}', end='')

        print(' \\\\ \n')

for filename in os.listdir(args.results_folder):
    file_path = os.path.join(args.results_folder, filename)
    with open(file_path, 'r') as json_file:
        res = json.load(json_file)
        params = res['hyper_params']

        ligand = params['ligand']
        tag = params['tag']

        results[ligand][tag] = res
        tags.add(tag)

tags = [base_tag] + args.tags

print(tags)
print_latex_header(tags)
print_the_contents(tags)
print_latex_bottom()