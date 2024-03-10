import json
import os
import table_printer as printer
import argparse
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import config.config as config

parser = argparse.ArgumentParser(description='Show comparisons table')
parser.add_argument('--filters', type=str, nargs='+', help='List of filters')
parser.add_argument('--cols', type=str, nargs='+', help='List of additional columns to display')

# python /home/brabecm4/diplomka/protein-binding-sites/res_presentation/model_comparison_table.py --filters ligand==MN --cols radius

args = parser.parse_args()
filters = {
    f.split('==')[0]: f.split('==')[1]
    for f in args.filters or []
}
add_report_cols = args.cols or []

def load_from_dir(directory_path):
    files = os.listdir(directory_path)
    json_files = [file for file in files if file.endswith('.json')]

    json_objects = []

    for json_file in json_files:
        file_path = os.path.join(directory_path, json_file)
        with open(file_path, 'r') as f:
            json_object = json.load(f)
            json_objects.append(json_object)

    return json_objects

def filter_reports(reports):
    def to_normal_str(val):
        return f'{val}'.replace(" ", "").replace('"', '').replace("'", '')
    
    for key in filters:
        reports = [r for r in reports if to_normal_str(r[key])==filters[key]]

    return reports

def get_avrg_mcc(scores):
    MCCs = [score['mcc'] for score in scores]

    return f'{sum(MCCs) / len(MCCs):.3f} ∈ [{max(MCCs):.3f} - {min(MCCs):.3f}]'

def get_line(report):
    models = report['models']

    return [
        report['ligand'],
        report['tag'],
        report['p-value'],
        report['t-statistic']] + [
        
        report[col] if col in report else '---'
        for col in add_report_cols

        ] + [
        models,
        [
            get_avrg_mcc(report['scores'][models[0]]),
            get_avrg_mcc(report['scores'][models[1]]),
        ],
        'models are different' if report['p-value'] <= 0.05 else 'models are same'
    ]

reports = load_from_dir(config.model_comparisons_folder)
reports = filter_reports(reports)

header = ['ligand', 'tag', 'p-value', 't-stats'] + add_report_cols + [ 'models', 'avrg mcc', 'result']

table = [header] + [get_line(report) for report in reports]

printer.print_table(table)