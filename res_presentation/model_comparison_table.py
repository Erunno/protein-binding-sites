import json
import os
import table_printer as printer
import argparse
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import config.config as config


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

def get_avrg_mcc(scores):
    MCCs = [score['mcc'] for score in scores]
    return sum(MCCs) / len(MCCs)

def get_line(report):
    models = report['models']

    return [
        report['ligand'],
        report['p-value'],
        report['t-statistic'],
        models,
        [
            get_avrg_mcc(report['scores'][models[0]]),
            get_avrg_mcc(report['scores'][models[1]]),
        ],
        'models are different' if report['p-value'] <= 0.05 else 'models are same'
    ]

reports = load_from_dir(config.model_comparisons_folder)
header = ['ligand', 'p-value', 't-stats', 'models', 'avrg mcc', 'result']

table = [header] + [get_line(report) for report in reports]

printer.print_table(table)