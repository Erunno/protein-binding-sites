import json
import math
import os
import re
import statistics
import sys
import curses
from scipy import stats
from contextlib import contextmanager

import table_printer as printer
import argparse
import culour

RESET = '\033[90m' # !!! not standard
# RESET = '\033[0m'
RED = '\033[91m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
MAGENTA = '\033[95m'

# RESET = ''
# RED = ''
# GREEN = ''
# YELLOW = ''
# BLUE = ''
# MAGENTA = ''

# python /home/brabecm4/diplomka/protein-binding-sites/res_presentation/compare_final.py --baseline-model-tag basic_v6 --results-folder /home/brabecm4/diplomka/protein-binding-sites/data/netw_results/final_runs 
# python /home/brabecm4/diplomka/protein-binding-sites/res_presentation/compare_final.py --baseline-model-tag one_prot_fst_v3_c --results-folder /home/brabecm4/diplomka/protein-binding-sites/data/netw_results/final_runs 

all_ligands = ['ADP', 'AMP', 'ATP', 'CA', 'DNA', 'FE', 'GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN']
samples = 10

tag_mappings = {
    'basic_v6': {
        'long': 'Baseline Model',
        'short': 'Baseline',
    },
    'neighboring_emb_5_v2_c': {
        'long': f'Big Network - 5 neighbors',
        'short': 'Big netw 5 nei.  ',
    },
    'nei_emb_5_avrg_v2_c': {
        'long': f'Small Network - 5 neighbors averaged',
        'short': 'Basic 5 nei. avrg',
    },
    'prot_all_c': {
        'long': 'All Protrusion Values in the First Layer',
        'short': 'All Protrusion   ',
    },
    'protrusion_bypass_v4_c': {
        'long': 'One Protrusion bypassing Hidden Layers',
        'short': 'Protrusion Bypass',
    },
    'SASA_bypassed_v2_c': {
        'long': 'One SASA value bypassing Hidden Layers',
        'short': 'SASA Bypass      ',
    },
    'one_prot_fst_v3_c': {
        'long': 'One Protrusion value in first layer',
        'short': 'One Prot in first',
    },
    'nei_emb_3_avrg_v2_c': {
        'long': 'Small Network - 3 neighbors averaged',
        'short': 'Basic 5 nei. avrg',
    },
    'nei_emb_3_v2_c': {
        'long': f'Big Network - 3 neighbors',
        'short': 'Big netw 3 nei.  ',
    },
    'SASA_fst_v2_c': {
        'long': 'SASA value in first layer',
        'short': 'SASA in first    ',
    },
}

parser = argparse.ArgumentParser(description='Description of your script')
parser.add_argument('--baseline-model-tag', help='Tag of the model to compare to')
parser.add_argument('--results-folder', help='File to save the hyperparameters')
parser.add_argument('--full-table', type=bool, default=False, help='Print full table results.')
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

def get_colored_cell(value, color):
    return f'{color}{value:.3f}{RESET}'

def get_colored_cell_if(value, threshold, below_is_good=True):
    if not value:
        return '---'
    
    is_green = value < threshold if below_is_good else value > threshold
    return get_colored_cell(value, GREEN if is_green else RED) 

def print_latex_table():
    pass

def print_the_table(tags, tag_aliases=None, to_compare_tag_index=None):
    if tag_aliases is None:
        tag_aliases = tag

    if to_compare_tag_index is not None:
        tags = [tags[0], tags[to_compare_tag_index]]
        tag_aliases = [tag_aliases[0], tag_aliases[to_compare_tag_index]]

    sub_cols =['mcc', 'avrg' if args.full_table else 'avrg          ', 'std dev', 'p-val']
    items_per_model = len(sub_cols)


    header = ['ligand'] + flatten([[tag] + [''] * (items_per_model - 1) for tag in tag_aliases ])
    subheader = [''] + sub_cols * len(tags)
    contents = []

    for ligand in all_ligands:
        line = [ligand]

        for tag in tags:
            if tag not in results[ligand]:
                line = line + ['---'] * items_per_model
                continue

            record = results[ligand][tag]
            mcc = record['final_result']['mcc']
            base_mcc = results[ligand][base_tag]['final_result']['mcc']

            p_val = get_p_value(ligand, tag)
            
            avrg = get_avrg(ligand, tag, comparing_tag=base_tag)
            base_avrg = get_avrg(ligand, base_tag, comparing_tag=tag)

            stddev = f'{BLUE}{get_stdv(ligand, tag, comparing_tag=base_tag):.3f}{RESET}'

            if base_mcc < mcc and p_val < 0.05:
                line.append(f'{YELLOW}{mcc:.3f} (+{mcc - base_mcc:.3f}){RESET}')
            else:
                line.append(get_colored_cell_if(mcc, threshold=base_mcc, below_is_good=False))

            if base_avrg < avrg and p_val < 0.05:
                line.append(f'{MAGENTA}{avrg:.3f} (+{avrg - base_avrg:.3f}){RESET}')
            else:
                line.append(get_colored_cell_if(avrg, threshold=base_avrg, below_is_good=False))

            line.append(stddev)
            line.append(get_colored_cell_if(p_val, 0.05))

        contents.append(line)
        
    printer.print_table([header, subheader, *contents])


def print_tags(menu_tags, highlighted_tag_index):
    print()
    for i, tag in enumerate(menu_tags):
        print(f'{YELLOW + "-> " if i == highlighted_tag_index else ""}{tag}{RESET}')
    print()

for filename in os.listdir(args.results_folder):
    file_path = os.path.join(args.results_folder, filename)
    with open(file_path, 'r') as json_file:
        res = json.load(json_file)
        params = res['hyper_params']

        ligand = params['ligand']
        tag = params['tag']

        results[ligand][tag] = res
        tags.add(tag)

tags = list(tags)
tags.sort(key=lambda x: x.lower() if x != base_tag else '\1')

tags_short = [tag_mappings[tag]['short'] if tag in tag_mappings else tag for tag in tags]
tags_long = [tag_mappings[tag]['long'] if tag in tag_mappings else tag for tag in tags]

print(tags_short)

if args.full_table:
    RESET = '\033[0m'
    print_the_table(tags)
    exit(0)


##### technical code displaying the interactive table ######

@contextmanager
def redirected_curses_output(window):
    original_stdout = sys.stdout
    
    class CursesWriter:
        def __init__(self, window):
            self.window = window
            self.y, self.x = window.getyx()

        def write(self, message):
            lines = message.split('\n')
            for i, line in enumerate(lines):
                if i > 0:
                    self.y += 1
                    self.x = 0
                try:
                    self.window.move(self.y, self.x)
                    parts = line.split('\033[')
                    for part in parts:
                        # Check if there is a color change request
                        if len(part) > 2 and 'm' == part[2]:  
                            code = part.split('m')[0]
                            try:
                                color_id = int(code)
                                self.window.attron(curses.color_pair(color_id))
                                text = 'm'.join(part.split('m')[1:])
                            except ValueError:
                                text = part
                        else:
                            text = part

                        try:
                            self.window.addstr(text)
                        except curses.error:
                            pass
                    
                except curses.error:
                    pass
                self.window.getyx()  # Update current cursor position

            self.y, self.x = window.getyx()

        def flush(self):
            pass  
    
    sys.stdout = CursesWriter(window)

    try:
        yield
    finally:
        sys.stdout = original_stdout

def setup_colors():
    curses.start_color()
    curses.init_pair(91, curses.COLOR_RED, curses.COLOR_BLACK)    # RED
    curses.init_pair(92, curses.COLOR_GREEN, curses.COLOR_BLACK)  # GREEN
    curses.init_pair(93, curses.COLOR_YELLOW, curses.COLOR_BLACK) # YELLOW
    curses.init_pair(94, curses.COLOR_BLUE, curses.COLOR_BLACK)   # BLUE
    curses.init_pair(95, curses.COLOR_MAGENTA, curses.COLOR_BLACK)# MAGENTA
    curses.init_pair(90, curses.COLOR_WHITE, curses.COLOR_BLACK)  # RESET (typically white or default)

def main(stdscr):
    curses.noecho()
    stdscr.keypad(True)
    setup_colors()

    current_tag_index = 1
    key = None

    while True:
        stdscr.clear()
        if key:
            offset = -1 if key == curses.KEY_UP else 1
            current_tag_index = (current_tag_index + offset + len(tags) - 1) % (len(tags) - 1)
        
        with redirected_curses_output(stdscr):
            print_tags(tags_long[1:], current_tag_index)  
            print_the_table(tags, tags_short, current_tag_index + 1)
            
        stdscr.refresh()
        key = stdscr.getch()

curses.wrapper(main)
