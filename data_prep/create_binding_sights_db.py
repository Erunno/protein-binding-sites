import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import config.config as config
import re
import json

# expected folder structure
#  root directory
#  | - Training sets (this can have different name if adjusted below)
#  |   - <one training set - ligand prefixed (eg. AMP_Training.txt)>
#  |   - ...
#  | - Testing sets (this can have different name if adjusted below)
#  |   - <one testing set - ligand prefixed (eg. AMP_Testing.txt)>
#  |   - ...

yu_ds_folder = config.yu_sequences_folder
biding_sights_output_file = config.binding_sights_file

##### format parameters #####

training_folder_name = 'Training sets'
testing_folder_name = 'Testing sets'

protein_id_column = 0
binding_sights_column = 4

#############################

def get_all_file_names_in(folder):
    return [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]

def get_training_and_testing_folder_paths():
    global training_folder_name
    global testing_folder_name
    global yu_ds_folder

    training_folder_path = os.path.join(yu_ds_folder, training_folder_name)
    testing_folder_path = os.path.join(yu_ds_folder, testing_folder_name)

    return training_folder_path, testing_folder_path 

def get_all_ligands():
    training_folder_path, _ = get_training_and_testing_folder_paths()
    all_file_names = get_all_file_names_in(training_folder_path)

    all_ligands = [re.split(r'\W+|_', fname)[0] for fname in all_file_names]
    
    return all_ligands

def get_ligand_file_in(folder, ligand):
    all_files = get_all_file_names_in(folder)
    return [s for s in all_files if s.startswith(ligand)][0]    

def parse_binding_sights_to_list(raw_binding_sights_string):
    return [
        int(binding_sight_string[1:]) - 1 # in the file indexing starts from 1
        for binding_sight_string in raw_binding_sights_string.split(' ')
        if binding_sight_string != ''
    ]

def get_biding_sights_from(folder_path, ligand):
    global protein_id_column
    global binding_sights_column

    ligand_data_file_name = get_ligand_file_in(folder_path, ligand)
    ligand_data_file_path = os.path.join(folder_path, ligand_data_file_name)

    with open(ligand_data_file_path, 'r') as file:
        protein_records = [line.split(';') for line in file.readlines() if line != '']

    return {
        record[protein_id_column].lower(): parse_binding_sights_to_list(record[binding_sights_column])
        for record in protein_records
    }

def get_binding_sights_for(ligand):
    training_folder_path, testing_folder_path = get_training_and_testing_folder_paths()

    training_biding_sights = get_biding_sights_from(training_folder_path, ligand)
    testing_biding_sights = get_biding_sights_from(testing_folder_path, ligand)

    return {
        **training_biding_sights,
        **testing_biding_sights
    }

################################
#            main              #
################################

all_ligands = get_all_ligands()

all_biding_sights = {
    ligand: get_binding_sights_for(ligand) 
    for ligand in all_ligands
}



with open(biding_sights_output_file, 'w') as output_file:
    json.dump(all_biding_sights, output_file, indent=2)
