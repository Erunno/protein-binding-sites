import json
import os
import sys 
import numpy as np

class DataLoader:

    def __init__(self, binding_sights_db_fname, dataset_by_ligands_db_fname, embeddings_folder, verbose=False):
        self.embeddings_folder = embeddings_folder

        with open(binding_sights_db_fname, 'r') as json_file:
            self.binding_sights_db = json.load(json_file)

        with open(dataset_by_ligands_db_fname, 'r') as json_file:
            self.datasets_by_ligands = json.load(json_file)

        if verbose:
            print('info: checking presence of all embeddings required by all datasets')

        for ligand in self.datasets_by_ligands:
            self.datasets_by_ligands[ligand] = self._assert_and_filter_not_found_proteins(self.datasets_by_ligands, ligand)

        if verbose:
            print('info: check completed')

    def get_data_set_for(self, ligand: str):
        dataset = self.datasets_by_ligands[ligand.upper()]
        all_proteins_of_ds = dataset['test'] + dataset['train']

        embeddings_for_proteins = {
            protein_id: self.load_protein_embeddings(protein_id) 
            for protein_id in all_proteins_of_ds
        }

        def get_Xy_for(dataType):
            X = np.array([
                one_embedding 
                for protein_id in dataset[dataType]
                for one_embedding in embeddings_for_proteins[protein_id]
            ])

            y = np.array([
                one_y_value

                for protein_id in dataset[dataType]
                for one_y_value in self._get_y_vector(
                    self.binding_sights_db[ligand][protein_id], 
                    len(embeddings_for_proteins[protein_id])
                )
            ])

            return X, y
        
        X_train, y_train = get_Xy_for('train')
        X_test,  y_test  = get_Xy_for('test')

        return X_train, y_train, X_test, y_test

    def _get_y_vector(self, binding_indexes, vector_len):
        return np.array([
            1 if idx in binding_indexes else 0
            for idx in range(vector_len)
        ])

    def _get_protein_embedding_fname(self, protein_id):
        return os.path.join(self.embeddings_folder, f'{protein_id.lower()}.npy') 

    def load_protein_embeddings(self, protein_id: str):
        protein_fname = self._get_protein_embedding_fname(protein_id)
        try:
            return np.load(protein_fname)      
        except:
            return None
        
    def _assert_and_filter_not_found_proteins(self, all_datasets, ligand):
        dataset = all_datasets[ligand]
        not_found = []
        
        for protein in dataset['train'] + dataset['test']:
            protein_fname = self._get_protein_embedding_fname(protein)
            embeddings_exists = os.path.exists(protein_fname)
            
            if not embeddings_exists:
                not_found.append(protein)

        if len(not_found) != 0:
            print(f'Err for ligand "{ligand}": following proteins not found in the folder: {", ".join(not_found)}\n\n', file=sys.stderr)

        return {
            dataset_type: [protein for protein in dataset[dataset_type] if protein not in not_found]
            for dataset_type in ['train', 'test']
        }
