import json
import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import config.config as config
from typing import Dict, List
import numpy as np
from numpy import ndarray

all_ligands = ['ADP', 'AMP', 'ATP', 'CA', 'DNA', 'FE', 'GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN']

class ChainRecord:
    def __init__(self, csv_line: str, 
                 embedding_folder = None,
                 protrusion_data = None):
        cols = ['id', 'chain', '<col_2>', 'ligand', 'binding sights', 'sequence'] 

        self.__original_line = csv_line
        
        self.__data = {}
        line_parts = csv_line.split(';')

        for i in range(len(cols)):
            self.__data[cols[i]] = line_parts[i]

        self.__embedding_folder = embedding_folder
        self.__embeddings_data = None
        self.__protrusion_record = self.__load_protrusion_record(protrusion_data)

    def protein_id(self) -> str:
        return self.__data['id'].lower()

    def sequence(self) -> str:
        return self.__data['sequence']

    def chain_id(self) -> str:
        return self.__data['chain'].upper()

    def binding_sights(self) -> ndarray[int]:
        sights_str = self.__data['binding sights'].split(' ')

        binding_indexes = [int(sight[1:]) - 1 for sight in sights_str]
        binding_residues = [sight[0] for sight in sights_str]

        self.__assert_binding_residues(binding_residues, binding_indexes)

        binding_flags = [0] * len(self.__data['sequence'])

        for binding_idx in binding_indexes:
            binding_flags[binding_idx] = 1

        return np.array(binding_flags)

    def ligand(self):
        return self.__data['ligand']

    def embeddings(self, embedder):
        if self.__embedding_folder is None:
            return None
        
        if self.__embeddings_data is not None:
            return self.embeddings
        
        emb_file = os.path.join(
            self.__embedding_folder, 
            embedder.upper(), 
            f'{self.full_id()}.npy')

        self.__embeddings_data = np.load(emb_file)

        return self.__embeddings_data
    
    def free_embeddings_from_RAM(self):
        self.__embeddings_data = None

    def __assert_binding_residues(self, binding_residues: List[str], binding_indexes: List[int]):
        for i in range(len(binding_residues)):
            binding_idx = binding_indexes[i]
            binding_res = binding_residues[i]

            assert self.__data['sequence'][binding_idx] == binding_res

    def full_id(self):
        return f'{self.protein_id()}{self.chain_id()}'

    def to_fasta_record(self):
        return f'>{self.full_id()}\n{self.sequence()}\n'
    
    def original_line(self):
        return self.__original_line
    
    def protrusion_vector_for(self, radius) -> List[int]:
        if self.__protrusion_record is None:
            return None
        
        for radius_record in self.__protrusion_record['results']:
            if radius_record['radius'] == radius:
                return np.array(radius_record['protrusion'])
            
        return None
    
    def get_all_protrusion_radii(self) -> List[float]:
        if self.__protrusion_record is None:
            return []
        
        all_radii = [ radius_record['radius']
                      for radius_record 
                      in self.__protrusion_record['results'] ]

        return list(all_radii)

    def has_protrusion_record(self):
        return self.__protrusion_record is not None

    def __load_protrusion_record(self, protrusion_data):
        if protrusion_data is None:
            return None
        
        if self.protein_id() not in protrusion_data:
            return None
        
        protein_record = protrusion_data[self.protein_id()]
        chain_record = next(
            filter(lambda x: x['chain'].upper() == self.chain_id(),
                   protein_record), None)

        return chain_record

class ProteinRecord:
    def __init__(self):
        self.chains = {}

    def add_chain(self, chain: ChainRecord):
        self.chains[chain.chain_id()] = chain

    def all_chains(self) -> List[ChainRecord]:
        return list(self.chains)

    def chain(self, chain_id) -> ChainRecord:
        return self.chains[chain_id.upper()]

class LigandDataset:
    def __init__(self, ligand, 
                 sequences_folder, protrusion_data, embedding_folder):
        self.ligand = ligand.upper()
        
        self.__training_data = []
        self.__testing_data = []
        self.__protrusion_data = protrusion_data
        self.__embedding_folder = embedding_folder

        self.__load_all_data(sequences_folder)

    def training(self) -> List[ChainRecord]:
        return self.__training_data

    def testing(self) -> List[ChainRecord]:
        return self.__testing_data

    def all(self) -> List[ChainRecord]:
        return self.__testing_data + self.__training_data

    def __load_all_data(self, sequence_folder):
        for subfolder in ['Training sets', 'Testing sets']:
            folder_path = os.path.join(sequence_folder, subfolder)
            file_name = [file for file in os.listdir(folder_path) if file.startswith(self.ligand)][0]
            file_path = os.path.join(folder_path, file_name)

            with open(file_path, "r") as file:
                chain_records = [line.strip() for line in file if line.strip()]

            chain_records = [ self.__construct_chain_record(rec)
                              for rec in chain_records]

            if (subfolder == 'Training sets'):
                self.__training_data = chain_records
            else:
                self.__testing_data = chain_records

    def get_all_radii(self):
        all_radii = set(self.all()[0].get_all_protrusion_radii())

        for chain in self.all():
            radii = set(chain.get_all_protrusion_radii())

            all_radii = all_radii.intersection(radii)

        return sorted(all_radii)
    
    def __construct_chain_record(self, record_line):
        return ChainRecord(record_line,
                           embedding_folder=self.__embedding_folder,
                           protrusion_data=self.__protrusion_data)

class SeqDatasetDb:
    def __init__(self, 
                 sequences_folder=config.yu_sequences_folder,
                 embeddings_folder=config.embeddings_folder):
         
        self.__sequences_folder = sequences_folder
        self.__protrusion_data = None
        self.__embedding_folder = embeddings_folder

    def get_data_set_for(self, ligand):
        return self.__construct_ligand_ds(ligand)
    
    def load_protrusion_data_file(self, file_path):
        with open(file_path, 'r') as file:
            data = json.load(file)

        self.__protrusion_data = data

    def get_all_chain_records(self) -> List[ChainRecord]:
        chains: List[ChainRecord] = []  

        for ligand in all_ligands:
            ds = self.__construct_ligand_ds(ligand)
            chains = chains + ds.all()

        dict: Dict[str, ChainRecord]= {}

        for chain in chains:
            id = chain.full_id()

            if id in dict:
                # sanity check
                assert dict[id].sequence() == chain.sequence()

            dict[id] = chain

        return dict.values()

    def __construct_ligand_ds(self, ligand):
        return LigandDataset(ligand,
                             sequences_folder=self.__sequences_folder,
                             embedding_folder=self.__embedding_folder,
                             protrusion_data=self.__protrusion_data)

class Helpers: 
    @staticmethod
    def filter_chains_with_protrusion(chains: List[ChainRecord]) -> List[ChainRecord]:
        return [chain for chain in chains if chain.has_protrusion_record()]
