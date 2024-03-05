import json
import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import config.config as config
from typing import Any, Callable, Dict, List, Union
import numpy as np
from numpy import ndarray
from Levenshtein import distance as lev_distance

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

    def ligand(self) -> str:
        return self.__data['ligand']

    def embeddings(self, embedder) -> Union[ndarray[ndarray[float]], None]:
        if self.__embedding_folder is None:
            return None
        
        if self.__embeddings_data is None:        
            emb_file = os.path.join(
                self.__embedding_folder, 
                embedder.upper(), 
                f'{self.full_id()}.npy')

            self.__embeddings_data = np.load(emb_file)

        return self.__embeddings_data
    
    def free_embeddings_from_RAM(self):
        self.__embeddings_data = None

    def full_id(self) -> str:
        return f'{self.protein_id()}{self.chain_id()}'

    def to_fasta_record(self) -> str:
        return f'>{self.full_id()}\n{self.sequence()}\n'
    
    def original_line(self) -> str:
        return self.__original_line
    
    def protrusion_vector_for(self, radius) -> Union[List[int], None]:
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

    def has_protrusion_record(self) -> bool:
        return self.__protrusion_record is not None
    
    def protrusion_sequence_matches_sequence(self) -> bool:
        dist_from_3d_sequence = self.get_lev_distance_of_3d_sequence_and_sequence()

        if dist_from_3d_sequence is None:
            return False
        
        return dist_from_3d_sequence == 0

    def get_lev_distance_of_3d_sequence_and_sequence(self) -> Union[int, None]:
        if self.__protrusion_record is None:
            return None

        _3d_structure_sequence = self.__protrusion_record['sequence']

        return lev_distance(_3d_structure_sequence, self.sequence())


    def __load_protrusion_record(self, protrusion_data):
        if protrusion_data is None:
            return None
        
        if self.protein_id() not in protrusion_data:
            return None
        
        protein_record = protrusion_data[self.protein_id()]
        chain_records = list(filter(lambda x: x['chain'].upper() == self.chain_id(),
                                    protein_record))
        
        if len(chain_records) == 0:
            return None
        
        if len(chain_records) == 1:
            return chain_records[0]
        
        return self.__get_best_matching_record(chain_records)

    def __get_best_matching_record(self, chain_records):
        best_rec = None
        best_dist = sys.maxsize

        for chain_rec in chain_records:
            _3d_structure_sequence = chain_rec['sequence']

            lev_dist = lev_distance(self.sequence(), _3d_structure_sequence)

            if lev_dist < best_dist:
                best_rec = chain_rec
                best_dist = lev_dist

            # prefer 'pdb' files over others
            if lev_dist == best_dist and chain_rec['origin_file'].endswith('pdb'):
                best_rec = chain_rec

        return best_rec
    
    def __assert_binding_residues(self, binding_residues: List[str], binding_indexes: List[int]):
        for i in range(len(binding_residues)):
            binding_idx = binding_indexes[i]
            binding_res = binding_residues[i]

            assert self.__data['sequence'][binding_idx] == binding_res

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

    def get_all_radii(self) -> List[float]:
        return Helpers.get_all_radii(self.all())
    
    def __construct_chain_record(self, record_line):
        return ChainRecord(record_line,
                           embedding_folder=self.__embedding_folder,
                           protrusion_data=self.__protrusion_data)
    
    def get_train_test_data(self, accessors, filters=[]):
        test, train = self.testing(), self.training()

        for filter in filters:
            train = filter(train)
            test = filter(test)

        X_train = Helpers.concat_chain_data(
            *accessors,
            chains=train
        )
        y_train = Helpers.concat_chain_data(
            DataAccessors.biding_sights_vect(),
            chains=train
        )

        X_test = Helpers.concat_chain_data(
            *accessors,
            chains=test
        )
        y_test = Helpers.concat_chain_data(
            DataAccessors.biding_sights_vect(),
            chains=test
        )

        return X_train, y_train, X_test, y_test

class SeqDatasetDb:
    @staticmethod
    def all_ligands() -> List[str]:
        return all_ligands
    
    def __init__(self, 
                 sequences_folder=config.yu_sequences_folder,
                 embeddings_folder=config.embeddings_folder):
         
        self.__sequences_folder = sequences_folder
        self.__protrusion_data = None
        self.__embedding_folder = embeddings_folder

    def get_dataset_for(self, ligand) -> LigandDataset:
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

        return list(dict.values())

    def __construct_ligand_ds(self, ligand):
        return LigandDataset(ligand,
                             sequences_folder=self.__sequences_folder,
                             embedding_folder=self.__embedding_folder,
                             protrusion_data=self.__protrusion_data)

class Helpers: 
    @staticmethod
    def filter_chains_with_protrusion(chains: List[ChainRecord]) -> List[ChainRecord]:
        return [chain for chain in chains if chain.has_protrusion_record()]

    @staticmethod
    def split_chains_to_valid_and_invalid_protrusion(chains: List[ChainRecord]) -> tuple[List[ChainRecord], List[ChainRecord]]:
        valid, invalid = [], []

        for chain in chains:
            array_to_extend = valid \
                              if chain.protrusion_sequence_matches_sequence() \
                              else invalid
            
            array_to_extend.append(chain) 

        return valid, invalid

    @staticmethod
    def filter_chains_with_valid_protrusion(chains: List[ChainRecord]) -> List[ChainRecord]:
        valid, _ = Helpers.split_chains_to_valid_and_invalid_protrusion(chains)
        return valid

    @staticmethod
    def concat_chain_data(*accessors, chains) -> ndarray:
        results = []

        for chain in chains:
            chain_vectors = accessors[0](chain)

            for accessor in accessors[1:]:
                part_vector = np.array(accessor(chain))
                chain_vectors = np.column_stack((chain_vectors, part_vector))

            for vect in chain_vectors:
                results.append(vect)
            
        return np.array(results)
        
    
    @staticmethod
    def get_all_radii(chains: List[ChainRecord]):
        all_radii = set(chains[0].get_all_protrusion_radii())

        for chain in chains:
            radii = set(chain.get_all_protrusion_radii())

            all_radii = all_radii.intersection(radii)

        return sorted(all_radii)

class DataAccessors:
    @staticmethod
    def biding_sights_vect():
        
        def get_biding_sights_vect(chain: ChainRecord):
            return chain.binding_sights()
        
        return get_biding_sights_vect
    
    @staticmethod
    def embeddings(embedder) -> Callable[[ChainRecord], Union[ndarray[ndarray[float]], None]]:
        
        def get_embeddings(chain: ChainRecord):
            return chain.embeddings(embedder)

        return get_embeddings

    @staticmethod
    def protrusion(*radii) -> Callable[[ChainRecord], Union[ndarray[ndarray[float]], None]]:
        
        def get_protrusion(chain: ChainRecord):
            result = np.array(chain.protrusion_vector_for(radii[0]))

            for radius in radii[1:]:
                protrusion = chain.protrusion_vector_for(radius)
                result = np.column_stack((result, protrusion))

            return result

        return get_protrusion