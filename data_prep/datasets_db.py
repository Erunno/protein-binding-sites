import json
import os
import sys

from data_prep.file_cache import use_cache
from data_prep.pdb_files_db import Chain3dStructure, PdbFilesDb
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import config.config as config
from typing import Any, Callable, Dict, List, Union
import numpy as np
from numpy import ndarray
from Levenshtein import distance as lev_distance

all_ligands = ['ADP', 'AMP', 'ATP', 'CA', 'DNA', 'FE', 'GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN']

class ChainRecord:
    pass

class ChainRecord:
    def __init__(self, csv_line: str, 
                 embedding_folder = None,
                 pdb_db: PdbFilesDb = None):
        
        cols = ['id', 'chain', 'binding_sight_ID', 'ligand', 'binding_sights', 'sequence'] 

        self.__original_line = csv_line
        
        self.__data = {}
        line_parts = csv_line.split(';')

        for i in range(len(cols)):
            self.__data[cols[i]] = line_parts[i]

        self.__embedding_folder = embedding_folder
        self.__embeddings_data = None

        self.__pdb_db = pdb_db
        self.__chain_structure = None

        self.__binding_sights = [self.__data['binding_sights']]

    def protein_id(self) -> str:
        return self.__data['id'].lower()

    def sequence(self) -> str:
        return self.__data['sequence']

    def chain_id(self) -> str:
        return self.__data['chain'].upper()
    
    def add_another_binding_sight_from(self, chain: ChainRecord):
        self.__binding_sights.append(chain.__data['binding_sights'])

    def binding_sights(self) -> ndarray[int]:
        return self.__compute_biding_sights(' '.join(self.__binding_sights))
    
    def original_binding_sights(self):
        return self.__compute_biding_sights(self.__data['binding_sights'])

    def __compute_biding_sights(self, bindings_string: str):
        sights_str = bindings_string.split(' ')

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
        with use_cache(self.get_chain_structure()) as chain_structure:
            return np.array(chain_structure.get_protrusion_vector(radius=radius))

    def protrusion_sequence_matches_sequence(self) -> bool:
        dist_from_3d_sequence = self.get_lev_distance_of_3d_sequence_and_sequence()

        if dist_from_3d_sequence is None:
            return False
        
        return dist_from_3d_sequence == 0

    def get_lev_distance_of_3d_sequence_and_sequence(self) -> Union[int, None]:
        if self.__pdb_db is not None:
            with use_cache(self.get_chain_structure()) as chain_structure:
                try:
                    _3d_structure_sequence = chain_structure.compute_sequence()
                except:
                    raise Exception(f'Protein {self.full_id()} does not have a corresponding 3D structure or the sequence cannot be read from the PDB file.')

        return lev_distance(_3d_structure_sequence, self.sequence())

    def get_chain_structure(self, loaded=False) -> Chain3dStructure:
        if self.__chain_structure is None:
            self.__chain_structure = self.__pdb_db.get_chain_structure(self.protein_id(), self.chain_id())

        if loaded:
            self.__chain_structure.load(preferred_sequence=self.sequence())

        return self.__chain_structure
    
    def get_SASA_vector(self):
        with use_cache(self.get_chain_structure()) as chain_structure:
            return chain_structure.get_SASA_vector()
  
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
                 sequences_folder, embedding_folder, pdb_db):
        self.ligand = ligand.upper()
        
        self.__training_data = []
        self.__testing_data = []

        self.__training_data_per_binding_sight = []
        self.__testing_data_per_binding_sight = []
        
        self.__embedding_folder = embedding_folder
        self.__pdb_db = pdb_db

        self.__load_all_data(sequences_folder)

    def training(self) -> List[ChainRecord]:
        return self.__training_data

    def testing(self) -> List[ChainRecord]:
        return self.__testing_data

    def training_per_binding_sight(self) -> List[ChainRecord]:
        return self.__training_data_per_binding_sight

    def testing_per_binding_sight(self) -> List[ChainRecord]:
        return self.__testing_data_per_binding_sight

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

            seen_chains = {}

            for chain in chain_records:
                id = chain.full_id()

                if id in seen_chains:
                    seen_chains[id].add_another_binding_sight_from(chain)
                else:
                    seen_chains[id] = chain

            if (subfolder == 'Training sets'):
                self.__training_data = list(seen_chains.values())
                self.__training_data_per_binding_sight = chain_records
            else:
                self.__testing_data = list(seen_chains.values())
                self.__testing_data_per_binding_sight = chain_records

        self.__training_data.sort(key=lambda ch: ch.full_id())
        self.__testing_data.sort(key=lambda ch: ch.full_id())

    def get_all_radii(self) -> List[float]:
        return Helpers.get_all_radii(self.all())
    
    def __construct_chain_record(self, record_line):
        return ChainRecord(record_line,
                           embedding_folder=self.__embedding_folder,
                           pdb_db=self.__pdb_db)
    
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
        self.__embedding_folder = embeddings_folder
        self.__pdb_db = None

    def get_dataset_for(self, ligand) -> LigandDataset:
        return self.__construct_ligand_ds(ligand)
    
    def set_pdb_db(self, pdb_db: PdbFilesDb):
        self.__pdb_db = pdb_db
    
    def get_all_chain_records_with_merged_binding_sites(self):
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

                # merge the binding sites of both chains
                chain.add_another_binding_sight_from(dict[id])

            dict[id] = chain

        return list(dict.values())

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
                             pdb_db=self.__pdb_db)

class Helpers: 
    @staticmethod
    def split_chains_to_valid_and_invalid_3D_file(chains: List[ChainRecord]) -> tuple[List[ChainRecord], List[ChainRecord]]:
        valid, invalid = [], []

        for chain in chains:
            array_to_extend = valid \
                              if chain.protrusion_sequence_matches_sequence() \
                              else invalid
            
            array_to_extend.append(chain) 

        return valid, invalid

    @staticmethod
    def filter_chains_with_valid_3D_file(chains: List[ChainRecord]) -> List[ChainRecord]:
        valid, _ = Helpers.split_chains_to_valid_and_invalid_3D_file(chains)
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
    
    @staticmethod
    def SASA_vector() -> Callable[[ChainRecord], Union[ndarray[ndarray[float]], None]]:
        def get_SASA_vector(chain: ChainRecord):
            return np.array(chain.get_SASA_vector())

        return get_SASA_vector

    @staticmethod
    def neighborhood_embeddings(embedder, neighbors_count: int) -> Callable[[ChainRecord], Union[ndarray[ndarray[float]], None]]:
        return DataAccessors.neighborhood_with_custom_embeddings(
            embedder, transform_embeddings_func=None, neighbors_count=neighbors_count)
    
    @staticmethod
    def neighborhood_with_custom_embeddings(embedder, transform_embeddings_func, neighbors_count: int) -> Callable[[ChainRecord], Union[ndarray[ndarray[float]], None]]:
        
        def get_neighborhood_embeddings(chain: ChainRecord):
            chain_structure = chain.get_chain_structure(loaded=False)
            
            # if cache is full this is not needed
            # chain_structure.load(preferred_sequence=chain.sequence())

            embeddings = chain.embeddings(embedder)
            
            if transform_embeddings_func is not None: 
                embeddings = transform_embeddings_func(embeddings)

            extended_embeddings = []  

            with use_cache(chain_structure) as cached_chain_structure:
                for i in range(cached_chain_structure.get_residue_count()):
                    neighbors = cached_chain_structure.get_nearest_residue_indexes(i)[:neighbors_count]

                    extended_embeddings.append(
                        np.array([embeddings[index] for index in neighbors]).flatten()
                    )

            return np.array(extended_embeddings)

        return get_neighborhood_embeddings
    
    @staticmethod
    def average_neighborhood_embeddings(embedder, neighbors_count: int) -> Callable[[ChainRecord], Union[ndarray[ndarray[float]], None]]:
        
        def get_neighborhood_embeddings(chain: ChainRecord):
            chain_structure = chain.get_chain_structure(loaded=False)
            
            # if cache is full this is not needed
            # chain_structure.load(preferred_sequence=chain.sequence())

            embeddings = chain.embeddings(embedder)
            extended_embeddings = []  

            with use_cache(chain_structure) as cached_chain_structure:
                for i in range(cached_chain_structure.get_residue_count()):
                    neighbors = cached_chain_structure.get_nearest_residue_indexes(i)[:neighbors_count]

                    sum_embeddings = embeddings[neighbors[0]]

                    for neighbor_index in neighbors[1:]:
                        sum_embeddings = sum_embeddings + embeddings[neighbor_index]

                    sum_embeddings = sum_embeddings / neighbors_count


                    extended_embeddings.append(
                        sum_embeddings
                    )

            return np.array(extended_embeddings)

        return get_neighborhood_embeddings