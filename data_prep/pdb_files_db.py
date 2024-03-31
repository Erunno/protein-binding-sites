import os
import sys
from typing import List, Type, Union
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import config.config as config
from Bio.PDB import PDBParser
from Bio.PDB import MMCIFParser
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.PDB.NeighborSearch import NeighborSearch
from Levenshtein import distance as lev_distance
from Bio.PDB.SASA import ShrakeRupley

class Chain3dStructure:
    def __init__(self, 
                 protein_id, chain_id,
                 structure_filename, 
                 load=False):
        
        self.__cache_key__ = f'{protein_id.lower()}{chain_id.upper()}'
        self.__cashable__ = [
            self.get_residue_count.__name__,
            self.get_nearest_residue_indexes.__name__,
            self.all_residues_has_alpha_carbon.__name__,
            self.compute_sequence.__name__,
            self.get_protrusion_vector.__name__,
            self.get_SASA_vector.__name__
        ]

        self.protein_id = protein_id
        self.chain_id = chain_id
        self.structure_file = structure_filename
        
        self._3d_structure: Union[Structure, None] = None
        self._3d_chain = None

        if load:
            self.load()

    def load(self, preferred_sequence=None):
        if (self.structure_file.endswith('cif')):
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)

        self._3d_structure = parser.get_structure(id, self.structure_file)

        found_chains = []

        for model in self._3d_structure:
            for chain in model.get_chains():
                if chain.id.lower() == self.chain_id.lower():
                    found_chains.append(chain)                    

        if len(found_chains) == 0:
            raise Exception(f"Chain with ID {self.chain_id} not found in the structure. Protein {self.protein_id} - file {self.structure_file}")

        if len(found_chains) == 1 or preferred_sequence is None:
            self._3d_chain = found_chains[0]
            return
        
        best_match_chain = found_chains[0]
        best_lev_dist = sys.maxsize
        
        for chain in found_chains:
            chain_seq = Chain3dStructure.compute_chain_sequence(chain)
            dist = lev_distance(chain_seq, preferred_sequence)

            if dist < best_lev_dist:
                best_lev_dist = dist
                best_match_chain = chain

        self._3d_chain = best_match_chain

    def get_residue_count(self) -> int:
        return len(self.get_residue_list())

    def get_nearest_residue_indexes(
            self, 
            residue_index: int, 
            
            # do not change these default values as some values may be stored in cache (alternatively delete cache)
            n_nearest: int = 10,
            distance_func='alpha_atoms') -> List[int]:

        if distance_func in DistancesCollection.names_to_functions:
            distance_func = DistancesCollection.names_to_functions[distance_func] 

        residue_list = self.get_residue_list()
        center_residue = residue_list[residue_index]

        distances = []

        for i in range(len(residue_list)):
            distances.append((
                i, distance_func(center_residue, residue_list[i])
            ))

        distances.sort(key=lambda x: x[1])

        return list([residue_index for residue_index, dist in distances[:n_nearest]])

    def free_memory(self):
        self._3d_structure = None
        self._3d_chain = None

    def get_residue_list(self) -> List[Residue]: 
        return Chain3dStructure.get_chain_residue_list(self._3d_chain) 

    def compute_sequence(self) -> str:
        return Chain3dStructure.compute_chain_sequence(self._3d_chain)

    def all_residues_has_alpha_carbon(self):
        for res in self.get_residue_list():
            if 'CA' not in res:
                return False
            
        return True

    def get_origin_filename(self):
        return self.structure_file

    def get_protrusion_vector(
            self, radius, 
            # do not change this default value as some values may be stored in cache (alternatively delete cache)
            protrusion_algorithm='atom_count_from_center_of_mass'):
        
        if protrusion_algorithm in ProtrusionFunctionsCollection.names_to_functions:
            protrusion_algorithm = ProtrusionFunctionsCollection.names_to_functions[protrusion_algorithm] 

        all_atoms = list(self._3d_chain.get_atoms())
        ns = NeighborSearch(all_atoms)

        protrusion = []

        for residue in self.get_residue_list():
            value = protrusion_algorithm(ns, residue, radius)
            protrusion.append(value)

        return protrusion
    
    def get_SASA_vector(self):
        sr = ShrakeRupley()
        sr.compute(self._3d_chain, level="R")
        
        return [
            res.sasa for res in self.get_residue_list()
        ]

    @staticmethod
    def get_chain_residue_list(biopython_chain) -> List[Residue]:
        # Exclude hetero residues or non-standard residues with insertion code ' '
        return list([residue for residue in biopython_chain if residue.get_id()[0] == " "])

    @staticmethod
    def compute_chain_sequence(biopython_chain) -> str:
        return ''.join([AminoAcidMapper.to_one_letter_code(
                            residue.get_resname())
                        for residue in Chain3dStructure.get_chain_residue_list(biopython_chain)])

allowed_file_types = ['pdb', 'ent', 'cif']

class PdbFilesDb:
    def __init__(self, storage_folder=config.pdbs_folder):
        all_files = self.__list_files(storage_folder)
        self.file_list = [file for file in all_files if file[-3:] in allowed_file_types]

        self.protein_id_parsers = self.__get_protein_id_parsers()

        self.protein_ids_to_file_mappings = self.__get_protein_id_to_files_mapping(self.file_list)

    def get_file_name_for(self, protein_id: str, chain_id) -> str: 
        potential_files = self.protein_ids_to_file_mappings[protein_id]

        full_id = f'{protein_id}{chain_id}'.lower()
        other_file = None

        for file in potential_files:
            if full_id in file.lower():
                return file
            
            if file.endswith('ent') or file.endswith('cif'):
                other_file = file

        return other_file

    def get_chain_structure(self, protein_id, chain_id) -> Chain3dStructure: 
        structure_file = self.get_file_name_for(protein_id, chain_id)
        return Chain3dStructure(protein_id=protein_id, 
                                chain_id=chain_id,
                                structure_filename=structure_file,
                                load=False)

    def __get_protein_id_parsers(self):
        def get_id_from_pdb(fname):
            return fname[-9:-5]
        def get_id_from_ent(fname):
            return fname[-8:-4]
        def get_id_from_cif(fname):
            return fname[-8:-4]

        return {
            'pdb': get_id_from_pdb,
            'ent': get_id_from_ent,
            'cif': get_id_from_cif,
        }
    
    def __get_prot_id(self, fname: str):
        id_parser = self.protein_id_parsers[fname[-3:]]
        return id_parser(fname)


    def __get_protein_id_to_files_mapping(self, files):
        mappings = {}

        for file in files: 
            id = self.__get_prot_id(file)

            if id not in mappings:
                mappings[id] = []

            mappings[id].append(file)

        return mappings


    def __list_files(self, folder):
        file_list = []
        for root, dirs, files in os.walk(folder):
            for file in files:
                file_list.append(os.path.join(root, file))
        return file_list


class AminoAcidMapper:
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
        'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
        'TYR': 'Y', 'VAL': 'V'
    }

    @staticmethod
    def to_one_letter_code(long_code) -> str:
        def map_name(name):
            if name in AminoAcidMapper.three_to_one:
                return AminoAcidMapper.three_to_one[name]
            return f'#({name})'
        
        return map_name(long_code)


class ResidueDistances:
    @staticmethod
    def distance_between_alpha_atoms_sq(residue_1: Residue, residue_2: Residue) -> float:
        if ('CA' not in residue_1) or ('CA' not in residue_2):
            return sys.maxsize

        coord_1 = residue_1['CA'].get_coord()
        coord_2 = residue_2['CA'].get_coord()

        return sum((a - b) * (a - b) for a, b in zip(coord_1, coord_2))

class ProtrusionFunctions:
    @staticmethod
    def compute_neighboring_atoms_from_center(ns, residue, radius):
        center = residue.center_of_mass()
        neighbors = ns.search(center=center, radius=radius, level='R') 
        return len(neighbors)

class DistancesCollection:
    names_to_functions = {
        'alpha_atoms': ResidueDistances.distance_between_alpha_atoms_sq
    }

class ProtrusionFunctionsCollection:
    names_to_functions = {
        'atom_count_from_center_of_mass': ProtrusionFunctions.compute_neighboring_atoms_from_center
    }
        