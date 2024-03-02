import re
import datasets_db
import traceback

# constant parameters
path_to_protrusion_file = '/home/brabecm4/diplomka/protein-binding-sites/data/3d_proc/protrusion.max-neighbors.big.json'
path_to_embeddings_folder = '/home/brabecm4/diplomka/protein-binding-sites/data/embedded_sequences'
path_to_sequences_folder = '/home/brabecm4/diplomka/protein-binding-sites/data/orig/yu_sequences'
tested_embedder = 'ESM'

# constants
ALL_LIGANS = ['ADP', 'AMP', 'ATP', 'CA', 'DNA', 'FE', 'GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN']

RESET = '\033[0m'
RED = '\033[91m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
SPACE = '                                 '

# init database

print(f'{YELLOW}Initializing...')
db = datasets_db.SeqDatasetDb(
    sequences_folder=path_to_sequences_folder, 
    embeddings_folder=path_to_embeddings_folder)

print(f'Loading protrusion...{RESET}')
db.load_protrusion_data_file(path_to_protrusion_file)

def can_load_all_datasets():
    all_chains = db.get_all_chain_records()

    if all_chains is None or len(all_chains) == 0:
        return False, 'No chains found'

    return True, None
    
def can_access_protein_id(): 
    all_chains = db.get_all_chain_records()

    for chain in all_chains:
        prot_id = chain.protein_id()

        if prot_id is None or prot_id.strip() == '':
            return False, f'Protein id was "{prot_id}"'
    
    return True, None

def can_access_chain_id(): 
    all_chains = db.get_all_chain_records()

    for chain in all_chains:
        ch_id = chain.chain_id()

        if ch_id is None or ch_id.strip() == '':
            return False, f'Chain id was "{ch_id}"'

    return True, None

def can_access_full_id(): 
    all_chains = db.get_all_chain_records()

    for chain in all_chains:
        full_id = chain.full_id()

        if full_id is None or full_id.strip() == '':
            return False, f'Chain id was "{full_id}"'

    return True, None

def protein_id_has_correct_format(): 
    all_chains = db.get_all_chain_records()

    for chain in all_chains:
        prot_id = chain.protein_id()

        if not re.match('^[a-z0-9]{4}$', prot_id):
            return False, f'Protein id had incorrect format: "{prot_id}"'

    return True, None


def chain_id_has_correct_format(): 
    all_chains = db.get_all_chain_records()

    for chain in all_chains:
        ch_id = chain.chain_id()

        if not re.match('^[A-Z]$', ch_id):
            return False, f'Chain id had incorrect format: "{ch_id}"'
        
    return True, None

def full_id_has_correct_format(): 
    all_chains = db.get_all_chain_records()

    for chain in all_chains:
        full_id = chain.full_id()

        if not re.match('^[a-z0-9]{4}[A-Z]$', full_id):
            return False, f'Full id had incorrect format: "{full_id}"'
        
    return True, None

def there_are_some_chains_with_protrusion_records():
    all_chains = db.get_all_chain_records()
    all_chains = datasets_db.Helpers.filter_chains_with_protrusion(all_chains)

    if len(all_chains) == 0:
        return False, 'No chains with protrusion record'

    return True, None 

def protrusion_radii_are_correct(): 
    all_chains = db.get_all_chain_records()
    all_chains = datasets_db.Helpers.filter_chains_with_protrusion(all_chains)

    for chain in all_chains:
        radii = chain.get_all_protrusion_radii()

        if radii is None or len(radii) == 0:
            return False, f'Unexpected radii list: "{radii}"'      

    return True, None

def all_ligand_datasets_not_empty():
    for ligand in ALL_LIGANS:
        ds = db.get_data_set_for(ligand)

        testing = ds.testing()
        training = ds.training()

        if testing is None or len(testing) == 0:
            return False, f'Unexpected testing data for ligand {ligand}, got "{testing}"'

        if training is None or len(training) == 0:
            return False, f'Unexpected testing data for ligand {ligand}, got "{training}"'

    return True, None

def protrusion_does_not_throw_exception(): 
    all_chains = db.get_all_chain_records()
    all_chains = datasets_db.Helpers.filter_chains_with_protrusion(all_chains)

    for chain in all_chains:
        radii = chain.get_all_protrusion_radii()

        for r in radii: 
            protrusion = chain.protrusion_vector_for(r)

            if protrusion is None or len(protrusion) == 0:
                return False, f'Unexpected protrusion values: "{protrusion}"'        

    return True, None

def all_radii_function_on_db_returns_correct_value():
    ligand_ds = db.get_data_set_for('AMP')
    all_chains = datasets_db.Helpers.filter_chains_with_protrusion(ligand_ds.all()) 

    expected_all_radii = sorted(all_chains[0].get_all_protrusion_radii())
    actual_all_radii = ligand_ds.get_all_radii()

    are_same = all(x == y for x, y in zip(expected_all_radii, actual_all_radii))

    if not are_same:
        return False, f'Get all radii failed - expected {expected_all_radii}, got {actual_all_radii}'

    return True, None

def can_load_all_embeddings():
    all_chains = db.get_all_chain_records()

    i = 0
    for chain in all_chains:
        embeddings = chain.embeddings(tested_embedder)
        expected_len = len(chain.sequence())

        i += 1
        if i % 20 == 0:
            print(f'\r  {BLUE}checking embeddings - {YELLOW}{i}{BLUE} checked out of {YELLOW}{len(all_chains)}{BLUE}{RESET}', end='', flush=True)

        if embeddings is None or len(embeddings) != expected_len:
            return False, f'Unexpected embeddings: len={len(embeddings)} should be {expected_len} \n"{embeddings}"'
        
        chain.free_embeddings_from_RAM()

    return True, None

# run test cases

def run(test_case):
    print(f'{BLUE}  ...running {test_case.__name__}{RESET}', end='', flush=True)
    try:
        success, err_message = test_case()
    except Exception as e:
        success = False
        err_message = f'Exception occurred: "{e}"\n\n{traceback.format_exc()}'

    if success:
        print(f'\r{GREEN}[✔] {test_case.__name__}{RESET}{SPACE}')
    else:
        print(f'\r{RED}[✘] {test_case.__name__}{RESET}{SPACE}')
        print(f'  Error: {err_message}')
        

print(f'\n{BLUE}Running tests...{RESET}\n')

run(can_load_all_datasets)
run(can_access_protein_id)
run(can_access_chain_id)
run(can_access_full_id)
run(protein_id_has_correct_format)
run(chain_id_has_correct_format)
run(full_id_has_correct_format)
run(all_ligand_datasets_not_empty)
run(there_are_some_chains_with_protrusion_records)
run(protrusion_radii_are_correct)
run(protrusion_does_not_throw_exception)
run(all_radii_function_on_db_returns_correct_value)
run(can_load_all_embeddings)

print()