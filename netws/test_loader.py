import data_loader as dl

binding_sights_db_filename = r'C:\Users\mbrabec\Desktop\MFF\diplomka\neural_netw_emb\data\binding_sights\binding_sights_by_ligand.json'
dataset_db_filename = r'C:\Users\mbrabec\Desktop\MFF\diplomka\neural_netw_emb\data\yu_datasets\proteins_by_datasets.json'
embeddings_folder = r'C:\Users\mbrabec\Desktop\MFF\diplomka\neural_netw_emb\data\embedded_sequences\unzipped'

data_loader = dl.DataLoader(
    binding_sights_db_fname=binding_sights_db_filename,
    dataset_by_ligands_db_fname=dataset_db_filename,
    embeddings_folder=embeddings_folder
)

X_train, y_train, X_test, y_test = data_loader.get_data_set_for('AMP')



print('loaded')