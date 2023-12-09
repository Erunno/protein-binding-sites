binding_sights_db_filename = r'C:\Users\mbrabec\Desktop\MFF\diplomka\neural_netw_emb\data\binding_sights\binding_sights_by_ligand.json'
dataset_db_filename = r'C:\Users\mbrabec\Desktop\MFF\diplomka\neural_netw_emb\data\yu_datasets\proteins_by_datasets.json'
embeddings_folder = r'C:\Users\mbrabec\Desktop\MFF\diplomka\neural_netw_emb\data\embedded_sequences\unzipped'
ligand = 'AMP'

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import data_loader as dl
from sklearn.metrics import matthews_corrcoef


data_loader = dl.DataLoader(
    binding_sights_db_fname=binding_sights_db_filename,
    dataset_by_ligands_db_fname=dataset_db_filename,
    embeddings_folder=embeddings_folder
)

X_train, y_train, X_test, y_test = data_loader.get_data_set_for('AMP')

X_train = torch.tensor(X_train, dtype=torch.float32)
y_train = torch.tensor(y_train, dtype=torch.float32)
 
X_test = torch.tensor(X_test, dtype=torch.float32)
y_test = torch.tensor(y_test, dtype=torch.float32)
 
# define the model
model = nn.Sequential(
    nn.Linear(len(X_train[0]), 1200),
    nn.ReLU(),
    nn.Linear(1200, 800),
    nn.ReLU(),
    nn.Linear(800, 1),
    nn.Sigmoid()
)
print(model)
 
# train the model
loss_fn   = nn.BCELoss()  # binary cross entropy
optimizer = optim.Adam(model.parameters(), lr=0.001)
 
n_epochs = 300
batch_size = 1000
 
for epoch in range(n_epochs):
    for i in range(0, len(X_train), batch_size):
        X_batch = X_train[i:i+batch_size]
        
        y_pred = model(X_batch)[:, 0]
        y_batch = y_train[i:i+batch_size]
        
        loss = loss_fn(y_pred, y_batch)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    print(f'Epoch: {epoch}, loss: {loss :.5f}')
 
with torch.no_grad():
    y_pred = model(X_test)[:, 0]
accuracy = (y_pred.round() == y_test).float().mean()


num_ones = np.count_nonzero(y_pred.round() == 1)
num_zeros = np.count_nonzero(y_pred.round() == 0)
print('zeros:', num_zeros, 'ones:', num_ones)

print(f"Accuracy {accuracy * 100:.5f}")


print(f'MCC: {matthews_corrcoef(y_test, y_pred.round())}')

# MCC