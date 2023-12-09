import torch
import torch.nn as nn
import torch.optim as optim
import data_loader as dl

binding_sights_db_filename = r'C:\Users\mbrabec\Desktop\MFF\diplomka\neural_netw_emb\data\binding_sights\binding_sights_by_ligand.json'
dataset_db_filename = r'C:\Users\mbrabec\Desktop\MFF\diplomka\neural_netw_emb\data\yu_datasets\proteins_by_datasets.json'
embeddings_folder = r'C:\Users\mbrabec\Desktop\MFF\diplomka\neural_netw_emb\data\embedded_sequences\unzipped'
ligand = 'AMP'

data_loader = dl.DataLoader(
    binding_sights_db_fname=binding_sights_db_filename,
    dataset_by_ligands_db_fname=dataset_db_filename,
    embeddings_folder=embeddings_folder
)

X_train, y_train, X_test, y_test = data_loader.get_data_set_for('AMP')


class ProteinClassifier(nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim):
        super(ProteinClassifier, self).__init__()
        self.hidden_dim = hidden_dim
        
        # Define layers
        self.fc1 = nn.Linear(input_dim, hidden_dim)  # Fully connected layer 1
        self.relu = nn.ReLU()  # Activation function
        self.fc2 = nn.Linear(hidden_dim, output_dim)  # Fully connected layer 2
        self.sigmoid = nn.Sigmoid()  # Activation function for binary classification

    def forward(self, x):
        x = self.fc1(x)
        x = self.relu(x)
        x = self.fc2(x)
        x = self.sigmoid(x)
        return x

# Define parameters
input_dim =  len(X_train[0][0])
hidden_dim = 128  # Number of neurons in the hidden layer
output_dim = 1  # Binary classification for binding site or not

# Initialize model
model = ProteinClassifier(input_dim, hidden_dim, output_dim)

# Define loss function and optimizer
criterion = nn.BCELoss()  # Binary Cross Entropy Loss
optimizer = optim.Adam(model.parameters(), lr=0.001)  # Adam optimizer

print('shape X: ', X_train.shape)
print('shape y: ', y_train.shape)

# Convert your data into tensors (Assuming 'embeddings' is your input and 'labels' is your target)
# Replace 'embeddings' and 'labels' with your actual data
inputs = torch.tensor(X_train[:10], dtype=torch.float32)
targets = torch.tensor(y_train[:10], dtype=torch.float32)

# Training loop
num_epochs = 10
for epoch in range(num_epochs):
    optimizer.zero_grad()  # Clear gradients
    
    # Forward pass
    outputs = model(inputs)
    loss = criterion(outputs, targets)
    
    # Backward pass and optimization
    loss.backward()
    optimizer.step()
    
    # Print loss every epoch
    print(f"Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}")