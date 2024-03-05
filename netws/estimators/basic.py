import numpy as np
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.utils.validation import check_X_y, check_array, check_is_fitted
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.init as init
from sklearn.preprocessing import StandardScaler

loading_bar_width = 25
UP_char = "\033[A"

class BasicNetwork(BaseEstimator, RegressorMixin):
    def __init__(self, 
                 input_size, hidden_sizes,
                 learning_rate,
                 epochs, batch_size):
        
        # Initialize the sequential model
        self.input_size = input_size
        self.hidden_sizes = hidden_sizes
        self.model = None
        self.optimizer = None
        self.criterion = nn.BCELoss()
        self.epochs = epochs
        self.batch_size = batch_size
        self.learning_rate = learning_rate

        self.device = 'cuda' if torch.cuda.is_available() else 'cpu'
        
        self.epoch_callback = None
        self.name = "BasicNetwork"
        self.verbose = False

    def set_name(self, name):
        self.name = name

    def set_verbose(self):
        self.verbose = True

    def get_name(self):
        return self.name

    def register_epoch_callback(self, callback):
        self.epoch_callback = callback

    def fit(self, X, y):
        if (self.verbose):
            print(f'info ({self.name}): fitting...', flush=True)

        # # Input validation
        X, y = check_X_y(X, y)

        # # Standardize input data
        scaler = StandardScaler()
        X = scaler.fit_transform(X)

        # Convert numpy arrays to PyTorch tensors
        X_tensor = torch.tensor(X, dtype=torch.float32).to(self.device)
        y_tensor = torch.tensor(y, dtype=torch.float32).view(-1, 1).to(self.device)

        # Define a PyTorch DataLoader for batching
        dataset = torch.utils.data.TensorDataset(X_tensor, y_tensor)
        dataloader = torch.utils.data.DataLoader(
            dataset, 
            batch_size=self.batch_size, shuffle=True)

        # Reset sequential model
        self._init_model()

        if (self.verbose):
            print('')

        # Training loop
        for epoch in range(self.epochs):
            if (self.verbose):
                bar_count = int((epoch + 1) * loading_bar_width / self.epochs)
                space_count = loading_bar_width - bar_count
                print(f'\r{UP_char} ({self.name}) epoch ({(epoch + 1):3d}/{self.epochs:3d}) [{"="*bar_count}{" "*space_count}]', flush=True)
            
            for batch_X, batch_y in dataloader:
                # Zero the gradients
                self.optimizer.zero_grad()

                # Forward pass
                predictions = self.model(batch_X)
                loss = self.criterion(predictions, batch_y)

                # Backward pass and optimization
                loss.backward()
                self.optimizer.step()

            if self.epoch_callback is not None:
                self.epoch_callback(epoch, loss, self)
        

        return self

    def predict(self, X):
        if (self.verbose):
            print(f'info ({self.name}): predicting...', flush=True)

        # Input validation
        X = check_array(X)

        # Standardize input data
        X = StandardScaler().fit_transform(X)
        X = torch.tensor(X, dtype=torch.float32).to(self.device)

        # Make predictions
        with torch.no_grad():
            predictions = self.model(X.to(self.device))[:, 0].cpu().round()

        return predictions

    def underling_model(self):
        return self.model
    
    def _init_model(self):
        self.model = self._create_simple_sequential_model(
            self.input_size, 1, self.hidden_sizes)
        
        self.optimizer = optim.Adam(self.model.parameters(), lr=self.learning_rate)

        self.model.to(self.device)

    def _create_simple_sequential_model(self, input_size, output_size, hidden_sizes):
        layers = []
        layer_sizes = [input_size] + hidden_sizes + [output_size]

        for i in range(len(layer_sizes) - 1):
            layers.append(nn.Linear(layer_sizes[i], layer_sizes[i + 1]))
            is_last_layer = i == len(layer_sizes) - 2

            if is_last_layer:
                layers.append(nn.Sigmoid())
            else:
                layers.append(nn.ReLU())

        return nn.Sequential(*layers)
