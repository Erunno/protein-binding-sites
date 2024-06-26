import numpy as np
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.metrics import matthews_corrcoef
from sklearn.utils.validation import check_X_y, check_array
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.preprocessing import StandardScaler

loading_bar_width = 25
UP_char = "\033[A"

class TrainerBase(BaseEstimator, RegressorMixin):
    def __init__(self,
                 learning_rate,
                 epochs, batch_size):
        
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
        self.best_threshold = None

    def set_name(self, name):
        self.name = name

    def set_verbose(self):
        self.verbose = True

    def get_name(self):
        return self.name

    def register_epoch_callback(self, callback):
        self.epoch_callback = callback

    def fit(self, X, y):
        orig_X = X
        if (self.verbose):
            print(f'info ({self.name}): fitting...', flush=True)

        # Input validation
        X, y = check_X_y(X, y)

        # Standardize input data
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

        self.optimize_threshold_for_mcc(orig_X, y)

        return self
    
    def optimize_threshold_for_mcc(self, X, y):
        y_pred_raw = self.predict_raw(X)
        best_mcc = -1

        for threshold in np.arange(0.0, 1.0, 0.05):
            pred = (np.array(y_pred_raw) >= threshold).astype(int)
            mcc = matthews_corrcoef(y, pred)

            if mcc > best_mcc:
                best_mcc = mcc
                self.best_threshold = threshold

    def predict(self, X):
        predictions = self.predict_raw(X)

        return (np.array(predictions) >= self.best_threshold).astype(int)

    def predict_raw(self, X):
        if (self.verbose):
            print(f'info ({self.name}): predicting...', flush=True)

        # Input validation
        X = check_array(X)

        # Standardize input data
        X = StandardScaler().fit_transform(X)
        X = torch.tensor(X, dtype=torch.float32).to(self.device)

        # Make predictions
        with torch.no_grad():
            predictions = self.model(X.to(self.device))[:, 0].cpu()

        return predictions

    def underling_model(self):
        return self.model
    
    def _init_model(self):
        self.model = self._get_model_to_train()
        self.optimizer = optim.Adam(self.model.parameters(), lr=self.learning_rate)
        self.model.to(self.device)

    def _get_model_to_train(self):
        raise Exception('Not implemented - should be implemented by the child class')
    