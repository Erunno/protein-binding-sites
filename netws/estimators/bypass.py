import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__))))
import numpy as np
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.utils.validation import check_X_y, check_array, check_is_fitted
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.init as init
from sklearn.preprocessing import StandardScaler
import trainer_base

class BypassedInputsModel(nn.Module):
    def __init__(self, input_size, output_size, hidden_sizes, bypassed_inputs_count):
        super(BypassedInputsModel, self).__init__()

        self.layers = nn.ModuleList()
        self.bypassed_inputs_count = bypassed_inputs_count

        layer_sizes = self._get_layer_sizes(input_size, output_size, hidden_sizes, bypassed_inputs_count)

        for i in range(len(layer_sizes) - 1):
            is_last_layer = i == len(layer_sizes) - 2

            if is_last_layer:
                self.layers.append(nn.Linear(layer_sizes[i] + bypassed_inputs_count, layer_sizes[i + 1]))
                self.layers.append(nn.Sigmoid())
            else:
                self.layers.append(nn.Linear(layer_sizes[i], layer_sizes[i + 1]))
                self.layers.append(nn.ReLU())

    def forward(self, x):
        saved_inputs = x[:, -self.bypassed_inputs_count:] 
        x = x[:, :-self.bypassed_inputs_count]

        for i, layer in enumerate(self.layers):
            is_last_layer = i == len(self.layers) - 2

            if is_last_layer:
                x = torch.cat((x, saved_inputs), dim=1)

            x = layer(x)

        return x
    
    def _get_layer_sizes(self, input_size, output_size, hidden_sizes, bypassed_inputs_count):
        layer_sizes = [input_size] + hidden_sizes + [output_size]

        layers_count = len(layer_sizes)

        layer_sizes[0] = layer_sizes[0] - bypassed_inputs_count
        layer_sizes[layers_count - 2] = layer_sizes[layers_count - 2] + bypassed_inputs_count

        return layer_sizes

class BypassedInputsNetwork(trainer_base.TrainerBase):
    def __init__(self, 
                 input_size, hidden_sizes,
                 bypassed_inputs,
                 learning_rate,
                 epochs, batch_size):
        
        super().__init__(learning_rate=learning_rate, epochs=epochs, batch_size=batch_size)

        self.bypassed_inputs = bypassed_inputs
        self.input_size = input_size
        self.hidden_sizes = hidden_sizes
        self.output_size = 1
        
    def _get_model_to_train(self):
        return BypassedInputsModel(
            self.input_size,
            self.output_size,
            self.hidden_sizes,
            bypassed_inputs_count=self.bypassed_inputs)