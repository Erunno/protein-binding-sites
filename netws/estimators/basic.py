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

loading_bar_width = 25
UP_char = "\033[A"

class BasicNetwork(trainer_base.TrainerBase):
    def __init__(self, 
                 input_size, hidden_sizes,
                 learning_rate,
                 epochs, batch_size):
        
        super().__init__(learning_rate=learning_rate, epochs=epochs, batch_size=batch_size)

        self.input_size = input_size
        self.hidden_sizes = hidden_sizes
        self.output_size = 1
        
    def _get_model_to_train(self):
        layers = []
        layer_sizes = [self.input_size] + self.hidden_sizes + [self.output_size]

        for i in range(len(layer_sizes) - 1):
            layers.append(nn.Linear(layer_sizes[i], layer_sizes[i + 1]))
            is_last_layer = i == len(layer_sizes) - 2

            if is_last_layer:
                layers.append(nn.Sigmoid())
            else:
                layers.append(nn.ReLU())

        return nn.Sequential(*layers)
