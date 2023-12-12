import numpy as np
import torch
import torch.nn as nn

def create_model(input_size, output_size, hidden_sizes):
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