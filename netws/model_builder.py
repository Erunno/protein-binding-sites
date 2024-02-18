import numpy as np
import torch
import torch.nn as nn

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

def create_model(input_size, output_size, hidden_sizes, bypassed_last_inputs_count = 0):
    return BypassedInputsModel(input_size, output_size, hidden_sizes, bypassed_last_inputs_count)


def create_simple_model(input_size, output_size, hidden_sizes):
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