import torch
import random
import numpy as np

def seed_all(seed_value):
    # Set seed for Python's random number generator
    random.seed(seed_value)

    # Set seed for NumPy
    np.random.seed(seed_value)

    # Set seed for PyTorch
    torch.manual_seed(seed_value)
    torch.cuda.manual_seed_all(seed_value)  # If using CUDA

    # Additional steps for ensuring deterministic behavior on CUDA devices
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
