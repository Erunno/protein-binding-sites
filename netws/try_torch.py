try:
    import torch
    
    x = torch.tensor([1, 2, 3])

    print("PyTorch is installed successfully.")
    print("Sample Tensor:", x)
except ImportError:
    print("PyTorch is not installed. Please install it using: pip install torch")
