try:
    import torch
    
    x = torch.tensor([1, 2, 3])

    print("PyTorch is installed successfully.")
    print(f"CUDA GPUs count: {torch.cuda.device_count()}")
    print("Sample Tensor:", x)
except ImportError:
    print("PyTorch is not installed. Install it using: pip install torch")
