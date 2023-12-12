import torch
import numpy as np
from sklearn.metrics import matthews_corrcoef

def get_statistics(model, X, y, epochs, print_res=False):
    with torch.no_grad():
        y_pred = model(X)[:, 0]

    accuracy = (y_pred.round() == y).float().mean() * 100
    mcc = matthews_corrcoef(y, y_pred.round())

    if (print_res):
        print(f"\nAccuracy: {accuracy:.5f}")
        print(f'MCC: {mcc}\n')


    return {
        'epochs': epochs,
        'acc': float(accuracy),
        'mcc': mcc
    }