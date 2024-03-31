import torch
import numpy as np
from sklearn.metrics import matthews_corrcoef

def get_statistics(model, X, y, epochs, print_res=False):
    y_pred = model.predict(X)
    
    y = np.array(y).round()
    y_pred = np.array(y_pred).round()
    accuracy = (y_pred == y).mean() * 100
    mcc = matthews_corrcoef(y, y_pred)

    if (print_res):
        print(f"\nAccuracy: {accuracy:.5f}")
        print(f'MCC: {mcc} (th: {model.best_threshold})')
        print(f'  zeros: {(y_pred == 0).sum()}, ones: {(y_pred == 1).sum()}\n')

    return {
        'epochs': epochs,
        'acc': float(accuracy),
        'mcc': mcc,
        'threshold': model.best_threshold
    }
