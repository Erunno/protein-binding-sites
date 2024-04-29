import torch
from estimators.basic import BasicNetwork
import torch.nn as nn
from sklearn.preprocessing import StandardScaler

def get_compressor_function(X, y, best_params_for_compressor):
    
    compressor = BasicNetwork(batch_size=best_params_for_compressor['batch_size'], 
                     input_size=len(X[0]),
                     hidden_sizes=best_params_for_compressor['layers'],
                     epochs=best_params_for_compressor['epochs'],
                     learning_rate=best_params_for_compressor['learning_rate'])
    compressor.set_verbose()

    print ('training compressor...')
    compressor.fit(X, y)

    print ('defining compressing function...')

    layers = list(compressor.underling_model().children())
    trained_embedder = nn.Sequential(*layers[:-2])

    def transform_embeddings(embeddings):
        X = StandardScaler().fit_transform(embeddings)
        X = torch.tensor(X, dtype=torch.float32).to('cuda')
        with torch.no_grad():
            return trained_embedder(X.to('cuda'))[:,:].cpu()
        
    return transform_embeddings
