import argparse


def get_args():
    allowed_ligands = ['ADP', 'AMP', 'ATP', 'CA', 'DNA', 'FE', 'GDP', 'GTP', 'HEME', 'MG', 'MN', 'ZN']
    allowed_embedder = ['BERT', 'ESM', 'T5']

    parser = argparse.ArgumentParser(description='Ligand binding sites neural network')
    
    parser.add_argument('--ligand', type=str, choices=allowed_ligands, help='Name of the ligand')
    parser.add_argument('--hidden-layers', type=int, nargs='+', help='List of hidden layer sizes')
    parser.add_argument('--batch-size', type=int, help='Size of one batch')
    parser.add_argument('--epochs', type=int, help='Number of epochs')
    parser.add_argument('--epoch-stats-interval', type=int, help='Interval of gathering statistics.')
    parser.add_argument('--seed', type=int, help='Seed of random.')
    parser.add_argument('--learning-rate', type=float, help='Learning rate')
    parser.add_argument('--verbose', type=bool, default=False, help='Print intermediate results.')
    parser.add_argument('--embedder', type=str, choices=allowed_embedder, help='Embedder to be used.')

    return parser.parse_args()
