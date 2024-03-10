import numpy as np
from sklearn.base import BaseEstimator, RegressorMixin

class PartialInputNetwork(BaseEstimator, RegressorMixin):
    def __init__(self, input_size, model):
        self.model = model
        self.input_size = input_size
        
    def set_name(self, name):
        self.model.set_name(name)

    def set_verbose(self):
        self.model.set_verbose()

    def get_name(self):
        return self.model.get_name()

    def register_epoch_callback(self, callback):
        self.model.register_epoch_callback(callback)

    def fit(self, X, y):
        X = np.array(X)
        y = np.array(y)

        X = X[:, :self.input_size]

        self.model.fit(X, y)
        
    def predict(self, X):
        X = np.array(X)
        X = X[:, :self.input_size]

        return self.model.predict(X)

    def underling_model(self):
        return self.model.underling_model()