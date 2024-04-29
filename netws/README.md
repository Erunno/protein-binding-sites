# netws

This folder contains the definitions of neural networks used in the thesis. Below are the main files that serve as executable scripts:

- [`network.v2.py`](./network.v2.py) is the primary entry point. It trains the network on training data and saves the results to a results store as defined in the config. The script needs to be modified to accommodate different models described in the thesis.
- [`network.composed.v2.py`](./network.composed.v2.py) serves as an entry to the "compressed" neural network design discussed in the thesis.

### Scripts for Testing and Validation:

- [`network.final_eval.py`](./network.final_eval.py) - Runs cross-validation for simpler network architectures.
- [`network.final_eval.composed.py`](./network.final_eval.composed.py) - Runs cross-validation for "compressed" network architectures.
- [`compare_estimators.py`](./compare_estimators.py) - Serves for ad-hoc model comparisons and was used during development.

## Estimators

The subfolder [`./estimators`](./estimators) contains the definitions of the neural networks and the training processes used.
