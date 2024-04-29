# p2rank_prep

This folder contains scripts that prepare data for the [P2Rank](https://github.com/rdk/p2rank) tool and then collects results from the training and testing phases.

- [`p2rank_dataset_preparation.py`](./p2rank_dataset_preparation.py) is used to convert the Yu dataset format into `.ds` files, as expected by [P2Rank](https://github.com/rdk/p2rank).
- [`show_results.py`](./show_results.py) collects results and stores them in the same format as the results collected from our models.

Scripts that have been used for the execution of these scripts are stored in [~/run_scripts/p2rank/](../run_scripts/p2rank/).
