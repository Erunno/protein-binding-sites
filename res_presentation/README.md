# res_presentation

This folder contains multiple scripts aimed at presenting the results from the model training and evaluation.

- [`create_result_table.py`](./create_result_table.py) generates tables with results based on different tags.
- [`compare_final.py`](./compare_final.py) provides a CLI tool for comparing models. It creates a table similar to the one in the thesis for model comparison.
- [`model_comparison_table.py`](./model_comparison_table.py) generates a model comparison table. This is an older version used during model development for ad-hoc comparisons.
- [`harvest_best_hp.py`](./harvest_best_hp.py) goes through all results produced by the training and harvests the best hyperparameters for each tag and ligand.
- [`latex_table_printer.py`](./latex_table_printer.py) prints LaTeX comparisons that are presented in the thesis.
- [`analyze_performance.py`](./analyze_performance.py) creates graphs that show MCC (Matthews Correlation Coefficient) over time.

### Helpers:

- [`results_loader.py`](./results_loader.py) loads all results that have been produced by the training.
- [`table_printer.py`](./table_printer.py) prints any table in a pretty format.
