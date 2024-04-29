# Using Combined Sequence and Structural Features to Predict Protein-Ligand Binding Sites

This repository contains the code for the thesis completed at Charles University titled "*Using Combined Sequence and Structural Features to Predict Protein-Ligand Binding Sites*".

## Folder Structure

The goal was not to implement a single monolith application. Instead, the repository contains an ensemble of scripts, each fulfilling distinct roles such as defining machine learning models, preparing data, and presenting results from actual experiments. The scripts are organized into function-based folders, each containing a README for detailed explanations. Below is a high-level overview:

- **config**: Contains configuration values for all scripts.
- **data_prep**: Includes multiple Python scripts for raw data preparation and creating abstractions for subsequent use.
- **netws** (short for *networks*): Houses model definitions and scripts to operate these models.
- **p2rank_prep**: Focuses on preparing data for the *P2Rank* tool and collecting results from its training and testing phases.
- **res_presentation** (short for *results presentation*): Provides scripts for presenting experimental results.
- **run_scripts**: Contains infrastructure code for automatic job distribution using Slurm, developed for the ParLab cluster at Charles University.
- **stat**: Offers scripts that present various statistics about the dataset used.

## Other Repositories

This repository serves as the main codebase for the project. Additional repositories include one for storing the data used in this project ([data repository](https://github.com/Erunno/protein-binding-sites-data)) and another for archiving results from all conducted experiments ([results repository](https://github.com/Erunno/protein-binding-sites-results)).
