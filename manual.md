# Manual

This manual describes the steps that need to be taken to train your own neural networks as described in the thesis.

## Libraries Required

These are the libraries you need to install:
- Pytorch
- SciPy
- numpy
- BioPython
- Levenshtein
- matplotlib (optional, only for plotting graphs)
- curses (optional, only for displaying interactive model comparison table)

## Training a Network Step by Step

Having installed all dependencies, we move on to setting up an environment to run the networks. Here is a high-level overview:

1. Get the Yu dataset
2. Get 3D protein structures
3. Get Protein Embeddings
4. Set up Config variables in [`config.py`](./config/config.py)
5. Fill the cache
6. Run tests
7. Train networks (or perform other analyses)
8. Inspect results

### Download Dataset and Other Data

The first step is to clone this main repository.

The second step is to download the entire Yu dataset, accessible at the [data repository](https://github.com/Erunno/protein-binding-sites-data). It contains sequences and needed pdb files (protein 3D structures).

### Obtaining Embeddings

Next, create a folder containing all embeddings of all chains present in the Yu dataset.

In the [data repository](https://github.com/Erunno/protein-binding-sites-data), there is a folder [`data repository/fasta`](https://github.com/Erunno/protein-binding-sites-data/tree/main/fasta) which contains a fasta file with all sequences of all protein chains used in the Yu dataset. This fasta file is used to obtain embeddings from the embedder of your choice (in our work, we used [ESM-2](https://github.com/facebookresearch/esm)). You have multiple options:

1. Run the script [`data_prep/run_embedder.sh`](./data_prep/run_embedder.sh) with the command: `sh data_prep/run_embedder.sh [esm|t5|bert] <directory containing fasta file(s)> <target directory of embeddings>`. This calls an API that runs the selected embedder automatically. All you need to do then is to unzip the obtained files and move them to one folder named after the embedder you selected.
2. Download and set up the embedder of your liking. The expected folder structure is `<folder with all your embeddings>/<your embedder name>/<Protein ID with a chain ID>.npy` (see example below). These files in the embedder folder should be numpy vectors.

### Loaded Data Checkpoint

At this point, you have all the data needed for our project. Your file system should look something like this:

```plaintext
root_folder:
    |- protein-binding-sites-data
        |- biolip_structures
            |- ...
        |- fasta
            |- all_chains_in_yu.fasta
        |- yu_dataset
            |- ...
        |- yu_sequences
            |- ...
    |- my_embeddings_store
        |- ESM
            |- 4sgbE.npy
            |- 7mdhA.npy
            |- 8icnA.npy
            |- ...
        |- T5
            |- ...
        |- ... (other embedders if you want to experiment...)
    |- protein-binding-sites
        |- config
            |- ...
        |- ...   
```

### Setup Config Variables

Next, update the [`config.py`](./config/config.py) script variables to point to the data we have downloaded and produced so far:

```python
data_top_folder = 'root/folder/of/your/project'

# present in the `protein-binding-sites-data` repository
yu_sequences_folder = f'{data_top_folder}/protein-binding-sites-data/yu_sequences'

# present in the `protein-binding-sites-data` repository
pdbs_folder = f'{data_top_folder}/protein-binding-sites-data/biolip_structures'

# has to be generated
embeddings_folder = f'{data_top_folder}/my_embeddings_store'
```

Also, reserve an empty folder for the cache:

```python
cache_folder = f'{data_top_folder}/cache_data_folder'
```

And the last variable we need to set is a folder where we are going to store reports (json objects) from the network training:

```python
networks_results_folder = f'{data_top_folder}/netw_results/netw_runs'
```

We do not have to concern ourselves with other variables for network training.

### Fill the Cache

Next, fill the cache with precomputed protrusion, SASA, and closest residues. Note that methods of the class [`ChainRecord`](./data_prep/datasets_db.py) concerning 3D structure won't work without a full cache—this is to avoid unwanted computation during network training.

To fill the cache, use the script [`data_prep/pdb_files_refresh_cache.py`](./data_prep/pdb_files_refresh_cache.py) without any parameters. See the actual methods of the class [`Chain3dStructure`](./data_prep/pdb_files_db.py) that are being cached in the method [`run_functions_to_cache`](./data_prep/pdb_files_refresh_cache.py), and add more if you want to experiment with more protrusion radii, for instance.

The [`data_prep/pdb_files_refresh_cache.py`](./data_prep/pdb_files_refresh_cache.py) script has multiple parameters like `--remove-old` that removes all records from the cache if any are present before. However, all we need is simply to run the script without any parameters to run the cache (by default, we precompute SASA, neighboring residues indexes, and protrusion for radii 1 to 10 in steps of 0.5).

### Run Tests

To see if everything worked, you should run the tests in files [`pdb_files_tests.py`](./pdb_files_tests.py) and [`datasets_tests.py`](./datasets_tests.py) using no parameters. The only modification that needs to be done is to add a path to the test cache in the file [`pdb_files_tests.py`](./pdb_files_tests.py)—this is done so that your production cache is not compromised by those tests.

### Run Network

Finally, you are ready to run the training of the network. Use the script [`netws/network.v2.py`](./netws/network.v2.py) with the desired hyperparameters. See the parameters that can be used in the script itself. Most of them are easy to follow; we will only mention the parameter `--tag`. This parameter tags writes the specific tag to the resulting report (json object produced by the script); this allows you to group your experiment under one tag. This mechanism is then used by the scripts displaying results like [`res_presentation/create_result_table.py`](./res_presentation/create_result_table.py), which allows you to see the best results obtained and compare them with other experiments under a different tag. A tag is semantically a group of the same experiment, and you can compare the performance of different tags.

Also, the [`netws/network.v2.py`](./netws/network.v2.py) as is runs the simple network with embeddings only (no 3D features, etc.). You can simply modify the script as you wish by uncommenting some of the data accessors (e.g., `dataset.DataAccessors.protrusion(3.5)` appends a protrusion value obtained from the radius 3.5 Å).

### See Results

The main script that prints a table with all results is [`res_presentation/create_result_table.py`](./res_presentation/create_result_table.py); it prints the best results. It has two modes: a mode with a flag `--compare`: the script provides a comparison table with results for each tag. Without it, the script prints the best results with details about hyperparameters used. You can also filter results for only selected tags using `--tag`, which takes a list of tags that is taken into consideration.

## More Functionality

Now that we have provided a foundational script, we are ready to explore more functionalities available in the repository. We assume the data are downloaded and the environment set up as described in the previous section.

### Network with Compression Layer (Composed Network)

As described in the thesis, we have implemented a network with a compression layer consisting of pre-trained networks. First, we need to identify the best baseline model serving as the pre-trained model. For this purpose, we implemented the script [`res_presentation/harvest_best_hp.py`](./res_presentation/harvest_best_hp.py) which harvests the best hyperparameters from the `results folder` and saves them to a JSON file. Both of these locations are defined by the [`config.py`](./config/config.py) script. We have already defined the `results folder`, and the JSON object is stored into a file defined by:

```python
# Define the proper place for best hyperparameters JSON file
best_HPs_file = f'{data_top_folder}/path/to/best_HPs.json'
```

Once you have defined the needed variable, simply run the [`res_presentation/harvest_best_hp.py`](./res_presentation/harvest_best_hp.py) script without parameters. We will need the file in other cases as well.

The main script for the layer with composed networks is [`netws/network.composed.v2.py`](./netws/network.composed.v2.py). The parameters are the same as those for the [`netws/network.v2.py`](./netws/network.v2.py) discussed above. You are only required to define a variable `tag_of_the_pretrained_model` in the script itself. This variable defines the tag of the model that is used as the pre-trained layer, that is to say, the best set of hyperparameters of this tag.

### Run Final Comparison

The final network training scripts are those used for the 5x2 cross-validation process described in the thesis. We had to divide these into two files: [`netws/network.final_eval.py`](./netws/network.final_eval.py) and [`netws/network.final_eval.composed.py`](./netws/network.final_eval.composed.py)—the first is for simpler networks and the other for the network with compression layers. However, both function similarly, i.e., they take the `--final-tag` parameter which defines for which tag the 5x2 cross-validation is performed (the best hyperparameters of the tag are used) and the `--ligand` which defines for which ligand we are running the model. Additionally, you are required to provide some values so that both scripts know how to define data based on the tag:

- [`netws/network.final_eval.py`](./netws/network.final_eval.py): variable `accessors` which defines data accessors for each tag.
- [`netws/network.final_eval.composed.py`](./netws/network.final_eval.composed.py): variable `neighbors` which defines the number of neighboring residues (minus one) used in the network.

Finally, all results from the comparisons are stored in the folder defined by the following variable in [`config.py`](./config/config.py):

```python
final_networks_results_folder = f'{data_top_folder}/path/to/folder/final_runs'
```

### See Results from the Final Comparison

Once you have run the 5x2 cross-validation, you may want to inspect the results. We implemented an interactive CLI tool present in the script [`res_presentation/compare_final.py`](./res_presentation/compare_final.py) that takes the parameter `--baseline-model-tag` which serves as the baseline model to which you can compare other models. The CLI is interactive (uses the package `curses`) so you can choose any model you like to compare with the one defined by `--baseline-model-tag` once you run it. Another option is to use the `--full-table` flag; then a static table is printed to the standard output.
