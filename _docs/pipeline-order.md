# Pipeline order

1. **Prepare fasta files**: The initial step involves generating FASTA files from an unconventional sequence format, comprising *protein IDs*, *binding sites*, and the *protein sequence* itself. The Yu dataset is distributed across multiple files throughout the base directory, necessitating the aggregation of all existing proteins. This task is accomplished by a script named `create_fasta_files_from_yu.py`. The script systematically processes all files within the Yu dataset, producing multiple output files.

2. **Execute Embedder**: The subsequent step entails embedding all harvested sequences using an API. The `run_embedder.sh` script orchestrates the API calls for all files within the specified folder. Usage: `./run_embedder.sh <input_files_folder> <output_folder>`. To modify the embedder in use, access the script and adjust the designated variable.

3. **Unzip Archives**: To proceed, the script `unzip_all_embeddings.sh` has been crafted to extract all retrieved archives. You can utilize it by employing the following command structure: `unzip_all_embeddings.sh <input_archives_folder> <target_folder>`.

4. **Create Binding Sites Database**: Transforming the Yu dataset into an easily parsed JSON format serves as a database for neural network training. create_binding_sights_db.py processes the dataset's folder structure and multiple CSV-like files, generating a JSON object. This object represents binding sites with a defined format. Adjust the parameters within the script itself.
   ```json
   {
        "<Ligand name>": {
            "<protein id>": [
                // binding sites (zero based indexing)
                42, 12 //...
            ]
            // ...
        },
        "<Ligand name>": {
            "<protein id>": [
                // binding sites (zero based indexing)
                69, 13 //...
            ]
            // ...
        },
        // ...
   }
   ```

5. **Create YU Datasets by Protein Database**: This step involves organizing proteins into training and testing databases based on the original data's folder structure. The create_yu_datasets_lists.py script generates a JSON database by parsing the root folder of the Yu dataset. The output comprises a structured JSON format for training and testing protein databases. Adjust parameters within the script as needed.
   ```json
   {
        "<Ligand name>": {
            "train": [
                "<protein id>", "<protein id>" //...
            ],
            "test": [
                "<protein id>", "<protein id>" //...
            ]
        },
        "<Ligand name>": {
            "train": [
                "<protein id>", "<protein id>" //...
            ],
            "test": [
                "<protein id>", "<protein id>" //...
            ]
        },
        // ...
   }
   ```
