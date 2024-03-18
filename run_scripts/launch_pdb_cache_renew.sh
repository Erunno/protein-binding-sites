#!/bin/bash

total_proteins=5570
batch_size=100

total_batches=$((total_proteins / batch_size + 1))

echo removing cache ...
rm /home/brabecm4/diplomka/protein-binding-sites/data/cache_data/production/*

echo launching $total_batches workers

for ((i = 0; i < total_batches; i++)); do
    start_idx=$((i * batch_size))
    sbatch /home/brabecm4/diplomka/protein-binding-sites/run_scripts/run_pdb_cache_renew.sh "$start_idx" "$batch_size"
done