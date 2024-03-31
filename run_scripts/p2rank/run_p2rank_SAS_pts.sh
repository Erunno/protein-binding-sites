#!/bin/bash

#SBATCH -p gpu-long
#SBATCH -A nprg058s
#SBATCH --cpus-per-task=2
#SBATCH --time=0-02:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --signal=B:USR1@30
#SBATCH --output=/home/brabecm4/diplomka/protein-binding-sites/data/logs/p2rank/p2rank%A_%a.out

ds_folder=/home/brabecm4/diplomka/protein-binding-sites/data/p2rank_ds/ds_per_chain
p2rank_path=/home/brabecm4/diplomka/protein-binding-sites/p2rank/p2rank

start_idx=$1
count=$2

cd $p2rank_path

files_to_process=$(find "$ds_folder" -type f | sort | tail -n +$((start_idx+1)) | head -n $count)

i=0
for file in $files_to_process; do
    echo "Processing file $((i++)) out of $count"
    echo $file

    ./prank.sh traineval \
        -t "$file" \
        -e "$file" \
        -loop 1 -delete_vectors 0 -sample_negatives_from_decoys 0 \
        -features '(chem,volsite,protrusion,bfactor,xyz)'
done