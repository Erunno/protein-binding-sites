#!/bin/bash

ligands=("ADP" "AMP" "ATP" "CA" "DNA" "FE" "GDP" "GTP" "HEME" "MG" "MN" "ZN")
# ligands=("ZN")

script="/home/brabecm4/diplomka/protein-binding-sites/run_scripts/p2rank/run_p2rank_traineval.sh"

for ligand in "${ligands[@]}"; do
    sbatch "$script" "$ligand"
done
