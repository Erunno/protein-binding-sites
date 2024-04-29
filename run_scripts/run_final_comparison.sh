#!/bin/bash

echo "Launching workers"

script_path="/home/brabecm4/diplomka/protein-binding-sites/run_scripts/launch_final_comparison.sh"

tags=("nei_compr_3_small_v4_c" "nei_compr_5_small_v4_c")
ligands=("GDP" "CA" "FE" "ADP" "DNA" "ZN" "MG" "AMP" "HEME" "MN" "ATP" "GTP")

# tags=("nei_compr_5_small_v4_c")
# ligands=("MG" "HEME")

radii=("1.0" "1.5" "2.0" "2.5" "3.0" "3.5" "4.0" "4.5" "5.5" "6.0" "6.5" "7.0" "7.5" "8.0" "8.5" "9.0" "9.5" "10.0" "10.5" )


for tag in "${tags[@]}"; do
    for ligand in "${ligands[@]}"; do
        sbatch "$script_path" "$tag" "$ligand"
    done
done

# for tag in "${tags[@]}"; do
#     for radius in "${radii[@]}"; do
#         sbatch "$script_path" "$tag" "ATP" "$radius"
#     done
# done

echo "Workers Launched"
