ligands=("ADP" "AMP" "ATP" "CA" "DNA" "FE" "GDP" "GTP" "HEME" "MG" "MN" "ZN")
# ligands=("AMP")
radius_values=("1.0" "2.0" "4.0" "6.0" "8.0" "10.0")

for ligand in "${ligands[@]}"; do
  for radius in "${radius_values[@]}"; do
    sbatch /home/brabecm4/diplomka/protein-binding-sites/run_scripts/run_comparison.sh "$ligand" "$radius"
  done
done