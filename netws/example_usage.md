```bash
python3 ./network.py --hidden-layers 200 50 --seed 42 --ligand AMP --learning-rate 0.01 --epochs 15 --batch-size 1000 --epoch-stats-interval 5 --verbose True --embedder BERT --protrusion-data-file ../data/3d_proc/protrusions.big.json --pdb-mappings-fname ../data/3d_proc/mappings_to_pdbs.json --tag not_specified
```

```bash
python /home/brabecm4/diplomka/protein-binding-sites/res_presentation/create_result_table.py --results-folder /home/brabecm4/diplomka/protein-binding-sites/data/netw_results --embedder-aliases ProtT5:-T5 ProtBert:-BERT --compare
```
