```bash
python .\network.py --hidden-layers 900 400 200 50 --seed 42 --ligand AMP --learning-rate 0.001 --epochs 300 --batch-size 1000 --epoch-stats-interval 10 --verbose True --embedder BERT --protrusion-data-file ..\data\3d_proc\protrusions.big.json --pdb-mappings-fname ..\data\3d_proc\mappings_to_pdbs.json --tag not_specified
```