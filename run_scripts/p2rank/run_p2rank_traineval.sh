#!/bin/sh

#SBATCH -p gpu-long
#SBATCH -A nprg058s
#SBATCH --gres=gpu:V100
#SBATCH --cpus-per-task=16
#SBATCH --time=0-10:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --signal=B:USR1@30
#SBATCH --output=/home/brabecm4/diplomka/protein-binding-sites/data/logs/p2rank/traineval_%A_%a.out

ligand=$1

echo computing $ligand ...

# analyze

types=("test" "train")

for type in "${types[@]}"; do
    
    ../p2rank/prank.sh -fail_fast false -c residues.groovy -out_subdir _$ligand analyze labeled-residues _$ligand/"$ligand"_"$type".ds
    
    file=/home/brabecm4/diplomka/protein-binding-sites/p2rank/p2rank-llmex-data/OUTPUT/_$ligand/analyze_labeled-residues_"$ligand"_"$type"/run.log
    not_working=$(grep -o '\[ERROR\] Dataset - error processing dataset item \[\([^]]*\)\]' $file | sed 's/\[ERROR\] Dataset - error processing dataset item \[\([^]]*\)\]/\1/' | tr '\n' ' ')

    echo not working: "$not_working"

    file=/home/brabecm4/diplomka/protein-binding-sites/p2rank/p2rank-llmex-data/_$ligand/"$ligand"_"$type".ds
    new_file=/home/brabecm4/diplomka/protein-binding-sites/p2rank/p2rank-llmex-data/_$ligand/"$ligand"_"$type"_clean.ds

    echo > $new_file

    while IFS= read -r line; do
        filename=$(echo "$line" | awk '{print $1}')
        pdb_id=$(basename "$filename")

        if [ -z "$pdb_id" ]; then
            echo "$line" >> $new_file
            continue
        fi

        if [[ "$not_working" == *"$pdb_id"* ]]; then
            echo "# " "$line" >> $new_file
        else
            echo "$line" >> $new_file
        fi
    done < "$file"
done

# train eval
cd /home/brabecm4/diplomka/protein-binding-sites/p2rank/p2rank-llmex-data
../p2rank/prank.sh -c residues.groovy -out_subdir _$ligand traineval -t _$ligand/"$ligand"_train_clean.ds -e _$ligand/"$ligand"_test_clean.ds 