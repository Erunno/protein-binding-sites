# run_scripts

This folder contains scripts that orchestrate work on the Parlab cluster at Charles University. The scripts fall into different categories as explained in the subsequent sections.

## P2Rank

The [./p2rank](./p2rank/) folder contains scripts that run the [P2Rank tool](https://github.com/rdk/p2rank).

## Training On Multiple Parlab Nodes

To leverage the capabilities of multiple Parlab nodes provided during the thesis work, we implemented a simple framework that uses a file with multiple commands in the following format:

```csv
DONE;python /path/to/repo/netws/network.final_eval.py --final-tag nei_5_comprs_v2_c --ligand GTP
NOT_RUN;python /path/to/repo/netws/network.final_eval.py --final-tag nei_3_comprs_v2_c --ligand DNA
RUNNING;python /path/to/repo/netws/network.final_eval.py --final-tag nei_3_comprs_v2_c --ligand MN
ERROR;python /path/to/repo/netws/network.final_eval.py --final-tag nei_3_comprs_v2_c --ligand ADP
NOT_RUN;python /path/to/repo/netws/network.final_eval.py --final-tag nei_3_comprs_v2_c --ligand MG
```

Each line represents a command with an attached state: `NOT_RUN`, `RUNNING`, `DONE`, `ERROR`.

This file can be created using one of the following scripts: [`create_compare_cases.py`](./create_compare_cases.py), [`create_run_cases.py`](./create_run_cases.py), and [`create_run_final_eval.py`](./create_run_final_eval.py).

Once you have defined this file with run cases, you can pass it to the [`run_commands.sh`](./run_commands.sh) script, which executes each command sequentially. The script can be launched on Parlab using [`launch_worker.sh`](./launch_worker.sh) (`sbatch launch_worker.sh <run cases file>`). You can also launch multiple such workers all working with the same file of run cases. Workers on Parlab cannot run indefinitely. The script [`keep_workers_alive.sh`](./keep_workers_alive.sh) runs elsewhere than on Parlab and periodically checks if these workers are still working, spawning more workers if necessary.

## Additional Launchers

Many tasks benefited from having multiple workers. For instance, we proactively filled the cache using [`run_pdb_cache_renew.sh](./run_pdb_cache_renew.sh) and [`launch_pdb_cache_renew.sh](./launch_pdb_cache_renew.sh).
