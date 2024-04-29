#!/bin/bash

desired_workers=4
wait_time_mins=5
wait_time_after_worker_launched_s=20

wait_time_s=$((wait_time_mins * 60))

while true; do
    current_waint_time=$wait_time_s

    echo
    date +"[%d/%m/%y %H:%M:%S]"

    ssh gpulab 'chmod 777 /home/brabecm4/diplomka/protein-binding-sites/run_scripts/tasks_stats.sh'
    ssh gpulab 'chmod 777 /home/brabecm4/diplomka/protein-binding-sites/data/run_data/run_cases'

    echo assesing work...
    stats=$(ssh gpulab 'sh /home/brabecm4/diplomka/protein-binding-sites/run_scripts/tasks_stats.sh /home/brabecm4/diplomka/protein-binding-sites/data/run_data/run_cases' | tr -d '\r')
    echo "$stats"

    not_run_count=$(ssh gpulab 'cat /home/brabecm4/diplomka/protein-binding-sites/data/run_data/run_cases | grep -c "NOT_RUN"')

    if ((not_run_count == 0)); then
        echo work is done
        echo ...exiting
        exit 0
    fi

    echo checking workers ...

    workers=$(ssh gpulab 'squeue | grep -c "brabecm4"')
    missing_workers=$((desired_workers - workers))

    echo $workers alive, $missing_workers are missing

    while [ $missing_workers -gt 0 ]; do
        echo \ \ \ launching new worker

        ssh gpulab 'sbatch /home/brabecm4/diplomka/protein-binding-sites/run_scripts/launch_worker.sh /home/brabecm4/diplomka/protein-binding-sites/data/run_data/run_cases'

        ((missing_workers--))
        current_waint_time=$wait_time_after_worker_launched_s
        sleep 5
    done

    echo sleeping for $current_waint_time seconds
    sleep $current_waint_time
done