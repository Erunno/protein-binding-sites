#!/bin/bash

# Check if a filename is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <filename with command>"
  exit 1
fi

filename="$1"

state=0
commands=0

# Function to handle cleanup and reset to NOT_RUN state
cleanup() {
  echo force exit ... cleaning up

  sed -i "s|RUNNING;$commands|NOT_RUN;$commands|" "$filename"
  exit 1
}

# Set up trap to call cleanup function on termination or interrupt
trap cleanup EXIT SIGTERM SIGINT SIGUSR1

while true; do

  command_found=0
  # Read the file line by line
  while IFS=';' read -r state commands; do
    # Trim leading and trailing whitespace
    state=$(echo "$state" | awk '{$1=$1};1')
    commands=$(echo "$commands" | awk '{$1=$1};1')

    # Check if the state is NOT_RUN
    if [ "$state" == "NOT_RUN" ]; then
      command_found=1
      break
    fi

  done < "$filename"

  if [ "$command_found" -eq 0 ]; then
    break
  fi

  # Change the state to RUNNING
  sed -i "s|$state;$commands|RUNNING;$commands|" "$filename"

  # Execute the commands
  eval "$commands"

  # Check the exit status of the commands
  if [ $? -eq 0 ]; then
    # If successful, change the state to DONE
    sed -i "s|RUNNING;$commands|DONE;$commands|" "$filename"
  else
    # If unsuccessful, change the state to ERROR
    sed -i "s|RUNNING;$commands|ERROR;$commands|" "$filename"
  fi

done