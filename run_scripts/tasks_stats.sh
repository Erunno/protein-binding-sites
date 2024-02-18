#!/bin/bash

if [ -z "$1" ]; then
  echo "Usage: $0 <filename with commands>"
  exit 1
fi

total=$(cat $1 | wc -l)
not_run=$(grep -c "NOT_RUN" < $1)
success=$(grep -c "DONE" < $1)
running=$(grep -c "RUNNING" < $1)
errors=$(grep -c "ERROR" < $1)
finished=$((success + errors))

echo stats: \[$finished/$total\] done
echo \ \ running=$running
echo \ \ not_run=$not_run 
echo \ \ success=$success
echo \ \ errors=$errors

