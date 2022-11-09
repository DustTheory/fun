#!/bin/bash
scriptDir=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
python3 "$scriptDir/../bench/0xplbenchmark.py" "$scriptDir/../test_data/portugal.txt"