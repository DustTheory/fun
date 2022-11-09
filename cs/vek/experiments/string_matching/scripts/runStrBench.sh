#!/bin/bash
scriptDir=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
cd $scriptDir
gcc "$scriptDir/../bench/str-bench.c" "$scriptDir/../build/stringmatching.o" -o "$scriptDir/../build/str-bench.o" && \ 
$scriptDir/../build/str-bench.o