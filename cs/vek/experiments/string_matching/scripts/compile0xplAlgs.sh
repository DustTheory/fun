#!/bin/bash
scriptDir=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
cd $scriptDir
g++ -O3 ../src/0x80plAlgo1.cpp -o ../build/0x80plAlgo1 
g++ -O3 -msse4 ../src/0x80plAlgo2.cpp -o ../build/0x80plAlgo2
g++ -O3 -msse4.2 ../src/0x80plAlgo3.cpp -o ../build/0x80plAlgo3 