#!/bin/bash
scriptDir=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
cd $scriptDir
nasm -O3 -f elf64 -o ../build/stringmatching_asm.o ../src/stringmatching_asm.asm && \
g++ -c -O3 -no-pie ../src/stringmatching.c ../build/stringmatching_asm.o -o ../build/stringmatching && \
../build/stringmatching 