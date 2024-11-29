#!/bin/bash

g++ -std=c++17 -pthread -I./parlaylib/include -Wall -g -o main main.cpp hypergraph.cpp Utility.cpp HashSet.cpp algorithms.cpp gs_index.cpp compress_index.cpp readhg.h #utils.h  run on gnu c++ compiler

declare -a dset=("po" "enron" "contact" "congress" "drug"  "ubuntu" "dblp" "aminer" "rpah") 
  
declare -a algorithms=("OI-construction" "LSBI-construction" "pHSACN" "SQuery" "PQuery") 

declare -a epision=(0.2 0.3 0.4 0.5 0.6 0.7 0.8)

declare -a miu=(2 5 10 15)

for dataset in "${dset[@]}"
do
    for algo in "${algorithms[@]}"
    do
        for mi in "${miu[@]}"
        do
            for ep in "${epision[@]}"
            do
                ./main $dataset $algo $mi $ep
                echo "------------"
            done
        done
    done 
done
