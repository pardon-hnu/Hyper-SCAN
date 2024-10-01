# Hyper-SCAN
This respository is created for the paper: Edge-centered Structural Graph Clustering Over Hypergraph

## datasets
enron and southern-women are shared in the datasets floder.
Entire datasets can be downloaded at https://drive.google.com/file/d/1Si69amFBcO3kpHxfp0f1M-TLYkY1aySb/view?usp=sharing.

## floder structure
```
----Hyper-SCAN
--------code: source code of pSCAN-adp/GS*-index/LSBI
--------datasets: example hypergraphs
--------SCAN: source code of SCAN for pairwise graph
--------tools: tools such as getting basic information about hypergraphs, compute the modularity of clustering results, etc.
```

## configure
```shell
g++ -std=c++11 -Wall -g -o main main.cpp hypergraph.cpp Utility.cpp HashSet.cpp algorithms.cpp gs_index.cpp compress_index.cpp utils.h readhg.h
```

## datasets map: name,filename
```python
dataset_2_filename = {
            {"enron" , "Enron.hyp"},
            {"congress" , "congress-bills.hyp"},
            {"contact" , "contact-primary-school.hyp"},
            {"dblp", "DBLP.hyp"},
            {"aminer","aminer.hyp"},          
            {"drug","NDC-substances.hyp"},
            {"ubuntu","threads-ask-ubuntu.hyp"},
            {"rpah","rpa_t130000000.hyp"}
        }
```

## excute
```shell
./main 1 $dataset $algorithm $iteration $log $mu $epsilon
```
dataset: enron, congress, contact,dblp,aminer,drug,ubuntu,rpah          
algorithm: pSCAN-adp+LI, GS*-index+LI, LSBI
iteration: Iterations to run each algorithm on each dataset
log: Activate logging to output core-numbers & iteration h-index statistics ?
mu: parameter mu
epsilon: parameter epsilon

An example:
```shell
./main 1 enron LSBI  1 1 5 0.6
1
enron
../datasets/Enron.hyp
hypergraph ready!
enron LSBI 1 1 5 0.6
LSBI(CI+SI+LI) 
start!
compute ok!
neicun 29870420
writing to: ../my_cpp12_hash/output_lab1_csv/cluster_optimize_compress_IIII_enron_5_0.600000.csv
50
Execution time= 0.002263: init_tm= 5.17134
insert time  0.009005
remove time  0.197032
```

