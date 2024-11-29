# Hyper-SCAN
This respository is created for the paper: Efficient Structural Clustering over Hypergraphs (submitted to ICDE 2025 round 2)

## datasets
po, enron and southern-women are shared in the datasets floder.
Entire datasets can be downloaded at https://drive.google.com/file/d/1Si69amFBcO3kpHxfp0f1M-TLYkY1aySb/view?usp=sharing.

## floder structure
```
----Hyper-SCAN
--------code: source code of pSCAN-adp/GS*-index/LSBI
--------datasets: example hypergraphs
--------SCAN: source code of SCAN for pairwise graph
--------tools: tools such as getting basic information about hypergraphs, compute the modularity of clustering results, etc.
--------output: the clustering results
```

## configure
```shell
g++ -std=c++17 -pthread -I./parlaylib/include -Wall -g -o main main.cpp hypergraph.cpp Utility.cpp HashSet.cpp algorithms.cpp gs_index.cpp compress_index.cpp readhg.h
```

## datasets map: name,filename
```python
dataset_2_filename = {
            {"po,""po.hyp"},
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
./main $dataset $algorithm $mu $epsilon
```
dataset: enron, congress, contact,dblp,aminer,drug,ubuntu,rpah          
algorithm: pHSCAN, OI-construction, LSBI-construction, SQuery,PQuery
mu: parameter mu
epsilon: parameter epsilon

An example:
```shell
./main enron SQuery 2 0.7
enron
../datasets/Enron.hyp
hypergraph ready!
enron SQuery 2 0.7
SQuery 
start!
compute ok!
neicun 29870420(Bytes)
writing to: ../output/cluster_SQuery_enron_2_0.700000.csv
Execution time= 0.00246434: init_tm= 5.01741
```

## remarks
PQuery depends on the exist lib: parlylib 
the details for compiling and running can refer to code/run.sh\
