# Hyper-SCAN
This respository is created for the paper: Edge-centered Structural Graph Clustering Over Hypergraph

## configure
```shell
g++ -std=c++11 -Wall -g -o main main.cpp hypergraph.cpp Utility.cpp HashSet.cpp algorithms.cpp gs_index.cpp compress_index.cpp utils.h readhg.h
```

## datasets map: name,filename
readhg.h:
```C++
std::map <std::string,std::string> dataset_to_filename = {
            {"enron" , "../datasets/Enron.hyp"},
            {"congress" , "../datasets/congress-bills.hyp"},
            {"contact" , "../datasets/contact-primary-school.hyp"},
            {"dblp", "../datasets/DBLP.hyp"},
            {"aminer","../datasets/aminer.hyp"},          
            {"drug","../datasets/NDC-substances.hyp"},
            {"threads-ask-ubuntu","../datasets/threads-ask-ubuntu.hyp"},
            {"rpah","../datasets/rpa_t130000000.hyp"}
        };
```

# waitting for upload main.cpp