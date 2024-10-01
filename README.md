# Hyper-SCAN
This respository is created for the paper: Edge-centered Structural Graph Clustering Over Hypergraph

## datasets
enron and southern-women are shared in the datasets file.
Entire datasets can be downloaded at https://drive.google.com/file/d/1Si69amFBcO3kpHxfp0f1M-TLYkY1aySb/view?usp=sharing.

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
            {"threads-ask-ubuntu","threads-ask-ubuntu.hyp"},
            {"rpah","rpa_t130000000.hyp"}
        }
```

# waitting for upload main.cpp