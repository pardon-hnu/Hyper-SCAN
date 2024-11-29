#ifndef ALGORITHM_H
#define ALGORITHM_H
//#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <string> 
#include <iostream>
#include <sstream>
#include <fstream>
#include <set>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <tuple>
#include "hypergraph.h"
#include <ctime>
#include <climits>
typedef  std::unordered_map<size_t, size_t> intIntMap;
typedef  std::map<std::string, size_t> strIntMap;
typedef  std::map<std::string, std::vector<size_t>> strvIntMap;
typedef  std::map<std::string, std::set<size_t>> strsIntMap;
typedef  std::map<size_t, std::set<std::string>> intsStrMap;
typedef  std::map<std::string, std::vector<std::string>> strvStrMap;
typedef  std::map<size_t, std::vector<std::string>> intvStrMap;
// typedef  std::unordered_map <size_t, std::vector<size_t>> uintvIntMap;
typedef  std::vector<std::string> strvec;
typedef  std::set<std::string> strset;
typedef  std::vector<size_t> intvec;
typedef std::vector<std::pair<std::string, std::string>> strstrprvec;
typedef std::map<std::string, std::string> strstrMap;
// typedef std::map <std::string,bool> strboolMap;
typedef std::unordered_map <size_t,bool> intboolMap;
typedef std::unordered_set<size_t> uintSet;
typedef std::vector<uintSet > uintsetvec;
typedef std::unordered_map<size_t, uintSet> intuSetintMap;
//typedef std::map<std::string, std::string> strstrMap;
typedef std::vector< intvec > intintvec;//
typedef std::pair<size_t,size_t> intpair;
typedef std::tuple<size_t,size_t,size_t> inttriplet;
typedef std::vector<inttriplet> vinttriplet;



typedef std::vector<std::pair<double,size_t>> intdouprvec;
typedef std::vector<intdouprvec> intintdouprvec;



typedef std::set<std::pair<double,size_t>>  douintset;      
typedef std::vector<douintset> douintsetvec;            
typedef std::pair<double,size_t> douintpair;



typedef std::tuple<size_t,size_t,size_t,double,double> compresstuple;
typedef std::vector<compresstuple> compresstuplevec;
typedef std::vector<compresstuplevec> compresstuplevecvec;



class Algorithm{
    Hypergraph hg;
    public:

    size_t Cluster_num=0;
    double exec_time = 0;
    strstrMap output;
    Algorithm( Hypergraph &H);
    ~Algorithm();
    
    void write_to_Cluster_num(size_t b);
    void writecluster(std::string folder,std::vector<unsigned int> edge_to_cluster,double sim,unsigned int u);
  
};

   void optimize_index(std::string dataset, intintvec& e_id_to_edge,intvec& init_nodes, intIntMap& node_index,Algorithm& a,double sim,unsigned int u);
  
#endif