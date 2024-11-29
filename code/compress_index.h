#ifndef COMPRESS_INDEX_H
#define COMPRESS_INDEX_H
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
#include "gbbs/bridge.h"
#include "gbbs/helpers/sparse_table.h"
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

typedef unsigned int uintE;
#define UINT_E_MAX UINT_MAX
constexpr unsigned int kUnclustered{UINT_E_MAX};
using Clustering = parlay::sequence<uintE>;

using DirectedEdge = std::pair<uintE, uintE>;
using VertexSet =
    gbbs::sparse_table<uintE, gbbs::empty, decltype(&parlay::hash64_2)>;


struct core_pair2{
   int count;
   size_t eid;
};
struct core_pair3
{
    /* data */
    double xiangsidu;
    size_t eid;
};

class compress_index{
    public:
    Hypergraph hg;
    intintvec tables;
    std::vector<intintvec> compress;
    std::vector<std::vector<core_pair2>> Cluster_Index;
    std::vector<std::vector<unsigned int>> segment_count;
    size_t Cluster_num=0;
    double exec_time = 0;
    double init_time=0;
    double updata_insert_time=0;
    double updata_remove_time=0;
    size_t neicun=0;
    strstrMap output;
    //std::vector< strstrMap > hnlog;
    //strstrMap timelogs;
    compress_index( Hypergraph &H);
    ~compress_index();
    void write_tables(intintvec tables);
    void write_indexII(std::vector<intintvec> compress);
    void write_indexIII(std::vector<intintvec> compress,std::vector<std::vector<unsigned int>> segment_count);
    void write_indexIIII(std::vector<std::vector<core_pair2>> Cluster_Index);
    void write_others(size_t Cluster_num);
    void writecluster(std::string folder,std::vector<unsigned int> edge_to_cluster,double sim,unsigned int u);
    void writecluster_p(std::string folder,parlay::sequence<unsigned int> clustering,double sim,unsigned int u);
   
};
 
   void compress_IIII_construct(double d,double d_segment_count,std::string dataset, intintvec& e_id_to_edge,intvec& init_nodes, intIntMap& node_index,compress_index& a);
 
   void compress_IIII_cluster_bianjie(double d,double d_segment_count,std::string dataset, intintvec& e_id_to_edge,intvec& init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u);
  void compress_IIII_cluster_parallel(double d,double d_segment_count,std::string dataset, intintvec& e_id_to_edge,intvec& init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u);


#endif