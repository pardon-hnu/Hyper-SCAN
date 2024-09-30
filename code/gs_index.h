#ifndef GS_INDEX_H
#define GS_INDEX_H
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


//for gs-index
typedef std::vector<std::pair<double,size_t>> intdouprvec;
typedef std::vector<intdouprvec> intintdouprvec;


//for gs-index  set   自动排序 默认按第一位   保证无重复？
typedef std::set<std::pair<double,size_t>>  douintset;      //单个core-order   单个neighbor-order
typedef std::vector<douintset> douintsetvec;            //neighbor order     多个core-order
typedef std::pair<double,size_t> douintpair;

struct core_pair
{
    /* data */
    double xiangsidu;
    size_t eid;
};

//const size_t MAX_SIZE=12345678;
//intintvec Cluster(MAX_SIZE,intvec{});
class Gs_index{
    public:
    Hypergraph hg;//实际上需要使用到node_index时直接(a.hg).node_index就可以  且在gs-index类里面  维护索引的时候可能有用
    std::vector<std::vector<core_pair>> core_order;
    std::vector<std::vector<core_pair>> neighbor_order;
    //intvec Edge_Cluster;//记录每条边属于哪个cluster
   // intintvec Cluster(MAX_SIZE,intvec{});  //纪录每个Cluster  , intintvec Cluster(size,intvec{})  ?????
    size_t Cluster_num=0;
    size_t max_neighbor=0;
    //std::vector <intvec> cluster(MAX_SIZE,intvec{});
    double exec_time = 0;
    double init_time=0;
    double LI_time=0;
    double insert_time=0;
    double remove_time=0;
    size_t neicun=0;
    size_t LI_neicun=0;
    strstrMap output;
    //std::vector< strstrMap > hnlog;
    //strstrMap timelogs;
    Gs_index( Hypergraph &H);
    ~Gs_index();
  
    void write_index(std::vector<std::vector<core_pair>> core_order);// 写入cluster_num  exectime inittime core_order neighbor_order max_neighbor等信息
    void write_others(size_t Cluster_num);
    void writecluster(std::string folder,intintvec Cluster,double sim,unsigned int u);
   
};

   void gs_index_construct(std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,Gs_index& a,double sim,unsigned int u);
   void gs_index_cluster(intintvec e_id_to_edge,Gs_index& a,double sim,unsigned int u);
   void gs_index_insert( Gs_index &a,intvec insert_hyperedge);
   void gs_index_remove(Gs_index &a,size_t eid);
   
#endif