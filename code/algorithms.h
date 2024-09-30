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


//for compress
typedef std::tuple<size_t,size_t,size_t,double,double> compresstuple;// 起始编号 结束编号 个数  相似度下界 相似度上界
typedef std::vector<compresstuple> compresstuplevec;
typedef std::vector<compresstuplevec> compresstuplevecvec;


//const size_t MAX_SIZE=12345678;
//intintvec Cluster(MAX_SIZE,intvec{});
class Algorithm{
    Hypergraph hg;//实际上需要使用到node_index时直接(a.hg).node_index就可以  且在gs-index类里面  维护索引的时候可能有用
    public:
    
    //intvec Edge_Cluster;//记录每条边属于哪个cluster
   // intintvec Cluster(MAX_SIZE,intvec{});  //纪录每个Cluster  , intintvec Cluster(size,intvec{})  ?????
    size_t Cluster_num=0;
    //std::vector <intvec> cluster(MAX_SIZE,intvec{});
    double exec_time = 0;
    //double core_exec_time = 0;
    //double correction_time = 0;
    //size_t nu_cu = 0;
   // size_t num_nbr_queries = 0;
    strstrMap output;
    //std::vector< strstrMap > hnlog;
    //strstrMap timelogs;
    Algorithm( Hypergraph &H);
    ~Algorithm();
    // void Peel(bool verbose = false);
    // void EPeel(bool verbose = false);
    // void local_core(bool log = false);
    // void local_core_opt_core_correct(bool log = false);
    // void local_core_omp(std::map<size_t, strvec >&, bool log = false);
    // size_t iterative_core_correct(Hypergraph& H, std::string u, size_t core_u, strIntMap &hn);
    // size_t iterative_core_correct_opt(Hypergraph& H, std::string u, size_t core_u);
   // void printcore();
   // bool write_results();
    //void write_to_Cluster(size_t i,size_t b);
    void write_to_Cluster_num(size_t b);
    void writecluster(std::string folder,intintvec Cluster,double sim,unsigned int u);
   // void writelog();
   // void writeNbrQ();
   //void writekdcore(std::string folder="../output/");
};
   void base_line(std::string dataset, intintvec e_id_to_edge, Algorithm& a,double sim,unsigned int u);
   
   void optimize_index(std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,Algorithm& a,double sim,unsigned int u);
  
   void optimize_indexII(std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,Algorithm& a,double sim,unsigned int u);
   
   void optimize_compress(std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,Algorithm& a,double sim,unsigned int u);

#endif