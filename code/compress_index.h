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
//const size_t MAX_SIZE=12345678;
//intintvec Cluster(MAX_SIZE,intvec{});
class compress_index{
    public:
    Hypergraph hg;//实际上需要使用到node_index时直接(a.hg).node_index就可以  且在gs-index类里面  维护索引的时候可能有用
    intintvec tables;
    std::vector<intintvec> compress;
    std::vector<std::vector<core_pair2>> Cluster_Index;
    std::vector<std::vector<unsigned int>> segment_count;
    //intvec Edge_Cluster;//记录每条边属于哪个cluster
   // intintvec Cluster(MAX_SIZE,intvec{});  //纪录每个Cluster  , intintvec Cluster(size,intvec{})  ?????
    size_t Cluster_num=0;
    //std::vector <intvec> cluster(MAX_SIZE,intvec{});
    double exec_time = 0;
    double init_time=0;
    double updata_insert_time=0;
    double updata_remove_time=0;
    size_t neicun=0;
    //double core_exec_time = 0;
    //double correction_time = 0;
    //size_t nu_cu = 0;
   // size_t num_nbr_queries = 0;
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
    void writecluster(std::string folder,intintvec Cluster,double sim,unsigned int u);
   // void writelog();
   // void writeNbrQ();
   //void writekdcore(std::string folder="../output/");
};
   void compress_II_construct(double d,std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u);
   void compress_II_cluster(double d,std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u);

   void compress_III_construct(double d,double d_segment_count,std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u);
   void compress_III_cluster(double d,double d_segment_count,std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u);
   
   void compress_IIII_construct(double d,double d_segment_count,std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u);
   void compress_IIII_cluster(double d,double d_segment_count,std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u);
   void compress_IIII_cluster_bianjie(double d,double d_segment_count,std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u);
  void compress_IIII_cluster_hash(double d,double d_segment_count,std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u);
   void compress_IIII_cluster_bianjie_hash(double d,double d_segment_count,std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u);
   void updata_insert( compress_index &a,size_t node,size_t num_hyperegde);
   void updata_remove( compress_index &a,size_t num_node,size_t num_hyperegde);

#endif