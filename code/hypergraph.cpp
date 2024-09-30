#include <vector>
#include <set>
#include <map>
#include <string> 
#include <iostream>
#include <set>
#include <algorithm>
#include <limits>
#include "hypergraph.h"

// typedef  std::map<std::string, size_t> strIntMap;
// typedef  std::map<std::string, std::vector<size_t>> strvIntMap;
// typedef  std::map<std::string, std::vector<std::string>> strvStrMap;
// typedef  std::map<size_t, std::vector<std::string>> intvStrMap;
// typedef  std::vector<std::string> strvec;
// typedef  std::vector<size_t> intvec;


    // std::map<size_t, strvec > e_id_to_edge; // # key => hyperedge_id, value => List of vertices in a hyperedge
   

    // public:
    //  // Auxiliary variables
    // std::map<std::string, std::vector<size_t> > inc_dict;  //# key => node, value = List of incident hyperedge ids.
    // std::map<std::string, strvec > init_nbr;  //# key => node id, value => List of Neighbours.
    // strIntMap init_nbrsize; // # initial nbrhood sizes. 
    // strvec init_nodes;
    // /* data structures that make the algorithm more efficient */
    // std::map<size_t, size_t > edge_min_hindex;
    // strIntMap lub;
    // strIntMap llb;
    Hypergraph::Hypergraph(){}
    Hypergraph::~Hypergraph(){}

    void Hypergraph::addEdge(size_t id, strvec edge){//每次只添加一条超边id
        // e_id_to_edge[id] = std::move(edge);
        intvec e(edge.size());
        for(size_t i = 0; i<edge.size(); i++) {
            sscanf(edge[i].c_str(), "%zu", &e[i]);//超边里面包含的顶点开始是按照字符串来存储其编号的，这里将字符串转换成无符号整数
        }  
        hyperedges.push_back(e);   //可以按照下标来访问每条超边，这个下标视为了超边的编号，而没有使用这里的id这一参数值？??
        //好像不是？这样的话e_id_to_edge是怎么知道边和编号对应关系的？  用的下标
        
    }
    void Hypergraph::initialise(){ //<顶点,编号>键值对 存在map node_index里面，init_nodes存储了所有顶点
        /*
        Initialise different variables e.g. number of neighbours, neighbour list etc.
        */
       std::unordered_set<size_t> V;
        for(auto elem: hyperedges){
            for(auto v_id: elem){
                // initalise init_nodes and incident dictionary
                V.insert(v_id);
            }
        }
        size_t i = 0;
        for(auto v: V)  {
            init_nodes.push_back(v);
            node_index[v] = i; //key和value值有相同不会影响，因为key和value各自都没有重复，在map里面按照key找value也不会找到不相干的定点去
            i+=1;
        }

    }
    // void Hypergraph::init_nbrs(){
    //     for(auto elem: e_id_to_edge){
    //         for(auto v_id: elem.second){
    //         // initialise number of neighbours and set of neighbours
    //             if ( init_nbr.find(v_id) == init_nbr.end() ) { // first insertion of v_id to init_nbr map
    //                 init_nbr[v_id] = strset();
    //                 for (auto u: elem.second){
    //                     if (u!=v_id){
    //                         init_nbr[v_id].insert(u);
    //                     }
    //                 }
    //                 init_nbrsize[v_id] =  init_nbr[v_id].size();
    //             }
    //             else{  // v_id exists in init_nbr map
    //                 for (auto u: elem.second){
    //                     if (u!=v_id){
    //                         init_nbr[v_id].insert(u);
    //                     }
    //                 }
    //                 init_nbrsize[v_id] = init_nbr[v_id].size();
    //             }
    //         }
    //     }   
    // }


    void Hypergraph::printHypergraph(){
        /* 
        Prints the edge list 
        */

        for (auto elem: hyperedges){
            // std::cout <<elem.first<<":";
            for(auto u: elem){
                std::cout << u <<" ";
            }
            std::cout<<"\n";
        }
      
    }
    void Hypergraph:: writeneighborhood(std::string file){
        if(file=="") file = "../output/log_"+dataset+"Nv.csv";
        std::cout <<"writing neighborhood dictionary: "<< file<<"\n";
        std::stringstream ss;
        std::unordered_map< size_t, std::unordered_set<size_t>> init_nbr;
        /* Compute neighbors first */
        for(auto elem: hyperedges){
            for(auto v_id: elem){
            // initialise number of neighbours and set of neighbours
                if ( init_nbr.find(v_id) == init_nbr.end() ) { // first insertion of v_id to init_nbr map
                    init_nbr[v_id] = std::unordered_set<size_t>();
                    for (auto u: elem){
                        if (u!=v_id){
                            init_nbr[v_id].insert(u);
                        }
                    }
                }
                else{  // v_id exists in init_nbr map
                    for (auto u: elem){
                        if (u!=v_id){
                            init_nbr[v_id].insert(u);
                        }
                    }
                }
            }
        }  
        for(auto pair: init_nbr){
            auto node = pair.first;
            ss<<node<<",";
            int _count = 0;
            int N = pair.second.size();
            for (auto nbr_v: pair.second){
                _count++;
                if(_count<N) ss<<nbr_v<<",";
                else    ss<<nbr_v<<"\n";
            }
        }
        std::ofstream out(file.c_str());
        if(out.fail())
        {
            out.close();
        }
        out << ss.str();
        out.close();
    }