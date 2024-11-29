#include <vector>
#include <set>
#include <map>
#include <string> 
#include <iostream>
#include <cmath>
#include <set>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <tuple>
#include "hypergraph.h"
#include "gs_index.h"
#include "HashSet.h"

void Gs_index::writecluster(std::string folder,std::vector<unsigned int> edge_to_cluster,double sim,unsigned int u){ 
    std::stringstream ss;
    for(size_t i = 0; i <hg.hyperedges.size(); i++){
        //for(auto elem: Cluster[i]){
            ss<<std::to_string(i)<<' '<<std::to_string(edge_to_cluster[i]);
       // }
        ss<<"\n";
    }
    std::string sim1=std::to_string(sim);
    std::string u1=std::to_string(u+1);
    std::string lujin="../output/";
    std::string file1=lujin+"cluster_"+output["algo"]+"_"+hg.dataset+"_"+u1+"_"+sim1+".csv";
    std::cout<<"writing to: "<<file1<<"\n";
    std::ofstream out(file1.c_str()); 
    if(out.fail())
    {
        out.close();
    }
    out << ss.str();
    out.close();
}


Gs_index::Gs_index(Hypergraph &H){
    hg = H;
    output["dataset"] = H.dataset;
}
Gs_index::~Gs_index(){}

bool comp(core_pair &a,core_pair &b)
{
    if(a.xiangsidu!=b.xiangsidu) return a.xiangsidu<b.xiangsidu;
    else return a.eid<b.eid;
}
 void Gs_index::write_others(size_t Cluster_num)
 {
        this->Cluster_num=Cluster_num;
 }

 void Gs_index::write_index(std::vector<std::vector<core_pair>> core_order){

    for(size_t i=0;i<core_order.size();i++)
    {
        (this->core_order).push_back(core_order[i]);
    }

 }


void generate_neighbors_by_tables(size_t i,intintvec &e_id_to_edge,intIntMap& node_index,intintvec &tables,std::set<size_t> &neighbors)
{
        for(auto v:e_id_to_edge[i])
        {
            auto j = node_index[v];
            for(size_t k=0;k<tables[j].size();k++)
            {
                if(tables[j][k]==i) continue;
                else neighbors.insert(tables[j][k]);
            }

        }
}

void generate_neighbors_by_enumerate(size_t i,intintvec &e_id_to_edge,std::set<size_t> &neighbors)
{
    for(size_t j=0;j<e_id_to_edge.size();j++)
    {
        if(j==i) continue;
        for (unsigned int k=0;k<e_id_to_edge[j].size();k++) {
                auto it = std::find(e_id_to_edge[i].begin(), e_id_to_edge[i].end(), e_id_to_edge[j][k]);
                if (it != e_id_to_edge[i].end()) {
                    neighbors.insert(j);
                    break;
                }
        }
    }

}



double geometricMean_gs1(std::vector<HashSet> &hash,const intintvec & Hyperedge,unsigned int num1,unsigned int num2,size_t e1,size_t e2) {        /**********得到相似度*/
    
    double result = sqrt(static_cast<double>(num1) * num2);
    double divisionResult = 0;  
    unsigned int count=0;

    size_t temp=0;
    if(num1>num2) {temp=e1;e1=e2;e2=temp;}
    
    for (unsigned int j=0;j<Hyperedge[e1].size();j++) {
              if (hash[e2].find(Hyperedge[e1][j])==1) {
                count++;
              }    
    }
    divisionResult = static_cast<double>(count) / result;
    return divisionResult;
}

double geometricMean_gs2(intIntMap& node_index,std::vector<short> &hash1,const intintvec & Hyperedge,unsigned int num1,unsigned int num2,size_t e1,size_t e2) {        /**********得到相似度*/
    
    double result = sqrt(static_cast<double>(num1) * num2);
    double divisionResult = 0;  
    unsigned int count=0;
    for (unsigned int j=0;j<Hyperedge[e1].size();j++) {
              auto i=node_index[Hyperedge[e1][j]];
              if (hash1[i]==1) {
                count++;
              } 
    }
    divisionResult = static_cast<double>(count) / result;
    return divisionResult;
}
double geometricMean_gs(const intintvec & Hyperedge,unsigned int num1,unsigned int num2,size_t e1,size_t e2) { 
    
    double result = sqrt(static_cast<double>(num1) * num2);
    double divisionResult = 0;  
   std::set<size_t> gongtongdingdian;
    for (unsigned int j=0;j<Hyperedge[e1].size();j++) {
              auto it = std::find(Hyperedge[e2].begin(), Hyperedge[e2].end(), Hyperedge[e1][j]);
              if (it != Hyperedge[e2].end()) {
                gongtongdingdian.insert(Hyperedge[e1][j]);
              }    
    }
    divisionResult = static_cast<double>(gongtongdingdian.size()) / result;
    return divisionResult;
}

void gs_index_construct(std::string dataset, intintvec& e_id_to_edge,intvec& init_nodes, intIntMap& node_index,Gs_index& a)
{
    
   
    clock_t start,e_tm;
    start = clock();
    std::vector<core_pair> Neighbor;
    size_t N = init_nodes.size();
    intintvec tables(N,intvec{});
    std::set<size_t> computed_neighbor;
    for(size_t i=0;i<e_id_to_edge.size();i++)
    {
        for (const auto& elem:e_id_to_edge[i]) {
             auto j = node_index[elem];
             tables[j].push_back(i);
        }

    }
    e_tm = clock();
    a.LI_time=double(e_tm - start) / double(CLOCKS_PER_SEC);
    //std::cout<<a.LI_time<<'\n';

    double xiangsidu=0;
    core_pair cp;
    for(size_t i=0;i<e_id_to_edge.size();i++)
    {

        generate_neighbors_by_tables(i,e_id_to_edge,node_index,tables,computed_neighbor);
        //generate_neighbors_by_enumerate(i,e_id_to_edge,computed_neighbor);
        for(auto v:computed_neighbor)
        {
           xiangsidu=geometricMean_gs(e_id_to_edge,e_id_to_edge[i].size(),e_id_to_edge[v].size(),i,v);  
            if(i!=v){
                cp.eid=v;
                cp.xiangsidu=-xiangsidu;
                Neighbor.push_back(cp);
            }
        }

        sort(Neighbor.begin(),Neighbor.end(),comp);
        (a.neighbor_order).push_back(Neighbor);
        Neighbor.clear();
        computed_neighbor.clear();
    }

    std::cout<<"compute ok!"<<'\n';
    size_t max_u=0;
    for(size_t i=0;i<(a.neighbor_order).size();i++)
    {
        if(a.neighbor_order[i].size()>max_u) max_u=a.neighbor_order[i].size();
    }
    std::vector<std::vector<core_pair>> core_order(max_u,std::vector<core_pair>{});
    
    for(size_t j=0;j<(a.neighbor_order).size();j++)
    {
        unsigned int i=1;
        for(auto x: a.neighbor_order[j])
        {
          cp.eid=j;
          cp.xiangsidu=x.xiangsidu;
          core_order[i-1].push_back(cp); 
          ++i;
        }
    }
    for(size_t j=0;j<max_u;j++)
    {
        sort(core_order[j].begin(),core_order[j].end(),comp);
    }
    std::cout<<"core-order ok!"<<'\n';
    a.write_index(core_order);
    e_tm = clock();
    a.init_time=double(e_tm - start) / double(CLOCKS_PER_SEC);
    a.output["inition time"]= std::to_string(a.init_time);
   
    a.max_neighbor=max_u;
    size_t neicun=0;
    for(size_t i=0;i<(a.neighbor_order).size();i++) {neicun+=a.neighbor_order[i].size()*16;}
    for(size_t i=0;i<core_order.size();i++) {neicun+=core_order[i].size()*16;}
    a.neicun=neicun;
    size_t LI_neicun=0;
    for(size_t i=0;i<tables.size();i++) {LI_neicun+=tables[i].size()*8;}
    a.LI_neicun=LI_neicun;
    std::cout<<"neicun:  "<<neicun<<"(Bytes)"<<'\n';
    
}

void gs_index_cluster(intintvec& e_id_to_edge,Gs_index& a,double sim,unsigned int u)
{
    a.output["algo"] = "SQuery-OI";
    clock_t start, end;
    start = clock();
    std::vector<unsigned int> edge_to_cluster(e_id_to_edge.size(),UINT_MAX);
    
    intvec cores;
    for(const auto &elem:a.core_order[u-1]){
        if(-elem.xiangsidu<sim) break;  
        cores.push_back(elem.eid);
    }

    sort(cores.begin(),cores.end());

    size_t Cluster_ID=1;
    for(const auto &elem:cores){
        if(edge_to_cluster[elem]!=UINT_MAX) continue;
        edge_to_cluster[elem]=elem;
        std::queue<size_t> Q;
        Q.push(elem);
        while(Q.empty()!=true){
            size_t temp=Q.front();
            Q.pop();
            for ( auto elem1:a.neighbor_order[temp]) {
                if(-elem1.xiangsidu<sim) break;  
                if(edge_to_cluster[elem1.eid]==UINT_MAX){
                    edge_to_cluster[elem1.eid]=elem;
                    if(a.neighbor_order[elem1.eid].size()>=u){ 
                        if(-(a.neighbor_order[elem1.eid][u-1].xiangsidu)>=sim) Q.push(elem1.eid);
                    }
                }
       
            }
        }
        Cluster_ID++;  
    }
    end = clock();
    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.write_others(Cluster_ID-1);
    a.writecluster("../output/",edge_to_cluster,sim,u);
   // std::cout<<Cluster_ID-1<<'\n';
}