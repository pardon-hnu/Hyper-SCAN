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
#include "algorithms.h"
//#include "utils.h"


void Algorithm::writecluster(std::string folder,std::vector<unsigned int> edge_to_cluster,double sim,unsigned int u){
    //std::string file = folder + "cluster_"+output["algo"]+"_"+hg.dataset+".csv";
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


Algorithm::Algorithm(Hypergraph &H){
    hg = H;
    output["dataset"] = H.dataset;
}
Algorithm::~Algorithm(){}


 void Algorithm::write_to_Cluster_num(size_t b)
 {
        Cluster_num=b;
 }


bool geometricMean2(double sim,intIntMap& node_index,const intintvec & Hyperedge,unsigned int num1,unsigned int num2,size_t e1,size_t e2) {        /**********得到相似度*/
    
    double result = sqrt(static_cast<double>(num1) * num2);
    double divisionResult = static_cast<double>(num1>num2?num2:num1) / result; 

    if(divisionResult<sim) return false;
    unsigned int count=0;
    divisionResult=0;
    std::set<size_t> gongtongdingdian;
    for (unsigned int j=0;j<Hyperedge[e1].size();j++) {
              auto it = std::find(Hyperedge[e2].begin(), Hyperedge[e2].end(), Hyperedge[e1][j]);
              if (it != Hyperedge[e2].end()) {
                gongtongdingdian.insert(Hyperedge[e1][j]);
                divisionResult = static_cast<double>(gongtongdingdian.size()) / result;
                if(divisionResult>=sim) break;   
              }    
    }
    if(divisionResult>=sim) return true;
    else return false;
}
bool checkcore_optimize_index(double sim,unsigned int u,size_t i,intintvec & tables,intIntMap& node_index,const intintvec & Hyperedge)
{
        std::set<size_t> neighbor;
        size_t count=0;
        for(auto v:Hyperedge[i]){
             auto j = node_index[v];
             for(size_t k=0;k<tables[j].size();k++)
             {
                        if(i==tables[j][k]) continue;
                        auto it = neighbor.find(tables[j][k]);
                        if (it != neighbor.end()) continue;
                        else{neighbor.insert(tables[j][k]);}

                        if(geometricMean2(sim,node_index,Hyperedge,Hyperedge[i].size(),Hyperedge[tables[j][k]].size(),i,tables[j][k]))
                         {
                            count+=1;
                         }
                          if(count>=u) return true;
             }
        }
    if(count>=u) return true;
    else return false;
}
void compute_optimize_index_n_sim(double sim,size_t i,intintvec & tables,intvec& n_sim,intIntMap& node_index,const intintvec & Hyperedge)
{

        std::set<size_t> neighbor;
        for(auto v:Hyperedge[i]){
             auto j = node_index[v];
             for(size_t k=0;k<tables[j].size();k++)
             {
                        if(i==tables[j][k]) continue;
                        auto it = neighbor.find(tables[j][k]);
                        if (it != neighbor.end()) continue;
                        else{neighbor.insert(tables[j][k]);}

                        if(geometricMean2(sim,node_index,Hyperedge,Hyperedge[i].size(),Hyperedge[tables[j][k]].size(),i,tables[j][k]))
                         {
                            n_sim.push_back(tables[j][k]);
                         }
                        
             }
        }
}
void optimize_index(std::string dataset, intintvec& e_id_to_edge,intvec& init_nodes, intIntMap& node_index,Algorithm& a,double sim,unsigned int u)
{
    
    a.output["algo"] = "pHSCAN";
    clock_t start, end;
    start = clock();
    size_t N = init_nodes.size();
    intintvec tables(N,intvec{});
    for(size_t i=0;i<e_id_to_edge.size();i++)
    {
        for (const auto& elem:e_id_to_edge[i]) {
             auto j = node_index[elem];
             tables[j].push_back(i);
        }
    }

    clock_t e_tm = clock();
    a.output["init_time"] = std::to_string(double(e_tm - start) / double(CLOCKS_PER_SEC));
   
    std::vector<unsigned int> edge_to_cluster(e_id_to_edge.size(),UINT_MAX);
    intvec n_sim;
    size_t Cluster_ID=1;
    for(size_t eid = 0; eid < e_id_to_edge.size(); eid++){
         if(checkcore_optimize_index(sim,u,eid,tables,node_index,e_id_to_edge)&&edge_to_cluster[eid]==UINT_MAX){
            edge_to_cluster[eid]=eid;
            std::queue<size_t> Q;
            Q.push(eid);
            while(Q.empty()!=true){
                size_t temp=Q.front();
                Q.pop();
                compute_optimize_index_n_sim(sim,temp,tables,n_sim,node_index,e_id_to_edge);
                for (const auto& elem1:n_sim) {
                    if(edge_to_cluster[elem1]==UINT_MAX){
                        
                        edge_to_cluster[elem1]=eid;

                        if(checkcore_optimize_index(sim,u,elem1,tables,node_index,e_id_to_edge)){
                            Q.push(elem1);
                        }
                    }
                }
                n_sim.clear();
            }
            Cluster_ID++;
         }
    }
    end = clock();
    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.write_to_Cluster_num(Cluster_ID-1);
    a.writecluster("../output/",edge_to_cluster,sim,u);
   // std::cout<<"cluster numbers: "<<Cluster_ID-1<<'\n';

}
