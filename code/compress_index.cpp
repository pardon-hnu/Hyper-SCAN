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
#include "compress_index.h"
//#include <climits>
#include "union_find_rules.h"
#include "gbbs/helpers/sparse_table.h"
#include "gbbs/bridge.h"
//#include "utils.h"
#include "HashSet.h"
#include <chrono>
#define epision 1e-10
typedef std::chrono::duration<double> tms;

void compress_index::writecluster(std::string folder,std::vector<unsigned int> edge_to_cluster,double sim,unsigned int u){
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


void compress_index::writecluster_p(std::string folder,parlay::sequence<unsigned int> clustering,double sim,unsigned int u){
    //std::string file = folder + "cluster_"+output["algo"]+"_"+hg.dataset+".csv";
    std::stringstream ss;
    for(size_t i = 0; i < hg.hyperedges.size(); i++){
       // for(auto elem: Cluster[i]){
           ss<<std::to_string(i)<<' '<<std::to_string(clustering[i]);
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

compress_index::compress_index(Hypergraph &H){
    hg = H;
    output["dataset"] = H.dataset;
}
compress_index::~compress_index(){}

void compress_index::write_tables(intintvec tables){
    for(size_t i=0;i<tables.size();i++){
        (this->tables).push_back(tables[i]);
    }
}
void compress_index::write_indexII(std::vector<intintvec> compress){
    for(size_t i=0;i<compress.size();i++){
        (this->compress).push_back(compress[i]);
    }

}
void compress_index::write_indexIII(std::vector<intintvec> compress,std::vector<std::vector<unsigned int>> segment_count){
    
    for(size_t i=0;i<compress.size();i++){
        (this->compress).push_back(compress[i]);
    }

    for(size_t i=0;i<segment_count.size();i++){
        (this->segment_count).push_back(segment_count[i]);
    }

}
void compress_index::write_indexIIII(std::vector<std::vector<core_pair2>> Cluster_Index){

    for(size_t i=0;i<Cluster_Index.size();i++){
        (this->Cluster_Index).push_back(Cluster_Index[i]);
    }
}
 void compress_index::write_others(size_t b)
 {
        Cluster_num=b;
 }

void tables_construct(intintvec &e_id_to_edge,intvec &init_nodes, intIntMap& node_index,compress_index& a)//建设完成之后，需要将tables写入到类里面
{
    size_t N = init_nodes.size();
    intintvec tables(N,intvec{});
    for(size_t i=0;i<e_id_to_edge.size();i++)
    {
        for (const auto& elem:e_id_to_edge[i]) {
             auto j = node_index[elem];
             tables[j].push_back(i);
        }
    }
    a.write_tables(tables);
}
void generate_neighbors_by_tables2(size_t i,intintvec &e_id_to_edge,intIntMap& node_index,compress_index& a,std::set<size_t> &neighbors)
{
        for(auto v:e_id_to_edge[i])
        {
            auto j = node_index[v];
            for(size_t k=0;k<a.tables[j].size();k++)
            {
                if(a.tables[j][k]==i) continue;
                else neighbors.insert(a.tables[j][k]);
            }

        }
}
void generate_neighbors_by_enumerate2(size_t i,intintvec &e_id_to_edge,std::set<size_t> &neighbors )
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

bool geometricMean2(double sim,const intintvec & Hyperedge,unsigned int num1,unsigned int num2,size_t e1,size_t e2) {       
    
    double result = sqrt(static_cast<double>(num1) * num2);
    if(result==0) return false;
    double divisionResult = static_cast<double>(num1>num2?num2:num1) / result;

    if(divisionResult<sim) return false;

    std::set<size_t> gongtongdingdian;
    divisionResult=0;
    for (unsigned int i=0;i<Hyperedge[e1].size();i++) {
              auto it = std::find(Hyperedge[e2].begin(), Hyperedge[e2].end(), Hyperedge[e1][i]);
              if (it != Hyperedge[e2].end()) {
               gongtongdingdian.insert(Hyperedge[e1][i]);
                divisionResult = static_cast<double>(gongtongdingdian.size()) / result;
                if(divisionResult>=sim) break;
            }
    }
    if(divisionResult>=sim) return true;
    else return false;
}
double geometricMean_compress_II(const intintvec & Hyperedge,unsigned int num1,unsigned int num2,size_t e1,size_t e2) {     
    
    double result = sqrt(static_cast<double>(num1) * num2);
    double divisionResult = 0;
    if(result==0) return 0;
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

bool comp2(core_pair2 &a,core_pair2 &b)
{
    if(a.count!=b.count) return a.count<b.count;
    else return a.eid<b.eid;
}
bool comp3(core_pair3 &a,core_pair3 &b)
{
    if(a.xiangsidu!=b.xiangsidu) return a.xiangsidu<b.xiangsidu;
    else return a.eid<b.eid;
}



double geometricMean_compress_IIII(std::vector<HashSet> &hash,const intintvec & Hyperedge,unsigned int num1,unsigned int num2,size_t e1,size_t e2) {        /**********得到相似度*/
    
    double result = sqrt(static_cast<double>(num1) * num2);
    double divisionResult = 0;  
    if(result==0) return 0;
    size_t temp=0;
    std::set<size_t> gongtongdingdian;
    if(num1>num2) {temp=e1;e1=e2;e2=temp;}
    
    for (unsigned int j=0;j<Hyperedge[e1].size();j++) {
              if (hash[e2].find(Hyperedge[e1][j])==1) {
               gongtongdingdian.insert(Hyperedge[e1][j]);
              }    
    }
    divisionResult = static_cast<double>(gongtongdingdian.size()) / result;
    return divisionResult;
}

void compress_IIII_construct(double d,double d_segment_count,std::string dataset, intintvec& e_id_to_edge,intvec& init_nodes, intIntMap& node_index,compress_index& a)
{
    
    clock_t start;
    start = clock();
    double d1[10]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
    
    std::vector<core_pair3> Neighbor;   
    std::set<size_t> computed_neighbor;

    unsigned int geshu1=(unsigned int)(1.0/d);
    std::vector<std::vector<size_t>> compress(geshu1,std::vector<size_t>{});
    
    unsigned int geshu=(unsigned int)(1.0/d_segment_count);
    std::vector<std::vector<core_pair2>> Cluster_Index(geshu,std::vector<core_pair2>{});

    std::vector<int> segment_count;
    std::vector< HashSet > hash;
    hash.reserve(e_id_to_edge.size());
    for(size_t i=0;i<e_id_to_edge.size();i++)
    {
        hash.emplace_back(e_id_to_edge[i].size());
    }
    
    std::cout<<"start!"<<'\n';
    tables_construct(e_id_to_edge,init_nodes,node_index,a);
    for(size_t i=0;i<e_id_to_edge.size();i++)
    {
        for (const auto& elem:e_id_to_edge[i]) {
             hash[i].insert(elem);
        }
    }
    core_pair2 cp;
    core_pair3 cp3;
    for(size_t i=0;i<e_id_to_edge.size();i++)
    {

        generate_neighbors_by_tables2(i,e_id_to_edge,node_index,a,computed_neighbor);
        for(auto v:computed_neighbor)
        {
            double xiangsidu=geometricMean_compress_IIII(hash,e_id_to_edge,e_id_to_edge[i].size(),e_id_to_edge[v].size(),i,v);
            if(i!=v){
                cp3.eid=v;
                cp3.xiangsidu=-xiangsidu;
                Neighbor.push_back(cp3);
            }
        }
        sort(Neighbor.begin(),Neighbor.end(),comp3);
        double t1=0;
        unsigned int count_1=0;

        for(auto p:Neighbor)
        {
            count_1=0;
            t1=-p.xiangsidu;
            while(t1>=d1[count_1]) 
            {
                count_1+=1;
               if(count_1>=geshu-1) break; 
            }          
            if(count_1>geshu-1) count_1=geshu-1;
            compress[geshu1-count_1-1].push_back(p.eid);
        }
       
        size_t segment_size=compress[0].size();
        for(unsigned int k=0;k<geshu;k++) {
            if(segment_size!=0) {
                cp.count=-segment_size;
                cp.eid=i;
                Cluster_Index[geshu-1-k].push_back(cp);
            }
            if(k!=geshu-1) segment_size+=compress[k+1].size();
        }
        (a.compress).push_back(compress);
        for(unsigned int k=0;k<geshu1;k++)
        {
            compress[k].clear();
        }
        computed_neighbor.clear();
        Neighbor.clear();

    }

    for(unsigned int k=0;k<geshu;k++) {
        sort(Cluster_Index[k].begin(),Cluster_Index[k].end(),comp2);
    }


    std::cout<<"compute ok!"<<'\n';
    
    a.write_indexIIII(Cluster_Index);
    clock_t e_tm = clock();
    a.init_time=double(e_tm - start) / double(CLOCKS_PER_SEC);
    a.output["inition time"]= std::to_string(a.init_time);
    size_t  neicun=0;
    for(size_t i=0;i<a.tables.size();i++) {neicun+=a.tables[i].size()*8;}
    for(size_t i=0;i<(a.compress).size();i++) {
        for(auto x:a.compress[i])
        {
            neicun+=x.size()*8;
        }
    }
    for(size_t i=0;i<Cluster_Index.size();i++) {neicun+=Cluster_Index[i].size()*12;}
    a.neicun=neicun;
    std::cout<<"neicun"<<' '<<neicun<<"(Bytes)"<<'\n';
    
}


void compress_IIII_cluster_bianjie(double d,double d_segment_count,std::string dataset, intintvec& e_id_to_edge,intvec& init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u)
{
    a.output["algo"] = "SQuery";
    
    const auto begin=std::chrono::steady_clock().now();
    std::vector<unsigned int> edge_to_cluster(a.compress.size(),kUnclustered);
    
    intvec cores;
    size_t Cluster_ID=1;
    intvec n_sim;
    double temp_sim=sim;
    unsigned int i=0;
    double d1[11]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}; 
    while (temp_sim>d1[i])
    {
        i+=1;
        if(i>=9) break;
    }

    unsigned int geshu=(unsigned int)(1.0/d_segment_count);

    for(const auto &x:a.Cluster_Index[i]){
        if(-x.count<u) break;
        cores.push_back(x.eid);
    }

   sort(cores.begin(),cores.end());

    for(const auto &x:cores){
        if(edge_to_cluster[x]!=kUnclustered) continue;
        
        edge_to_cluster[x]=x;
        std::queue<size_t> Q;
        Q.push(x);
        while(!Q.empty()){
            size_t temp=Q.front();
            Q.pop();
            for(unsigned k=0;k<geshu-i;k++)
            {
                    for ( auto elem1:a.compress[temp][k]) {
                        if(edge_to_cluster[elem1]==kUnclustered){
                            
                            edge_to_cluster[elem1]=x;
                            unsigned int count_sim=0;
                            for(unsigned int w=0;w<geshu-i;w++)
                            {
                                count_sim+=a.compress[elem1][w].size();
                                if(count_sim>=u) break;
                            }
                            if(count_sim>=u){
                                Q.push(elem1);
                            }
                        }
            
                    }
                
            }
        }
       
        Cluster_ID++;
       
    }
     const auto end=std::chrono::steady_clock().now();
    tms total_time=end-begin;
     a.exec_time = total_time.count();
    a.write_others(Cluster_ID-1);
    a.writecluster("../output/",edge_to_cluster,sim,u);
   // std::cout<<Cluster_ID-1<<'\n';

}

VertexSet MakeVertexSet(const size_t capacity) {
  return gbbs::make_sparse_table<uintE, gbbs::empty,
                                 decltype(&parlay::hash64_2)>(
      capacity, {UINT_E_MAX, gbbs::empty{}}, parlay::hash64_2);
}

void ClusterCores(
                  const intvec& cores,
                  const intintvec& core_similar_neighbors,
                  Clustering* clustering) {
  VertexSet cores_set{MakeVertexSet(cores.size())};
  parlay::parallel_for(0, cores.size(), [&](const size_t i) {
    cores_set.insert(std::make_pair(cores[i], gbbs::empty{}));
  });
  
  parlay::parallel_for(0, cores.size(), [&](const size_t i) {
    const size_t core{cores[i]};
    (*clustering)[core] = core;
  });

  constexpr auto find{gbbs::find_variants::find_compress};
  auto unite{gbbs::unite_variants::Unite<decltype(find)>{find}};

  parlay::parallel_for(0, cores.size(), [&](const size_t i) {
    const unsigned int core{cores[i]};
    const auto& neighbors{core_similar_neighbors[i]};
    constexpr bool kParallelizeInnerLoop{false};
    parlay::parallel_for(0, neighbors.size(),
                 [&](const size_t j) {
                   const unsigned int neighbor{neighbors[j]};
                   if (core > neighbor && cores_set.contains(neighbor)) {
                     unite(core, neighbor, *clustering);
                   }
                 },
                 kParallelizeInnerLoop);
  });
  parlay::parallel_for(0, cores.size(), [&](const size_t i) {
    const unsigned int core{cores[i]};
    (*clustering)[core] = find(core, *clustering);
  });

}

void AttachNoncoresToClusters(const intvec& cores,
                              const intintvec& core_similar_neighbors,
                              Clustering* clustering) {
  parlay::parallel_for(0, cores.size(), [&](const size_t i) {
    const unsigned int core{cores[i]};
    const unsigned int core_cluster{(*clustering)[core]};
    const auto& neighbors{core_similar_neighbors[i]};
    constexpr bool kParallelizeInnerLoop{false};
    parlay::parallel_for(0,neighbors.size(),
                 [&](const size_t j) {
                   const unsigned int neighbor{neighbors[j]};
                   auto* neighbor_cluster_address{&(*clustering)[neighbor]};
                   if (*neighbor_cluster_address == kUnclustered) {
                     gbbs::atomic_compare_and_swap(neighbor_cluster_address,
                                                   kUnclustered, core_cluster);
                   }
                 },
                 kParallelizeInnerLoop);
  });

}



void compress_IIII_cluster_parallel(double d,double d_segment_count,std::string dataset, intintvec& e_id_to_edge,intvec& init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u)
{
    a.output["algo"] = "PQuery";
    
    const auto begin=std::chrono::steady_clock().now();
    parlay::sequence<unsigned int> clustering(a.compress.size(),kUnclustered);//e_id_to_edge
    intvec cores;
    double temp_sim=sim;
    
    unsigned int i=0;
    unsigned geshu=(unsigned int)(1.0/d_segment_count);
    double d1[11]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
    while (temp_sim>d1[i])
    {
        i+=1;
        if(i>=9) break;
    }

    for(auto x:a.Cluster_Index[i]){ 
        if(-x.count<u) break;
        cores.push_back(x.eid);
    }

    intintvec core_similar_neighbors(cores.size(),intvec{});
    parlay::parallel_for(0, cores.size(), [&](const size_t j) {
        size_t temp=cores[j];
        for(unsigned k=0;k<geshu-i;k++)
        {
            for( auto elem1:a.compress[temp][k]){
                core_similar_neighbors[j].push_back(elem1);
            }         
        }
        
    });

   


    
    ClusterCores(cores,core_similar_neighbors,&clustering);
    AttachNoncoresToClusters(cores,core_similar_neighbors,&clustering);
    

    const auto end=std::chrono::steady_clock().now();
    tms total_time=end-begin;
     a.exec_time = total_time.count();
    //a.output["execution time"]= std::to_string(a.exec_time);
     //printf("Total time: %lf \n",total_time.count());
    // a.write_others(Cluster_ID-1);
     a.writecluster_p("../output/",clustering,sim,u);
    //std::cout<<Cluster_ID-1<<'\n';

}
