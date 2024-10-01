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
#include "utils.h"

//--------------------------------------------------------------------- Utility functions ------------------------------------------------------------------------------

void Algorithm::writecluster(std::string folder,intintvec Cluster,double sim,unsigned int u){
    //std::cout << "core: \n";把core保存到csv文件中，core是运行过程中的变量，只存在于运行过程，我们要看到他的话要么print，要么写道csv文件里再看
    std::string file = folder + "cluster_"+output["algo"]+"_"+hg.dataset+".csv";// output是什么？ 好像就是打了个表 用关键词去对应另一个
  //  std::cout<<"writing to: "<<file<<"\n";
    std::stringstream ss;
    for(size_t i = 0; i < Cluster_num; i++){
        for(auto elem: Cluster[i]){
            ss<< std::to_string(elem)<<',';
        }
        ss<<"\n";
    }
    // for(auto elem: Cluster)
    // {
    //     // std::cout<<elem.first<<","<<elem.second<<"\n";
    //    // ss << std::to_string(elem.first) << "," << std::to_string(elem.second) << "\n";
    // }
    std::string sim1=std::to_string(sim);
    std::string u1=std::to_string(u+1);
    std::string lujin="../output/";//
    std::string file1=lujin+"cluster_"+output["algo"]+"_"+hg.dataset+"_"+u1+"_"+sim1+".csv";
    std::cout<<"writing to: "<<file1<<"\n";
    std::ofstream out(file1.c_str());  //file中的../output/字段使得该文件无法创建
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

// ------------------------------------------------------------------------ Peel ------------------------------------------------------------------------------------



 void Algorithm::write_to_Cluster_num(size_t b)
 {
        Cluster_num=b;
 }

bool geometricMean(double sim,const intintvec & Hyperedge,unsigned int num1,unsigned int num2,size_t e1,size_t e2) {        /**********得到相似度*/
    
    double result = sqrt(static_cast<double>(num1) * num2);
    double divisionResult = static_cast<double>(num1>num2?num2:num1) / result; //?有问题？

    if(divisionResult<sim) return false;

    if(((double)num1<sim*sim*num2)||((double)num2<sim*sim*num1)) return false;//?有问题？
    std::set<size_t> gongtongdingdian;
   // unsigned int count=0;
    divisionResult=0;
    for (unsigned int i=0;i<Hyperedge[e1].size();i++) {
             // if((static_cast<double>(count+Hyperedge[e1].size()-i) / result)<sim) break;
              auto it = std::find(Hyperedge[e2].begin(), Hyperedge[e2].end(), Hyperedge[e1][i]);
              if (it != Hyperedge[e2].end()) {
                //count++;// 这里有提前终止的可能？
                gongtongdingdian.insert(Hyperedge[e1][i]);
                divisionResult = static_cast<double>(gongtongdingdian.size()) / result;
                if(divisionResult>=sim) break;
            }
    }
   // divisionResult = static_cast<double>(count) / result;
    if(divisionResult>=sim) return true;
    else return false;
}


//纯裸的base_line,这里n_sim不该存在,在过程中check_core,以及计算n_sim
bool checkcore_base_line(double sim,unsigned int u,size_t i,const intintvec & Hyperedge)
{
    
    size_t count=0;
    for(size_t j=0;j<Hyperedge.size();j++)
    {
        if(i!=j)
        {
            if(geometricMean(sim,Hyperedge,Hyperedge[i].size(),Hyperedge[j].size(),i,j)) count+=1;
        }
        if(count>=u) return true;
    }
    if(count>=u) return true;
    else return false;
}
void compute_base_line_n_sim(double sim,size_t i,intvec & n_sim,const intintvec & Hyperedge)
{

    for(size_t j=0;j<Hyperedge.size();j++)
    {
        if(i!=j)
        {
            if(geometricMean(sim,Hyperedge,Hyperedge[i].size(),Hyperedge[j].size(),i,j)) n_sim.push_back(j);
        }
    }

}
void base_line(std::string dataset, intintvec e_id_to_edge, Algorithm& a,double sim,unsigned int u)
{

    a.output["algo"] = "base_line";
    clock_t start, end;
    a.output["init_time"] = std::to_string(double(0) / double(CLOCKS_PER_SEC));
    start = clock();
    /**********处理过程*********/
    intvec n_sim;
    intvec edge_to_cluster(e_id_to_edge.size(),0);//判断超边是否已经被加入过当前cluster
    intintvec Cluster(e_id_to_edge.size(),intvec{});//保存cluster 
    //  还是要把cluster里面超边的连接信息保存下来，不能只保存一个超边编号，会导致只知道他们在一个cluster里面而不知道他们的连接关系(实际上是可以知道的 不过要凭超边编号去原图中找连接关系)
    // 和普通图还是区别很大的 普通图中没有一条边和有一条边是不同的子图   超图中超边被加入cluster之后，由于共享节点的关系 不止会和引入它的超边有连接关系
    size_t Cluster_ID=1;
    for(size_t eid = 0; eid < e_id_to_edge.size(); eid++){
         if(checkcore_base_line(sim,u,eid,e_id_to_edge)&&edge_to_cluster[eid]==0){
            edge_to_cluster[eid]=Cluster_ID;
            Cluster[Cluster_ID-1].push_back(eid); //??????
            std::queue<size_t> Q;
            Q.push(eid);
            while(Q.empty()!=true){
                size_t temp=Q.front();
                Q.pop();
                compute_base_line_n_sim(sim,temp,n_sim,e_id_to_edge);
                for (const auto& elem1:n_sim) {
                    if(edge_to_cluster[elem1]!=Cluster_ID){//还没被加入过当前cluster
                        
                        edge_to_cluster[elem1]=Cluster_ID;
                        Cluster[Cluster_ID-1].push_back(elem1);  //??????    直接在这里定义一个局部变量  最后调用函数将结果保存到csv文件中

                        if(checkcore_base_line(sim,u,elem1,e_id_to_edge)){
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
    a.writecluster("../output/",Cluster,sim,u);
    std::cout<<Cluster_ID-1<<'\n';
}

//base_line+反向超边索引，这里也不该用n_sim,在过程中check_core,以及计算n_sim
bool geometricMean2(double sim,intIntMap& node_index,const intintvec & Hyperedge,unsigned int num1,unsigned int num2,size_t e1,size_t e2) {        /**********得到相似度*/
    
    double result = sqrt(static_cast<double>(num1) * num2);
    double divisionResult = static_cast<double>(num1>num2?num2:num1) / result;  // 这里如果有重复节点 将会被统计进下界里面  但是只会多不会少啊

    if(divisionResult<sim) return false;  //opt1  这里导致的和gs index差异？
    //if(((double)num1<sim*sim*num2)||((double)num2<sim*sim*num1)) return false; //opt4   这里导致的和gs index差异？
    unsigned int count=0;
    divisionResult=0;
    std::set<size_t> gongtongdingdian;
    for (unsigned int j=0;j<Hyperedge[e1].size();j++) {
             // if((static_cast<double>(count+Hyperedge[e1].size()-j) / result)<sim) break;  //opt2
              auto it = std::find(Hyperedge[e2].begin(), Hyperedge[e2].end(), Hyperedge[e1][j]);
              if (it != Hyperedge[e2].end()) {
               // count++;// 这里有提前终止的可能？如果超边里面有重复顶点 这里直接用count统计的话有被统计多次的可能 同时也有只统计一次的时候 会造成差异？
                gongtongdingdian.insert(Hyperedge[e1][j]);
                divisionResult = static_cast<double>(gongtongdingdian.size()) / result;
                if(divisionResult>=sim) break;   //opt3
              }    
    }
   // divisionResult = static_cast<double>(count) / result;
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
void optimize_index(std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,Algorithm& a,double sim,unsigned int u)
{
    
   a.output["algo"] = "optimize_index";
    clock_t start, end;
    start = clock();
    /***************得到Nv邻域*******************/  
    //intintvec N_Sim(e_id_to_edge.size(),intvec{});// 保存各自亲密度高的边
    size_t N = init_nodes.size();//顶点数
    intintvec tables(N,intvec{});//索引表 有多少个顶点就建多少个表
    for(size_t i=0;i<e_id_to_edge.size();i++)//建表
    {
        for (const auto& elem:e_id_to_edge[i]) {
             auto j = node_index[elem];
             tables[j].push_back(i);
        }
    }

    clock_t e_tm = clock();
    a.output["init_time"] = std::to_string(double(e_tm - start) / double(CLOCKS_PER_SEC));
    //start = clock();
    /**********处理过程*********/
    intvec edge_to_cluster(e_id_to_edge.size(),0);//判断超边是否已经被加入过当前cluster
    intvec n_sim;
    intintvec Cluster(e_id_to_edge.size(),intvec{});//保存cluster 
    //  还是要把cluster里面超边的连接信息保存下来，不能只保存一个超边编号，会导致只知道他们在一个cluster里面而不知道他们的连接关系(实际上是可以知道的 不过要凭超边编号去原图中找连接关系)
    // 和普通图还是区别很大的 普通图中没有一条边和有一条边是不同的子图   超图中超边被加入cluster之后，由于共享节点的关系 不止会和引入它的超边有连接关系
    size_t Cluster_ID=1;
    //std::cout<<N_Sim[462][0]<<' '<<N_Sim[462][1]<<' '<<N_Sim[462][2]<<' '<<N_Sim[462][3]<<' '<<N_Sim[462][4]<<' '<<N_Sim[462][5]<<'\n';
    for(size_t eid = 0; eid < e_id_to_edge.size(); eid++){
         if(checkcore_optimize_index(sim,u,eid,tables,node_index,e_id_to_edge)&&edge_to_cluster[eid]==0){
            edge_to_cluster[eid]=Cluster_ID;
            Cluster[Cluster_ID-1].push_back(eid); //??????
            std::queue<size_t> Q;
            Q.push(eid);
            while(Q.empty()!=true){
                size_t temp=Q.front();
                Q.pop();
                compute_optimize_index_n_sim(sim,temp,tables,n_sim,node_index,e_id_to_edge);
                for (const auto& elem1:n_sim) {
                    if(edge_to_cluster[elem1]!=Cluster_ID){//还没被加入过当前cluster
                        
                        edge_to_cluster[elem1]=Cluster_ID;
                        Cluster[Cluster_ID-1].push_back(elem1);  //??????    直接在这里定义一个局部变量  最后调用函数将结果保存到csv文件中

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
    a.writecluster("../output/",Cluster,sim,u);
    std::cout<<Cluster_ID-1<<'\n';

}

//没啥用版本
bool geometricMean3(double sim,const intintvec & Hyperedge,unsigned int num1,unsigned int num2,size_t e1,size_t e2) {        /**********得到相似度*/
    
    double result = sqrt(static_cast<double>(num1) * num2);
    double divisionResult = static_cast<double>(num1>num2?num2:num1) / result;  

    if(divisionResult<sim) return false;  //opt1

    unsigned int count=0;
    divisionResult=0;
    for (unsigned int j=0;j<Hyperedge[e1].size();j++) {
              if((static_cast<double>(count+Hyperedge[e1].size()-j) / result)<sim) break;  //opt2
              auto it = std::find(Hyperedge[e2].begin(), Hyperedge[e2].end(), Hyperedge[e1][j]);
              if (it != Hyperedge[e2].end()) {
                count++;// 这里有提前终止的可能？
                divisionResult = static_cast<double>(count) / result;
                if(divisionResult>=sim) break;   //opt3
              }    
    }
   // divisionResult = static_cast<double>(count) / result;
    if(divisionResult>=sim) return true;
    else return false;
}
void optimize_indexII(std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,Algorithm& a,double sim,unsigned int u)
{
    
   a.output["algo"] = "optimize_indexII";
    clock_t start, end;
    start = clock();
    /***************得到Nv邻域*******************/  
    //复杂度太高？  快速得到超边的邻接超边方法？ 构建索引将包含某个顶点的超边用一个表组织在一起 HGMatch中有这个构建方法？
    intintvec N_Sim(e_id_to_edge.size(),intvec{});// 保存各自亲密度高于sim的边
    intintvec Neighbor(e_id_to_edge.size(),intvec{});
    size_t N = init_nodes.size();//顶点数
    intintvec tables(N,intvec{});//索引表 有多少个顶点就建多少个表
    for(size_t i=0;i<e_id_to_edge.size();i++)//建表
    {
        for (const auto& elem:e_id_to_edge[i]) {
             auto j = node_index[elem];
             tables[j].push_back(i);
        }
    }

    for(size_t i=0;i<N;i++)
    {
        for(size_t j=0;j<tables[i].size();j++)
        {

            for(size_t k=j+1;k<tables[i].size();k++)
            {
                //已经算过了就直接continue 这样的话需要将所有亲密度大于0的邻接超边都保存 而不是只保存大于sim的到n_sim中去
                //这样的话时间 空间开销 都会增大  耗时变大十倍
                auto item2 = std::find(Neighbor[tables[i][j]].begin(),Neighbor[tables[i][j]].end(),tables[i][k]);
                if(item2!=Neighbor[tables[i][j]].end()) continue;

                if(geometricMean3(sim,e_id_to_edge,e_id_to_edge[tables[i][j]].size(),e_id_to_edge[tables[i][k]].size(),tables[i][j],tables[i][k])) {
                 
                 if(tables[i][j]!=tables[i][k]){
                     N_Sim[tables[i][j]].push_back(tables[i][k]);
                     N_Sim[tables[i][k]].push_back(tables[i][j]);
                  }
                  
                }  
                 Neighbor[tables[i][j]].push_back(tables[i][k]);
                 Neighbor[tables[i][k]].push_back(tables[i][j]);
            }
          
        }
    }

    clock_t e_tm = clock();
    a.output["init_time"] = std::to_string(double(e_tm - start) / double(CLOCKS_PER_SEC));
    start = clock();
    /**********处理过程*********/
    intvec edge_to_cluster(e_id_to_edge.size(),0);//判断超边是否已经被加入过当前cluster
    intintvec Cluster(e_id_to_edge.size(),intvec{});//保存cluster 
    //  还是要把cluster里面超边的连接信息保存下来，不能只保存一个超边编号，会导致只知道他们在一个cluster里面而不知道他们的连接关系(实际上是可以知道的 不过要凭超边编号去原图中找连接关系)
    // 和普通图还是区别很大的 普通图中没有一条边和有一条边是不同的子图   超图中超边被加入cluster之后，由于共享节点的关系 不止会和引入它的超边有连接关系
    size_t Cluster_ID=1;
    //std::cout<<N_Sim[462][0]<<' '<<N_Sim[462][1]<<' '<<N_Sim[462][2]<<' '<<N_Sim[462][3]<<' '<<N_Sim[462][4]<<' '<<N_Sim[462][5]<<'\n';
    for(size_t eid = 0; eid < e_id_to_edge.size(); eid++){
         if(N_Sim[eid].size()>=u&&edge_to_cluster[eid]==0){
            edge_to_cluster[eid]=Cluster_ID;
            Cluster[Cluster_ID-1].push_back(eid); //??????
            std::queue<size_t> Q;
            Q.push(eid);
            while(Q.empty()!=true){
                size_t temp=Q.front();
                Q.pop();
                for (const auto& elem1:N_Sim[temp]) {
                    if(edge_to_cluster[elem1]!=Cluster_ID){//还没被加入过当前cluster
                        
                        edge_to_cluster[elem1]=Cluster_ID;
                        Cluster[Cluster_ID-1].push_back(elem1);  //??????    直接在这里定义一个局部变量  最后调用函数将结果保存到csv文件中

                        if(N_Sim[elem1].size()>=u){
                            Q.push(elem1);
                        }
                    }

                       
                }
            }
            Cluster_ID++;
         }
    }
    end = clock();
    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.write_to_Cluster_num(Cluster_ID-1);
    a.writecluster("../output/",Cluster,sim,u);
    std::cout<<Cluster_ID-1<<'\n';

}


//compress index,optimize_compress中前面计算用geometricMean_compress,后面聚类用geometricMean
double geometricMean_compress(const intintvec & Hyperedge,unsigned int num1,unsigned int num2,size_t e1,size_t e2) {        /**********得到相似度*/
    
    double result = sqrt(static_cast<double>(num1) * num2);
    double divisionResult = 0;  
    unsigned int count=0;
    for (unsigned int j=0;j<Hyperedge[e1].size();j++) {
              auto it = std::find(Hyperedge[e2].begin(), Hyperedge[e2].end(), Hyperedge[e1][j]);
              if (it != Hyperedge[e2].end()) {
              //  if(node_index[Hyperedge[e1][j]]<i) {return -1;}//已经在前面的表中算过了
                count++;
              }    
    }
    divisionResult = static_cast<double>(count) / result;
    return divisionResult;
}
bool checkcore(double sim,unsigned int u,size_t i,intintdouprvec & tables,const compresstuplevec &compress,intIntMap& node_index,const intintvec & Hyperedge)
{
    
    std::set<size_t> n_sim;
    for(auto tup:compress)
    {
        if(n_sim.size()>=u) return true;
        if(std::get<4>(tup)<sim) break;
        for(auto v:Hyperedge[i]){
                auto j = node_index[v];
                for(size_t k=0;k<tables[j].size();k++)
                {
                    if(tables[j][k].first<sim) continue;
                    if(tables[j][k].second>=std::get<0>(tup)&&tables[j][k].second<=std::get<1>(tup)){
                        if(geometricMean(sim,Hyperedge,Hyperedge[i].size(),Hyperedge[tables[j][k].second].size(),i,tables[j][k].second))
                         {
                            if(i!=tables[j][k].second) n_sim.insert(tables[j][k].second);//
                         }
                          if(n_sim.size()>=u) return true;
                    }
                }
        }
    }
    if(n_sim.size()>=u) return true;
    else return false;
}
void compute_N_SIM(double sim,size_t i,intintdouprvec & tables,const compresstuplevec &compress,std::set<size_t> & N_SIM,intIntMap& node_index,const intintvec & Hyperedge)
{

    for(auto tup:compress)
    {
        if(std::get<4>(tup)<sim) break;
        for(auto v:Hyperedge[i]){
                auto j = node_index[v];
                for(size_t k=0;k<tables[j].size();k++)
                {
                    if(tables[j][k].first<sim) continue;
                    if(tables[j][k].second>=std::get<0>(tup)&&tables[j][k].second<=std::get<1>(tup)){
                        if(geometricMean(sim,Hyperedge,Hyperedge[i].size(),Hyperedge[tables[j][k].second].size(),i,tables[j][k].second))
                         {
                            if(i!=tables[j][k].second)  N_SIM.insert(tables[j][k].second);//
                         }
                    }
                }
        }
    }
}
void optimize_compress(std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,Algorithm& a,double sim,unsigned int u)
{
    
   a.output["algo"] = "optimize_compress";
    clock_t start, end;
    start = clock();
    /***************得到Nv邻域*******************/  
   
    douintset Neighbor;   //单个neighbor
    compresstuplevecvec compress(e_id_to_edge.size(),compresstuplevec{});  //存储压缩之后的索引结构
    size_t N = init_nodes.size();//顶点数
    intintdouprvec tables(N,intdouprvec{});//索引表 有多少个顶点就建多少个表
    std::cout<<"start!"<<'\n';

    for(size_t i=0;i<e_id_to_edge.size();i++)//建表
    {
        for (const auto& elem:e_id_to_edge[i]) {
             auto j = node_index[elem];
             tables[j].push_back(std::make_pair(0,i));
        }
    }
    //std::cout<<"2!"<<'\n';
    for(size_t i=0;i<e_id_to_edge.size();i++)
    {

        for(auto v:e_id_to_edge[i])
        {
            auto j = node_index[v];
            for(size_t k=0;k<tables[j].size();k++)
            {
                auto it = std::find_if(Neighbor.begin(), Neighbor.end(), [&](const std::pair<double, size_t>& p) {
                    return p.second == tables[j][k].second;
                });
                if(it!=Neighbor.end()) 
                {
                   if(-it->first>tables[j][k].first) tables[j][k].first=-it->first; //把负的转回来
                    continue;//已经算过了不用再算了
                }
                double xiangsidu=geometricMean_compress(e_id_to_edge,e_id_to_edge[i].size(),e_id_to_edge[tables[j][k].second].size(),i,tables[j][k].second);
                if(i!=tables[j][k].second){
                    if(xiangsidu>tables[j][k].first)  tables[j][k].first=xiangsidu;
                    Neighbor.insert(std::make_pair(-xiangsidu,tables[j][k].second));  //1-
                }
            }

        }
        
        //压缩一下存储到compress  先5个一压缩吧
        //for(auto tup:compress[i])  {  tup[0]   tup[1]}  是这样访问？
        //compresstuple  tuple_5 = std::make_tuple(1, 2, 3, 4, 5);
        unsigned int temp=0;
        double max_xiangsidu=0;
        double min_xiangsidu=0;
        size_t min_e=0;
        size_t max_e=0;
        for(auto p:Neighbor)
        {
            temp+=1;
            if(temp==1) {
                max_e=p.second;
                min_e=p.second;
                max_xiangsidu=-p.first;
                min_xiangsidu=-p.first;
            }
            else{
                if(p.second>max_e) max_e=p.second;
                if(p.second<min_e) min_e=p.second;
                min_xiangsidu=-p.first;
            }
            if(temp==7) {
                compress[i].push_back(std::make_tuple(min_e,max_e,temp,min_xiangsidu,max_xiangsidu));
                temp=0;
            }
        }
        if(temp!=0)
        {
            compress[i].push_back(std::make_tuple(min_e,max_e,temp,min_xiangsidu,max_xiangsidu));//把负的相似度调回来了
        }
       //  std::cout<<"start!"<<'\n';
        //清空Neighbor
        Neighbor.clear();

    }
    std::cout<<"compute ok!"<<'\n';
    

    clock_t e_tm = clock();
    a.output["init_time"] = std::to_string(double(e_tm - start) / double(CLOCKS_PER_SEC));

// 统计一下内存大小
    size_t  neicun=0;
    for(size_t i=0;i<tables.size();i++) {neicun+=tables[i].size()*16;}
    for(size_t i=0;i<compress.size();i++) {neicun+=compress[i].size()*40;}
    std::cout<<"neicun"<<' '<<neicun<<'\n';


    start = clock();
    /**********处理过程*********/
    intvec edge_to_cluster(e_id_to_edge.size(),0);//判断超边是否已经被加入过当前cluster
    intintvec Cluster(e_id_to_edge.size(),intvec{});//保存cluster 
    //  还是要把cluster里面超边的连接信息保存下来，不能只保存一个超边编号，会导致只知道他们在一个cluster里面而不知道他们的连接关系(实际上是可以知道的 不过要凭超边编号去原图中找连接关系)
    // 和普通图还是区别很大的 普通图中没有一条边和有一条边是不同的子图   超图中超边被加入cluster之后，由于共享节点的关系 不止会和引入它的超边有连接关系
    size_t Cluster_ID=1;

    std::set<size_t> N_SIM;//聚类时临时存储 与i相似超边编号，在每次结束时记得先把里面的内容清空

    for(size_t i=0;i<e_id_to_edge.size();i++){
        if(edge_to_cluster[i]!=0) continue;
        
        size_t count=0;
        for(auto x:compress[i]) {
            if(std::get<4>(x)>=sim) count+=std::get<2>(x);
            else break;
        }
        if(count<u) continue;

        if(!checkcore(sim,u,i,tables,compress[i],node_index,e_id_to_edge)) continue;
        
        edge_to_cluster[i]=Cluster_ID;
        Cluster[Cluster_ID-1].push_back(i); //??????
        std::queue<size_t> Q;
        Q.push(i);
        while(!Q.empty()){
            size_t temp=Q.front();
            Q.pop();
            //计算N_SIM_V
            compute_N_SIM(sim,temp,tables,compress[temp],N_SIM,node_index,e_id_to_edge);
            for ( auto elem1:N_SIM) {
                //double xiangsidu=geometricMean(sim,e_id_to_edge,,,,);
                if(edge_to_cluster[elem1]!=Cluster_ID){//还没被加入过当前cluster
                    
                    edge_to_cluster[elem1]=Cluster_ID;
                    Cluster[Cluster_ID-1].push_back(elem1);  //??????    直接在这里定义一个局部变量  最后调用函数将结果保存到csv文件中
                   
                    if(checkcore(sim,u,elem1,tables,compress[elem1],node_index,e_id_to_edge)){ //?
                        Q.push(elem1);
                    }
                }
       
            }
            N_SIM.clear();
        }
        
       // N_SIM.clear();
        Cluster_ID++;
       
    }
    end = clock();

    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.write_to_Cluster_num(Cluster_ID-1);
    a.writecluster("../output/",Cluster,sim,u);
    std::cout<<Cluster_ID-1<<'\n';

}

