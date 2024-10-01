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
#include "utils.h"
#include "HashSet.h"
#define epision 1e-10
//--------------------------------------------------------------------- Utility functions ------------------------------------------------------------------------------

void compress_index::writecluster(std::string folder,intintvec Cluster,double sim,unsigned int u){
    //std::cout << "core: \n";把core保存到csv文件中，core是运行过程中的变量，只存在于运行过程，我们要看到他的话要么print，要么写道csv文件里再看
    std::string file = folder + "cluster_"+output["algo"]+"_"+hg.dataset+".csv";// output是什么？ 好像就是打了个表 用关键词去对应另一个
   // std::cout<<"writing to: "<<file<<"\n";
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
    std::string lujin="../output/";
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


compress_index::compress_index(Hypergraph &H){
    hg = H;
    output["dataset"] = H.dataset;
}
compress_index::~compress_index(){}

// ------------------------------------------------------------------------ Peel ------------------------------------------------------------------------------------

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
    size_t N = init_nodes.size();//顶点数
    intintvec tables(N,intvec{});//索引表 有多少个顶点就建多少个表
    for(size_t i=0;i<e_id_to_edge.size();i++)//建表
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
                    //是邻居 加入neighbor
                    neighbors.insert(j);
                    break;
                }
        }
    }
    
}

bool geometricMean2(double sim,const intintvec & Hyperedge,unsigned int num1,unsigned int num2,size_t e1,size_t e2) {        /**********得到相似度*/
    
    double result = sqrt(static_cast<double>(num1) * num2);
    if(result==0) return false;//remove
    double divisionResult = static_cast<double>(num1>num2?num2:num1) / result;//?有问题？

    if(divisionResult<sim) return false;

    //if(((double)num1<sim*sim*num2)||((double)num2<sim*sim*num1)) return false;//?有问题？ 会导致差异？
    std::set<size_t> gongtongdingdian;
    //unsigned int count=0;
    divisionResult=0;
    for (unsigned int i=0;i<Hyperedge[e1].size();i++) {
             // if((static_cast<double>(count+Hyperedge[e1].size()-i) / result)<sim) break;
              auto it = std::find(Hyperedge[e2].begin(), Hyperedge[e2].end(), Hyperedge[e1][i]);
              if (it != Hyperedge[e2].end()) {
               // count++;// 这里有提前终止的可能？
               gongtongdingdian.insert(Hyperedge[e1][i]);
                divisionResult = static_cast<double>(gongtongdingdian.size()) / result;
                if(divisionResult>=sim) break;
            }
    }
   // divisionResult = static_cast<double>(count) / result;
    if(divisionResult>=sim) return true;
    else return false;
}
double geometricMean_compress_II(const intintvec & Hyperedge,unsigned int num1,unsigned int num2,size_t e1,size_t e2) {        /**********得到相似度*/
    
    double result = sqrt(static_cast<double>(num1) * num2);
    double divisionResult = 0;
    if(result==0) return 0;//remove  
    //unsigned int count=0;
    std::set<size_t> gongtongdingdian;
    for (unsigned int j=0;j<Hyperedge[e1].size();j++) {
              auto it = std::find(Hyperedge[e2].begin(), Hyperedge[e2].end(), Hyperedge[e1][j]);
              if (it != Hyperedge[e2].end()) {
              //  if(node_index[Hyperedge[e1][j]]<i) {return -1;}//已经在前面的表中算过了
               // count++;
               gongtongdingdian.insert(Hyperedge[e1][j]);
              }    
    }
    divisionResult = static_cast<double>(gongtongdingdian.size()) / result;
    return divisionResult;
}


//压缩算法，按照相似度分组，并且存储所有编号（按照与该边相似度大小顺序）或者掩码
bool checkcore_compress_II(double d,double sim,unsigned int u,size_t i,const intintvec &compress,const intintvec & Hyperedge)
{
    unsigned int count=0;
    // for(auto x:compress)
    // {
    //     if(count>=u) return true;
    //     if(!geometricMean2(sim,Hyperedge,Hyperedge[i].size(),Hyperedge[x[0]].size(),i,x[0])) break;
    //     if(geometricMean2(sim,Hyperedge,Hyperedge[i].size(),Hyperedge[x[x.size()-1]].size(),i,x[x.size()-1])){
    //         count=count+x.size();
    //         continue;
    //     }
    //     for(auto z:x)
    //     {
    //         if(geometricMean2(sim,Hyperedge,Hyperedge[i].size(),Hyperedge[z].size(),i,z)) count+=1;
    //         else break;
    //         if(count>=u) return true;
    //     }
    // }
    double d1[11]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
    for(unsigned k=0;k<(unsigned int)(1.0/d);k++)
    {
        if(count>=u) return true;
        //double temp=(9-k)*d;
        double temp=d1[9-k];
       //  if(temp+d-sim<epision) break;   //???
       // if(d1[10-k]-sim<epision) break;//这里用小于0还是小于epision？  先去掉吧
        if(temp>=sim) {count+=compress[k].size();continue;} //有等于号的时候不改？
        for(auto z:compress[k])
        {
            if(geometricMean2(sim,Hyperedge,Hyperedge[i].size(),Hyperedge[z].size(),i,z)) {
                count+=1;
                if(count>=u) return true;
            }
            else break;
            //if(count>=u) return true;
        }
        break;

    }
    if(count>=u) return true;
    else return false;
    
}
void compute_compress_II_n_sim(double d,double sim,size_t i,const intintvec &compress,intvec &n_sim,const intintvec & Hyperedge)
{
//   for(auto x:compress)
//   {
//     //if(x.empty()) continue;
//     if(!geometricMean2(sim,Hyperedge,Hyperedge[i].size(),Hyperedge[x[0]].size(),i,x[0])) break;
//     if(geometricMean2(sim,Hyperedge,Hyperedge[i].size(),Hyperedge[x[x.size()-1]].size(),i,x[x.size()-1])){
//         for(auto y:x) {n_sim.push_back(y);}
//         continue;
//     }
//     for(auto z:x)
//     {
//         if(geometricMean2(sim,Hyperedge,Hyperedge[i].size(),Hyperedge[z].size(),i,z)) n_sim.push_back(z);
//         else break;
//     }
//   }
    double d1[11]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};// 变化tau需要打新的表
    for(unsigned k=0;k<(unsigned int)(1.0/d);k++)
    {
       // double temp=(9-k)*d;
        double temp=d1[9-k];//9要变成4 9 49 99 等
        //if(temp+d-sim<epision) break;   //???
        //if(d1[10-k]-sim<epision) break;   //???
        if(temp>=sim) { 
            for(auto y:compress[k]) {n_sim.push_back(y);}
            continue;
        }
        for(auto z:compress[k])
        {
            if(geometricMean2(sim,Hyperedge,Hyperedge[i].size(),Hyperedge[z].size(),i,z)) {
                n_sim.push_back(z);
            }
            else break;
        }
        break;//

    }
   
}
void compress_II_construct(double d,std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u)
{
    
    clock_t start;
    start = clock();
    double d1[10]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
    /***************得到Nv邻域*******************/  
   
    douintset Neighbor;   //单个neighbor（包含相似度以及编号）
    std::set<size_t> computed_neighbor;//已经计算过了的邻居（不包含相似度）
    //compresstuplevecvec compress(e_id_to_edge.size(),compresstuplevec{});  //存储压缩之后的索引结构
    //size_t N = init_nodes.size();//顶点数
    //intintvec tables(N,intvec{});//索引表 有多少个顶点就建多少个表
    //std::vector<intintvec> compress(e_id_to_edge.size(),intintvec{});
    //double  d=0.1; //划分间隙  固定区间划分还是根据第一个的相似度  固定长度划分
    unsigned int geshu1=(unsigned int)(1.0/d);
    std::vector<std::vector<std::vector<size_t>>> compress(e_id_to_edge.size(),std::vector<std::vector<size_t>>(geshu1));
    
    std::cout<<"start!"<<'\n';

    // for(size_t i=0;i<e_id_to_edge.size();i++)//建表
    // {
    //     for (const auto& elem:e_id_to_edge[i]) {
    //          auto j = node_index[elem];
    //          tables[j].push_back(i);
    //     }
    // }
    tables_construct(e_id_to_edge,init_nodes,node_index,a);
   // std::cout<<"2!"<<'\n';
    for(size_t i=0;i<e_id_to_edge.size();i++)
    {
        generate_neighbors_by_tables2(i,e_id_to_edge,node_index,a,computed_neighbor);
        //generate_neighbors_by_enumerate2(i,e_id_to_edge,computed_neighbor);
        for(auto v:computed_neighbor)
        {
            double xiangsidu=geometricMean_compress_II(e_id_to_edge,e_id_to_edge[i].size(),e_id_to_edge[v].size(),i,v);
            if(i!=v){
                Neighbor.insert(std::make_pair(-xiangsidu,v));  //1-
            }
        }
        // for(auto v:e_id_to_edge[i])
        // {
        //     auto j = node_index[v];
        //     for(size_t k=0;k<tables[j].size();k++)
        //     {
        //         auto it = computed_neighbor.find(tables[j][k]);
        //         if(it!=computed_neighbor.end()) 
        //         {
        //             continue;//已经算过了不用再算了
        //         }
        //         else {computed_neighbor.insert(tables[j][k]);}

        //         double xiangsidu=geometricMean_compress_II(e_id_to_edge,e_id_to_edge[i].size(),e_id_to_edge[tables[j][k]].size(),i,tables[j][k]);
        //         if(i!=tables[j][k]){
        //             // if(xiangsidu>tables[j][k].first)  tables[j][k].first=xiangsidu;
        //             Neighbor.insert(std::make_pair(-xiangsidu,tables[j][k]));  //1-
        //         }
        //     }

        // }
        
        double t=0;
        unsigned int count_1=0;

        for(auto p:Neighbor)
        {
            count_1=0;
            t=-p.first;
            while(t>=d1[count_1])//t>=d t-d>=epision  改成这样会不会边界情况好了?  没有好    有些从对变错  有些从错变对
            {
                count_1+=1;
                //t=t-d;
                if(count_1>=9) break;
            }
            if(count_1>9) count_1=9;
            compress[i][geshu1-count_1-1].push_back(p.second);
        }
        computed_neighbor.clear();
        Neighbor.clear();

    }
    std::cout<<"compute ok!"<<'\n';
    

    clock_t e_tm = clock();
    a.init_time=double(e_tm - start) / double(CLOCKS_PER_SEC);
    a.output["inition time"]= std::to_string(a.init_time);
    a.write_indexII(compress);
// 统计一下内存大小
    size_t  neicun=0;
    for(size_t i=0;i<a.tables.size();i++) {neicun+=a.tables[i].size()*8;}
    for(size_t i=0;i<compress.size();i++) {
        for(auto x:compress[i])
        {
            neicun+=x.size()*8;
        }
    }
    a.neicun=neicun;
    std::cout<<"neicun"<<' '<<neicun<<'\n';

}
void compress_II_cluster(double d,std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u)
{
    clock_t start, end;
    a.output["algo"] = "optimize_compress_II";
    start = clock();
    /**********处理过程*********/
    intvec edge_to_cluster(e_id_to_edge.size(),0);//判断超边是否已经被加入过当前cluster
    intintvec Cluster(e_id_to_edge.size(),intvec{});//保存cluster 
    //  还是要把cluster里面超边的连接信息保存下来，不能只保存一个超边编号，会导致只知道他们在一个cluster里面而不知道他们的连接关系(实际上是可以知道的 不过要凭超边编号去原图中找连接关系)
    // 和普通图还是区别很大的 普通图中没有一条边和有一条边是不同的子图   超图中超边被加入cluster之后，由于共享节点的关系 不止会和引入它的超边有连接关系
    size_t Cluster_ID=1;

    intvec n_sim;//聚类时临时存储 与i相似超边编号，在每次结束时记得先把里面的内容清空

    for(size_t i=0;i<e_id_to_edge.size();i++){
        if(a.compress[i].empty()) continue;
        if(edge_to_cluster[i]!=0) continue;
        if(!checkcore_compress_II(d,sim,u,i,a.compress[i],e_id_to_edge)) continue;
        
        edge_to_cluster[i]=Cluster_ID;
        Cluster[Cluster_ID-1].push_back(i); //??????
        std::queue<size_t> Q;
        Q.push(i);
        while(!Q.empty()){
            size_t temp=Q.front();
            Q.pop();
            //计算N_SIM_V
            compute_compress_II_n_sim(d,sim,temp,a.compress[temp],n_sim,e_id_to_edge);
            for ( auto elem1:n_sim) {
                //double xiangsidu=geometricMean(sim,e_id_to_edge,,,,);
                if(edge_to_cluster[elem1]!=Cluster_ID){//还没被加入过当前cluster
                    
                    edge_to_cluster[elem1]=Cluster_ID;
                    Cluster[Cluster_ID-1].push_back(elem1);  //??????    直接在这里定义一个局部变量  最后调用函数将结果保存到csv文件中
                   
                    if(checkcore_compress_II(d,sim,u,elem1,a.compress[elem1],e_id_to_edge)){ //?
                        Q.push(elem1);
                    }
                }
       
            }
            n_sim.clear();
        }
        
       // N_SIM.clear();
        Cluster_ID++;
       
    }
    end = clock();

    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.write_others(Cluster_ID-1);
    a.writecluster("../output/",Cluster,sim,u);
    std::cout<<Cluster_ID-1<<'\n';
}

//改进压缩索引，增加cluster_index
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

// bool checkcore_compress_IIII(double d,unsigned int k, std::map<size_t,int>&Cluster_i_jia_1,double d_segment_count,double sim,unsigned int u,size_t i,const intintvec &compress,const intintvec & Hyperedge)
// {
//     //使用map时
//     // if(k<9){
//     //     auto it=Cluster_Index[k+1].find(i);
//     //     if(it!=Cluster_Index[k+1].end()&&it->second>=u) return true;
//     // }
    
//     //使用set时
    
//     auto it=Cluster_i_jia_1.find(i);
//     if(it!=Cluster_i_jia_1.end()&&it->second<u) return false;
//     double d1[11]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
//     unsigned int count=0;
//     for(unsigned k=0;k<(unsigned int)(1.0/d);k++)
//     {
//         if(count>=u) return true;
//         //double temp=(9-k)*d;
//         double temp=d1[9-k];
//         if(temp+d-sim<epision) break;
//         if(temp>=sim) {count+=compress[k].size();continue;}
//         for(auto z:compress[k])
//         {
//             if(geometricMean2(sim,Hyperedge,Hyperedge[i].size(),Hyperedge[z].size(),i,z)) {
//                 count+=1;
//                 if(count>=u) return true;
//             }
//             else break;
//             //if(count>=u) return true;
//         }

//     }
//     if(count>=u) return true;
//     else return false;
    
// }
bool checkcore_compress_IIIII(double d,double sim,unsigned int u,size_t i,const intintvec &compress,const intintvec & Hyperedge)
{
    unsigned int count=0;
    double d1[11]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};//变化tau需要打新的表
    for(unsigned k=0;k<(unsigned int)(1.0/d);k++)
    {
        if(count>=u) return true;
        //double temp=(9-k)*d;
        double temp=d1[9-k];//9要变成geshu-1 每次手动改变成4 9 49 99 199也行
       // if(temp+d-sim<epision) break;
        if(temp>=sim) {count+=compress[k].size();continue;}//
        for(auto z:compress[k])
        {
            if(geometricMean2(sim,Hyperedge,Hyperedge[i].size(),Hyperedge[z].size(),i,z)) {
                count+=1;
                if(count>=u) return true;
            }
            else break;
            //if(count>=u) return true;
        }
        break;

    }
    if(count>=u) return true;
    else return false;
    
}

double geometricMean_compress_IIII(std::vector<HashSet> &hash,const intintvec & Hyperedge,unsigned int num1,unsigned int num2,size_t e1,size_t e2) {        /**********得到相似度*/
    
    double result = sqrt(static_cast<double>(num1) * num2);
    double divisionResult = 0;  
   // unsigned int count=0;
    if(result==0) return 0;// remove
    size_t temp=0;
    std::set<size_t> gongtongdingdian;
    if(num1>num2) {temp=e1;e1=e2;e2=temp;}// e1代表两条超边中规模较小的那条
    
    for (unsigned int j=0;j<Hyperedge[e1].size();j++) {
              //auto j=node_index[x];
             // auto it = hash[e2].find(Hyperedge[e1][j]); //hash[e2][j] 会在没有找到的时候自动创建一个？
              if (hash[e2].find(Hyperedge[e1][j])==1) {
               // count++;
               gongtongdingdian.insert(Hyperedge[e1][j]);
              }    
    }
    divisionResult = static_cast<double>(gongtongdingdian.size()) / result;
    return divisionResult;
}

void compress_IIII_construct(double d,double d_segment_count,std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u)
{
    
    clock_t start;
    start = clock();
    /***************得到Nv邻域*******************/  
    //double d1[5]={0.2,0.4,0.6,0.8,1.0};
    double d1[10]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};//变化tau需要打新的表
    // double d1[50]={0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,
    // 0.48, 0.5,0.52,0.54,0.56,0.58,0.6,0.62,0.64,0.66,0.68,0.7,0.72,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,
    // 0.98,1.0};
//     double d1[100]={0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,
// 0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,
// 0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,
// 0.6,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.7,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,
// 0.8,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.0};
//     double d1[200]={0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.075,0.08,0.085,0.09,0.095,
// 0.1,0.105,0.11,0.115,0.12,0.125,0.13,0.135,0.14,0.145,0.15,0.155,0.16,0.165,0.17,0.175,0.18,0.185,0.19,0.195,
// 0.2,0.205,0.21,0.215,0.22,0.225,0.23,0.235,0.24,0.245,0.25,0.255,0.26,0.265,0.27,0.275,0.28,0.285,0.29,0.295,
// 0.3,0.305,0.31,0.315,0.32,0.325,0.33,0.335,0.34,0.345,0.35,0.355,0.36,0.365,0.37,0.375,0.38,0.385,0.39,0.395,
// 0.4,0.405,0.41,0.415,0.42,0.425,0.43,0.435,0.44,0.445,0.45,0.455,0.46,0.465,0.47,0.475,0.48,0.485,0.49,0.495,
// 0.5,0.505,0.51,0.515,0.52,0.525,0.53,0.535,0.54,0.545,0.55,0.555,0.56,0.565,0.57,0.575,0.58,0.585,0.59,0.595,
// 0.6,0.605,0.61,0.615,0.62,0.625,0.63,0.635,0.64,0.645,0.65,0.655,0.66,0.665,0.67,0.675,0.68,0.685,0.69,0.695,
// 0.7,0.705,0.71,0.715,0.72,0.725,0.73,0.735,0.74,0.745,0.75,0.755,0.76,0.765,0.77,0.775,0.78,0.785,0.79,0.795,
// 0.8,0.805,0.81,0.815,0.82,0.825,0.83,0.835,0.84,0.845,0.85,0.855,0.86,0.865,0.87,0.875,0.88,0.885,0.89,0.895,
// 0.9,0.905,0.91,0.915,0.92,0.925,0.93,0.935,0.94,0.945,0.95,0.955,0.96,0.965,0.97,0.975,0.98,0.985,0.99,0.995,1.0};

    std::vector<core_pair3> Neighbor;   //单个neighbor（包含相似度以及编号）
    std::set<size_t> computed_neighbor;//已经计算过了的邻居（不包含相似度）
    //compresstuplevecvec compress(e_id_to_edge.size(),compresstuplevec{});  
    //size_t N = init_nodes.size();//顶点数
    //intintvec tables(N,intvec{});//索引表 有多少个顶点就建多少个表   这里该怎么不用N限制大小？同时又能用下标访问？（为了insert考虑）
    //函数里的tables可以限定大小，类自带的tables不限定大小，函数里面的tables计算完成之后再复制给类里面的tables
    //core-order  tables  cluster-index compress  neighbor 这几个都可以这样处理  在函数里面先算，最后再写到类里面的对应变量去存储
    //std::vector<intintvec> compress(e_id_to_edge.size(),intintvec{});//存储压缩之后的索引结构
    unsigned int geshu1=(unsigned int)(1.0/d);
    std::vector<std::vector<size_t>> compress(geshu1,std::vector<size_t>{});
    //double  d=0.1; //划分间隙  固定区间划分还是根据第一个的相似度  固定长度划分?

    //double d_segment_count=0.1;  //分段统计间隙
    unsigned int geshu=(unsigned int)(1.0/d_segment_count);//会不会出问题？
    //std::vector<std::map<size_t,int>> Cluster_Index(geshu,std::map<size_t,int>{});//用map嘛这里
    std::vector<std::vector<core_pair2>> Cluster_Index(geshu,std::vector<core_pair2>{});

    std::vector<int> segment_count;//初始化为空
    //std::vector<std::vector<unsigned int>> segment_count(e_id_to_edge.size(),std::vector<unsigned int>{})
    std::vector< HashSet > hash;
    hash.reserve(e_id_to_edge.size());
    //std::vector<std::vector<short>> hash1(e_id_to_edge.size(),std::vector<short>(5000));
    for(size_t i=0;i<e_id_to_edge.size();i++)
    {
        hash.emplace_back(e_id_to_edge[i].size());
    }
    
    std::cout<<"start!"<<'\n';
    tables_construct(e_id_to_edge,init_nodes,node_index,a);
    for(size_t i=0;i<e_id_to_edge.size();i++)//建表
    {
        //HashSet h1(e_id_to_edge[i].size()+2);
        for (const auto& elem:e_id_to_edge[i]) {
            // auto j = node_index[elem];
             hash[i].insert(elem);  //超边中有相同顶点只被添加一次 按理说不该有相同值
          //  h1.insert(elem);
        }
       // hash.push_back(h1);
    }
   // std::cout<<"2!"<<'\n';
    core_pair2 cp;
    core_pair3 cp3;
    for(size_t i=0;i<e_id_to_edge.size();i++)
    {

        generate_neighbors_by_tables2(i,e_id_to_edge,node_index,a,computed_neighbor);
        //generate_neighbors_by_enumerate2(i,e_id_to_edge,computed_neighbor);
        for(auto v:computed_neighbor)
        {
            double xiangsidu=geometricMean_compress_IIII(hash,e_id_to_edge,e_id_to_edge[i].size(),e_id_to_edge[v].size(),i,v);//有问题 用了之后结果不一样
          //  double xiangsidu=geometricMean_compress_II(e_id_to_edge,e_id_to_edge[i].size(),e_id_to_edge[v].size(),i,v);
            if(i!=v){
                // Neighbor.insert(std::make_pair(-xiangsidu,v));  //1-
                cp3.eid=v;
                cp3.xiangsidu=-xiangsidu;
                Neighbor.push_back(cp3);  //1-
            }
        }
        sort(Neighbor.begin(),Neighbor.end(),comp3);
        double t1=0;
        unsigned int count_1=0;

        for(auto p:Neighbor)
        {
            count_1=0;
            t1=-p.xiangsidu;
            while(t1>=d1[count_1]) //t1-d>=epision?t>=d 改成这样会不会边界情况好了 但是其他情况错了？ 在边界时候还是没有好
            {
                count_1+=1;
               // t1=t1-d;
               if(count_1>=geshu-1) break;  //9也要变成和tau相对应的
            }                     //这里算错了导致存的位置不对 sim=  0.6   0.8的时候只有一个能对
            if(count_1>geshu-1) count_1=geshu-1;
            compress[geshu1-count_1-1].push_back(p.eid);
            // if(t1<0.1) {compress[9].push_back(p.second);continue;}
            // if(t1<0.2) {compress[8].push_back(p.second);continue;}
            // if(t1<0.3) {compress[7].push_back(p.second);continue;}
            // if(t1<0.4) {compress[6].push_back(p.second);continue;}
            // if(t1<0.5) {compress[5].push_back(p.second);continue;}
            // if(t1<0.6) {compress[4].push_back(p.second);continue;}
            // if(t1<0.7) {compress[3].push_back(p.second);continue;}
            // if(t1<0.8) {compress[2].push_back(p.second);continue;}
            // if(t1<0.9) {compress[1].push_back(p.second);continue;}
            // compress[0].push_back(p.second);
        }
       
        
        // for(unsigned int k=0;k<geshu;k++) segment_count.push_back(0); 
        // for(unsigned int k=0;k<geshu;k++)
        // {
        //    segment_count[k]=compress[k].size();
        // }
        // for(unsigned int k=1;k<geshu;k++) segment_count[k]+=segment_count[k-1];// 总和而不是区间内部个数
        size_t segment_size=compress[0].size();
        for(unsigned int k=0;k<geshu;k++) {
            //大的存在后面 小的存在前面
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
       // segment_count.clear();
        Neighbor.clear();

    }

    for(unsigned int k=0;k<geshu;k++) {
        sort(Cluster_Index[k].begin(),Cluster_Index[k].end(),comp2);
    }


    std::cout<<"compute ok!"<<'\n';
    
    a.write_indexIIII(Cluster_Index);//这个时间需要不需要统计？
    clock_t e_tm = clock();
    a.init_time=double(e_tm - start) / double(CLOCKS_PER_SEC);
    a.output["inition time"]= std::to_string(a.init_time);
    //a.write_indexII(compress);
    
    // 统计一下内存大小
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
    std::cout<<"neicun"<<' '<<neicun<<'\n';
    // std::cout<<segment_count[4714][0]<<' '<<segment_count[4714][1]<<' '<<segment_count[4714][2]<<' '<<segment_count[4714][3]<<' '<<segment_count[4714][4]<<' '<<'\n';
    // std::cout<<segment_count[4714][5]<<' '<<segment_count[4714][6]<<' '<<segment_count[4714][7]<<' '<<segment_count[4714][8]<<' '<<segment_count[4714][9]<<' '<<'\n';
    
}
void compress_IIII_cluster(double d,double d_segment_count,std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u)
{
    a.output["algo"] = "optimize_compress_IIII";
    clock_t start, end;

    start = clock();
    /**********处理过程*********/
    intvec edge_to_cluster(e_id_to_edge.size(),0);//判断超边是否已经被加入过当前cluster
    intintvec Cluster(e_id_to_edge.size(),intvec{});//保存cluster 
    //  还是要把cluster里面超边的连接信息保存下来，不能只保存一个超边编号，会导致只知道他们在一个cluster里面而不知道他们的连接关系(实际上是可以知道的 不过要凭超边编号去原图中找连接关系)
    // 和普通图还是区别很大的 普通图中没有一条边和有一条边是不同的子图   超图中超边被加入cluster之后，由于共享节点的关系 不止会和引入它的超边有连接关系
    size_t Cluster_ID=1;
    intvec n_sim;//聚类时临时存储 与i相似超边编号，在每次结束时记得先把里面的内容清空

    //判断一下sim值是不是能够被d_segment_count整除
    double temp_sim=sim;
    unsigned int i=0;
    //std::map<size_t,int> cluster_i_jia_1={};
    double d1[11]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
    while (temp_sim>d1[i])// 这里要不要也改成打表的形式  避免多次与d相减？  只考虑非边界情况就可以了  这里好像阴差阳错都可以
    {
        i+=1;
        if(temp_sim<d1[i]) {temp_sim=-1;break;}
        if(i>=9) break;
    }

    // for(auto x:a.Cluster_Index[i-1])
    // {
    //     cluster_i_jia_1[x.second]=-x.first;
    // }


    for(auto x:a.Cluster_Index[i]){ //改Cluster index的数据结构和gs index一样
        //if(compress[x.second].empty()) continue;
        if(-x.count<u) break;
        if(edge_to_cluster[x.eid]!=0) continue;
        if(temp_sim>0){//epision而不是0，为了防止边界情况temo_sim是一个大于0的极小值的情况
            if(!checkcore_compress_IIIII(d,sim,u,x.eid,a.compress[x.eid],e_id_to_edge)) continue;   
        }//对于sim大于0.9的时候专门设置的
        edge_to_cluster[x.eid]=Cluster_ID;
        Cluster[Cluster_ID-1].push_back(x.eid); //??????
        std::queue<size_t> Q;
        Q.push(x.eid);
        while(!Q.empty()){
            size_t temp=Q.front();
            Q.pop();
            //计算N_SIM_V
            compute_compress_II_n_sim(d,sim,temp,a.compress[temp],n_sim,e_id_to_edge);//瓶颈在于这里core节点的邻居扩散？
            for ( auto elem1:n_sim) {
                //double xiangsidu=geometricMean(sim,e_id_to_edge,,,,);
                if(edge_to_cluster[elem1]!=Cluster_ID){//还没被加入过当前cluster
                    
                    edge_to_cluster[elem1]=Cluster_ID;
                    Cluster[Cluster_ID-1].push_back(elem1);  //??????    直接在这里定义一个局部变量  最后调用函数将结果保存到csv文件中
                    // if(checkcore_compress_IIIII(d,i,cluster_i_jia_1,d_segment_count,sim,u,elem1,a.compress[elem1],e_id_to_edge)){ //?
                    //     Q.push(elem1);
                    // }
                    if(checkcore_compress_IIIII(d,sim,u,elem1,a.compress[elem1],e_id_to_edge)){ //?
                        Q.push(elem1);
                    }
                }
       
            }
            n_sim.clear();
        }
        
       // N_SIM.clear();
        Cluster_ID++;
       
    }


    if(temp_sim<0)//
    {
        for(auto x:a.Cluster_Index[i-1]){
            //if(compress[x.second].empty()) continue;
            if(-x.count<u) break;
            if(edge_to_cluster[x.eid]!=0) continue;
            if(!checkcore_compress_IIIII(d,sim,u,x.eid,a.compress[x.eid],e_id_to_edge)) continue;   

            edge_to_cluster[x.eid]=Cluster_ID;
            Cluster[Cluster_ID-1].push_back(x.eid); //??????
            std::queue<size_t> Q;
            Q.push(x.eid);
            while(!Q.empty()){
                size_t temp=Q.front();
                Q.pop();
                //计算N_SIM_V
                compute_compress_II_n_sim(d,sim,temp,a.compress[temp],n_sim,e_id_to_edge);//瓶颈？
                for ( auto elem1:n_sim) {
                    //double xiangsidu=geometricMean(sim,e_id_to_edge,,,,);
                    if(edge_to_cluster[elem1]!=Cluster_ID){//还没被加入过当前cluster
                        
                        edge_to_cluster[elem1]=Cluster_ID;
                        Cluster[Cluster_ID-1].push_back(elem1);  //??????    直接在这里定义一个局部变量  最后调用函数将结果保存到csv文件中
                    
                        if(checkcore_compress_IIIII(d,sim,u,elem1,a.compress[elem1],e_id_to_edge)){ //?
                            Q.push(elem1);
                        }
                    }
        
                }
                n_sim.clear();
            }
            
            // N_SIM.clear();
            Cluster_ID++;
       
        }
    }

    end = clock();

    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.write_others(Cluster_ID-1);
    a.writecluster("../output/",Cluster,sim,u);
    std::cout<<Cluster_ID-1<<'\n';

}


void compress_IIII_cluster_bianjie(double d,double d_segment_count,std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u)
{
    a.output["algo"] = "optimize_compress_IIII";
    clock_t start, end;

    start = clock();
    /**********处理过程*********/
    intvec edge_to_cluster(e_id_to_edge.size(),0);//判断超边是否已经被加入过当前cluster
    intintvec Cluster(e_id_to_edge.size(),intvec{});//保存cluster 
    //  还是要把cluster里面超边的连接信息保存下来，不能只保存一个超边编号，会导致只知道他们在一个cluster里面而不知道他们的连接关系(实际上是可以知道的 不过要凭超边编号去原图中找连接关系)
    // 和普通图还是区别很大的 普通图中没有一条边和有一条边是不同的子图   超图中超边被加入cluster之后，由于共享节点的关系 不止会和引入它的超边有连接关系
    size_t Cluster_ID=1;
    intvec n_sim;//聚类时临时存储 与i相似超边编号，在每次结束时记得先把里面的内容清空

    //判断一下sim值是不是能够被d_segment_count整除
    double temp_sim=sim;
    unsigned int i=0;
    double d1[11]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}; //
    //std::map<size_t,int> cluster_i_jia_1={};
    while (temp_sim>d1[i])//  这里要不要也改成打表的形式  避免多次与d相减？    !=0
    {
        i+=1;
        if(i>=9) break;
    }

    // for(auto x:a.Cluster_Index[i-1])
    // {
    //     cluster_i_jia_1[x.second]=-x.first;
    // }
    unsigned geshu=(unsigned int)(1.0/d_segment_count);
   //double temp2=0;
    for(auto x:a.Cluster_Index[i]){
        //if(compress[x.second].empty()) continue;
        if(-x.count<u) break;
        if(edge_to_cluster[x.eid]!=0) continue;
        
        edge_to_cluster[x.eid]=Cluster_ID;
        Cluster[Cluster_ID-1].push_back(x.eid); //??????
        std::queue<size_t> Q;
        Q.push(x.eid);
        while(!Q.empty()){
            size_t temp=Q.front();
            Q.pop();
            //计算N_SIM_V
            //compute_compress_II_n_sim(d,sim,temp,a.compress[temp],n_sim,e_id_to_edge);//瓶颈在于这里core节点的邻居扩散？
            for(unsigned k=0;k<geshu-i;k++)
            {
              //  temp2=d1[9-k];// 9->geshu
                //if(temp2<d1[i]) break; // sim最大0.9
                    for ( auto elem1:a.compress[temp][k]) {
                        //double xiangsidu=geometricMean(sim,e_id_to_edge,,,,);
                        if(edge_to_cluster[elem1]!=Cluster_ID){//还没被加入过当前cluster
                            
                            edge_to_cluster[elem1]=Cluster_ID;
                            Cluster[Cluster_ID-1].push_back(elem1);  //??????    直接在这里定义一个局部变量  最后调用函数将结果保存到csv文件中
                            unsigned int count_sim=0;
                            for(unsigned int w=0;w<geshu-i;w++)
                            {
                                // temp2=d1[9-w];
                                count_sim+=a.compress[elem1][w].size();
                                if(count_sim>=u) break; // temp2+d-sim<epision ?
                            }
                            if(count_sim>=u){ //?
                                Q.push(elem1);
                            }
                        }
            
                    }
                
            }
           // n_sim.clear();
        }
        
       // N_SIM.clear();
        Cluster_ID++;
       
    }


    end = clock();

    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.write_others(Cluster_ID-1);
    a.writecluster("../output/",Cluster,sim,u);
    std::cout<<Cluster_ID-1<<'\n';

}

void updata_insert( compress_index &a,size_t node,size_t num_hyperegde)
{
    //首先判断是不是新顶点 新超边  之后更新LI如果需要   接着按照推理11更新si ci  
    //insert一条超边需要改(a.hg).hyperedges  .init_nodes  .node_index
    //还需要改li si ci
    clock_t start,end;
    start=clock();
    double d1[10]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
    bool is_new_edge=false;
    bool is_new_node=false;
    std::set<size_t> old_neighbors;
    std::set<size_t> new_neighbors; 
    intvec insert_hyperedge;
    intvec insert_table;
    size_t size_edge=(a.hg.hyperedges).size();
    if(num_hyperegde>=size_edge) 
    {
        is_new_edge=true;//大于size_edge时 新加入超边的编号也是size_edge
        insert_hyperedge.push_back(node);
        (a.hg.hyperedges).push_back(insert_hyperedge);
    }
    else
    {
        generate_neighbors_by_tables2(num_hyperegde,(a.hg).hyperedges,(a.hg).node_index,a,old_neighbors);
        (a.hg).hyperedges[num_hyperegde].push_back(node);// 随机插入的顶点如果已经存在 随机插入的超边中了呢？
    }

    size_t size_node=((a.hg).init_nodes).size();
    if(((a.hg).node_index).find(node)==((a.hg).node_index).end())
    {
        is_new_node=true;
        ((a.hg).init_nodes).push_back(node);
        (a.hg).node_index[node]=size_node;
    }
    
    //更新LI
    if(is_new_node)
    {
        if(is_new_edge) insert_table.push_back(size_edge);
        else insert_table.push_back(num_hyperegde);
        (a.tables).push_back(insert_table);
    }
    else
    {
        auto j=(a.hg).node_index[node];
        if(is_new_edge) a.tables[j].push_back(size_edge);
        else a.tables[j].push_back(num_hyperegde);
    }
    if(is_new_edge) generate_neighbors_by_tables2(size_edge,(a.hg).hyperedges,(a.hg).node_index,a,new_neighbors);
    else generate_neighbors_by_tables2(num_hyperegde,(a.hg).hyperedges,(a.hg).node_index,a,new_neighbors);
    
    //更新si CI
    double old_xiangsidu=0;
    double new_xiangsidu=0;
    size_t i=0;
    if(is_new_edge) i=size_edge;
    else i=num_hyperegde;
    std::vector<std::vector<size_t>> compress(10,std::vector<size_t>{});
    std::vector<core_pair3> Neighbor;
    // if(is_new_edge)
    // {
        double t1=0;
        unsigned int count_1=0;
        unsigned int count_0=0;
        unsigned int pos=0;
        double xiangsidu=0;
        core_pair3 cp3;
        core_pair2 cp22;
       // std::set<size_t>::iterator it;
    // }
    // else
    // {
       // std::cout<<new_neighbors.size()<<'\n';
        for(auto v:new_neighbors)
        {
            auto it=old_neighbors.find(v);
            if(old_neighbors.empty()||it==old_neighbors.end()) //是新邻居
            {
                old_xiangsidu=0;
                double  temp= sqrt(static_cast<double>(((a.hg).hyperedges[i]).size()) * ((a.hg).hyperedges[v]).size());
                new_xiangsidu=1.0/temp;
                
                count_1=0;
                t1=new_xiangsidu;
                while(t1>=d1[count_1]) //t1-d>=epision?t>=d 改成这样会不会边界情况好了 但是其他情况错了？ 在边界时候还是没有好
                {
                    count_1+=1;
                // t1=t1-d;
                if(count_1>=9) break;  //9也要变成和tau相对应的
                }                     //这里算错了导致存的位置不对 sim=  0.6   0.8的时候只有一个能对
                if(count_1>9) count_1=9;
                // if(is_new_edge) {

                // }
                // else a.compress[i][9-count_1].push_back();
                pos=0;
                for(auto y:a.compress[v][9-count_1])
                {
                    //按顺序插入 所以每个的相似度都要计算
                    xiangsidu=geometricMean_compress_II((a.hg).hyperedges,(a.hg).hyperedges[v].size(),(a.hg).hyperedges[y].size(),v,y);
                    if(xiangsidu>new_xiangsidu) pos+=1;
                    else break;
                }
                if(pos<(a.compress[v][9-count_1]).size()) (a.compress[v][9-count_1]).insert((a.compress[v][9-count_1]).begin()+pos,i);
                else (a.compress[v][9-count_1]).push_back(i);
                for(size_t k=0;k<=count_1;k++)
                {
                    for(size_t h=0;h<a.Cluster_Index[k].size();h++)
                    {
                        if(a.Cluster_Index[k][h].eid==v) {a.Cluster_Index[k][h].count+=1;break;}
                        if(h==a.Cluster_Index[k].size()-1) {cp22.eid=v;cp22.count=1;a.Cluster_Index[k].push_back(cp22);}
                    }
                }
            }
            else//不是新邻居
            {
               // std::cout<<"111"<<'\n';
                new_xiangsidu=geometricMean_compress_II((a.hg).hyperedges,(a.hg).hyperedges[i].size(),(a.hg).hyperedges[v].size(),i,v);
                if(std::find(((a.hg).hyperedges[v]).begin(),((a.hg).hyperedges[v]).end(),node)!=((a.hg).hyperedges[v]).end())
                {
                    double temp1=sqrt(static_cast<double>(((a.hg).hyperedges[i]).size()-1.0) * ((a.hg).hyperedges[v]).size());
                    double temp2=sqrt(static_cast<double>(((a.hg).hyperedges[i]).size()) * ((a.hg).hyperedges[v]).size())-1.0;
                    old_xiangsidu=new_xiangsidu*temp2/temp1;
                }
                else
                {
                    double temp1=sqrt(static_cast<double>(((a.hg).hyperedges[i]).size())-1.0);
                    double temp2=sqrt(static_cast<double>(((a.hg).hyperedges[i]).size()));
                    old_xiangsidu=new_xiangsidu*temp2/temp1;
                }
                //更新v的si 以及ci
                count_1=0;
                t1=new_xiangsidu;
                while(t1>=d1[count_1]) //t1-d>=epision?t>=d 改成这样会不会边界情况好了 但是其他情况错了？ 在边界时候还是没有好
                {
                    count_1+=1;
                // t1=t1-d;
                if(count_1>=9) break;  //9也要变成和tau相对应的
                }                     //这里算错了导致存的位置不对 sim=  0.6   0.8的时候只有一个能对
                if(count_1>9) count_1=9;

                count_0=0;
                t1=old_xiangsidu;
                while(t1>=d1[count_0]) //t1-d>=epision?t>=d 改成这样会不会边界情况好了 但是其他情况错了？ 在边界时候还是没有好
                {
                    count_0+=1;
                // t1=t1-d;
                if(count_0>=9) break;  //9也要变成和tau相对应的
                }                     //这里算错了导致存的位置不对 sim=  0.6   0.8的时候只有一个能对
                if(count_0>9) count_0=9;

                pos=0;
                for(auto y:a.compress[v][9-count_0])
                {
                    //按顺序插入 所以每个的相似度都要计算
                    //xiangsidu=geometricMean_compress_II((a.hg).hyperedges,(a.hg).hyperedges[v].size(),(a.hg).hyperedges[y].size(),v,y);
                    //if(xiangsidu>new_xiangsidu) pos+=1;
                    if(y!=i) pos+=1;
                    else break;
                }
                if(!a.compress[v][9-count_0].empty()){
                    if(pos==a.compress[v][9-count_0].size()) pos--;
                    if(a.compress[v][9-count_0][pos]==i) (a.compress[v][9-count_0]).erase((a.compress[v][9-count_0]).begin()+pos);
                }
                //std::cout<<"222"<<'\n';
                pos=0;
                for(auto y:a.compress[v][9-count_1])
                {
                    //按顺序插入 所以每个的相似度都要计算
                    xiangsidu=geometricMean_compress_II((a.hg).hyperedges,(a.hg).hyperedges[v].size(),(a.hg).hyperedges[y].size(),v,y);
                    if(xiangsidu>new_xiangsidu) pos+=1;
                    else break;
                }
                //(a.compress[v][9-count_1]).insert((a.compress[v][9-count_1]).begin()+pos,i);
                if(pos<(a.compress[v][9-count_1]).size()) (a.compress[v][9-count_1]).insert((a.compress[v][9-count_1]).begin()+pos,i);
                else (a.compress[v][9-count_1]).push_back(i);
                if(count_1>count_0)
                {
                    for(size_t k=count_0+1;k<=count_1;k++)
                    {
                        for(size_t h=0;h<a.Cluster_Index[k].size();h++)
                        {
                            if(a.Cluster_Index[k][h].eid==v) {a.Cluster_Index[k][h].count+=1;break;}
                            if(h==a.Cluster_Index[k].size()-1) {cp22.eid=v;cp22.count=1;a.Cluster_Index[k].push_back(cp22);}
                        }
                    }
                }
                if(count_1<count_0)
                {
                    for(size_t k=count_1+1;k<=count_0;k++)
                    {
                        for(size_t h=0;h<a.Cluster_Index[k].size();h++)
                        {
                            if(a.Cluster_Index[k][h].eid==v) {a.Cluster_Index[k][h].count-=1;break;}
                            //if(h==a.Cluster_Index[k].size()-1) {cp22.eid=v;cp22.count=1;a.Cluster_Index[k].push_back(cp22);}
                        }
                    }
                }
                
            }
           // std::cout<<"333"<<'\n';
             if(i!=v){
                // Neighbor.insert(std::make_pair(-xiangsidu,v));  //1-
                cp3.eid=v;
                cp3.xiangsidu=-new_xiangsidu;
                Neighbor.push_back(cp3);  //1-
             }
            
            
        }
        sort(Neighbor.begin(),Neighbor.end(),comp3);
        
        for(auto p:Neighbor)
        {
            count_1=0;
            t1=-p.xiangsidu;
            while(t1>=d1[count_1]) //t1-d>=epision?t>=d 改成这样会不会边界情况好了 但是其他情况错了？ 在边界时候还是没有好
            {
                count_1+=1;
               // t1=t1-d;
               if(count_1>=9) break;  //9也要变成和tau相对应的
            }                     //这里算错了导致存的位置不对 sim=  0.6   0.8的时候只有一个能对
            if(count_1>9) count_1=9;
            compress[9-count_1].push_back(p.eid);
        }

        if(is_new_edge) (a.compress).push_back(compress);
        else 
        {
            (a.compress)[i].clear();//有无问题
            for(auto y:compress)
            {
                (a.compress[i]).push_back(y);
            }
        }
        
   // }
      //更新i的ci   以及所有ci排序
      size_t sum=0;
      
      for(size_t k=0;k<10;k++)
      {
            sum+=a.compress[9-k].size();
            if(sum==0) continue;
            if(!is_new_edge){
                for(size_t h=0;h<a.Cluster_Index[k].size();h++)
                {
                    if(a.Cluster_Index[k][h].eid==i) {a.Cluster_Index[k][h].count=sum;break;}
                    if(h==a.Cluster_Index[k].size()-1) {cp22.eid=i;cp22.count=sum;a.Cluster_Index[k].push_back(cp22);}
                }
            }
            else {cp22.eid=i;cp22.count=sum;a.Cluster_Index[k].push_back(cp22);}
      }


    for(unsigned int k=0;k<10;k++) {
        sort(a.Cluster_Index[k].begin(),a.Cluster_Index[k].end(),comp2);
    }

    end=clock();
    a.updata_insert_time=double(end - start) / double(CLOCKS_PER_SEC);
    // size_t max_u=0;
    // for(size_t i=0;i<a.neighbor_order.size();i++)
    // {
    //     if(a.neighbor_order[i].size()>max_u) max_u=a.neighbor_order[i].size();
    // }
    // a.max_neighbor=max_u;
    std::cout<<"insert time  "<<a.updata_insert_time<<'\n';
}
void updata_remove( compress_index &a,size_t num_node,size_t num_hyperegde)
{
    clock_t start,end;
    start=clock();
    double old_xiangsidu=0;
    double new_xiangsidu=0;
    std::vector<std::vector<size_t>> compress(10,std::vector<size_t>{});
    std::vector<core_pair3> Neighbor;
    double t1=0;
    unsigned int count_1=0;
    unsigned int count_0=0;
    double xiangsidu=0;
    core_pair3 cp3;
    core_pair2 cp22;
    std::set<size_t> old_neighbors;
    std::set<size_t> new_neighbors; 
    double d1[10]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};

    generate_neighbors_by_tables2(num_hyperegde,(a.hg).hyperedges,(a.hg).node_index,a,old_neighbors);
    size_t node;
    if((a.hg).hyperedges[num_hyperegde].size()!=0){
        node=(a.hg).hyperedges[num_hyperegde][num_node];
        (a.hg).hyperedges[num_hyperegde].erase((a.hg).hyperedges[num_hyperegde].begin()+num_node);
    }
    auto j=(a.hg).node_index[node];
    size_t pos=0;
    for(size_t k=0;k<a.tables[j].size();k++)
    {
        if(a.tables[j][k]!=num_hyperegde) pos+=1;
        else break;
    }
    if(!a.tables[j].empty())
    {
        if(pos==a.tables[j].size()) pos-=1;
        if(a.tables[j][pos]==num_hyperegde) {a.tables[j].erase(a.tables[j].begin()+pos);}
    }
    if((a.hg).hyperedges[num_hyperegde].size()!=0) generate_neighbors_by_tables2(num_hyperegde,(a.hg).hyperedges,(a.hg).node_index,a,new_neighbors);
    size_t i=num_hyperegde;
    for(auto v:old_neighbors)
    {
        // if((a.hg).hyperedges[num_hyperegde].size()==0)
        // {
        //     //单独处理删除后为空
        // }
        // else
        // {
        //     generate_neighbors_by_tables2(num_hyperegde,(a.hg).hyperedges,(a.hg).node_index,a,new_neighbors);
        //     //三种情况处理

        // }
            auto it=new_neighbors.find(v);
            if(new_neighbors.empty()||it==new_neighbors.end()) //
            {
                new_xiangsidu=0;
                double  temp= sqrt(static_cast<double>(((a.hg).hyperedges[i]).size()+1) * ((a.hg).hyperedges[v]).size());
                old_xiangsidu=1.0/temp;

                count_1=0;
                t1=old_xiangsidu;
                while(t1>=d1[count_1]) //t1-d>=epision?t>=d 改成这样会不会边界情况好了 但是其他情况错了？ 在边界时候还是没有好
                {
                    count_1+=1;
                // t1=t1-d;
                if(count_1>=9) break;  //9也要变成和tau相对应的
                }                     //这里算错了导致存的位置不对 sim=  0.6   0.8的时候只有一个能对
                if(count_1>9) count_1=9;
                // if(is_new_edge) {

                // }
                // else a.compress[i][9-count_1].push_back();
                pos=0;
                for(auto y:a.compress[v][9-count_1])
                {
                    //按顺序插入 所以每个的相似度都要计算
                   // xiangsidu=geometricMean_compress_II((a.hg).hyperedges,(a.hg).hyperedges[v].size(),(a.hg).hyperedges[y].size(),v,y);
                    if(y!=i) pos+=1;
                    else break;
                }
                // if(pos<(a.compress[v][9-count_1]).size()) (a.compress[v][9-count_1]).insert((a.compress[v][9-count_1]).begin()+pos,i);
                // else (a.compress[v][9-count_1]).push_back(i);
                if(pos<(a.compress[v][9-count_1]).size()) (a.compress[v][9-count_1]).erase((a.compress[v][9-count_1]).begin()+pos);
                //else {(a.compress[v][9-count_1]).erase((a.compress[v][9-count_1]).begin()+pos);}
                for(size_t k=0;k<=count_1;k++)
                {
                    for(size_t h=0;h<a.Cluster_Index[k].size();h++)
                    {
                        if(a.Cluster_Index[k][h].eid==v) {a.Cluster_Index[k][h].count-=1;break;}
                        //if(h==a.Cluster_Index[k].size()-1) {cp22.eid=v;cp22.count=1;a.Cluster_Index[k].push_back(cp22);}
                    }
                }
            }

            else
            {
                new_xiangsidu=geometricMean_compress_II((a.hg).hyperedges,(a.hg).hyperedges[i].size(),(a.hg).hyperedges[v].size(),i,v);
                if(std::find(((a.hg).hyperedges[v]).begin(),((a.hg).hyperedges[v]).end(),node)!=((a.hg).hyperedges[v]).end())
                {
                    double temp1=sqrt(static_cast<double>(((a.hg).hyperedges[i]).size()+1.0) * ((a.hg).hyperedges[v]).size());
                    double temp2=sqrt(static_cast<double>(((a.hg).hyperedges[i]).size()) * ((a.hg).hyperedges[v]).size())+1.0;
                    old_xiangsidu=new_xiangsidu*temp2/temp1;
                }
                else
                {
                    double temp1=sqrt(static_cast<double>(((a.hg).hyperedges[i]).size())+1.0);
                    double temp2=sqrt(static_cast<double>(((a.hg).hyperedges[i]).size()));
                    old_xiangsidu=new_xiangsidu*temp2/temp1;
                }
            

                count_1=0;
                t1=new_xiangsidu;
                while(t1>=d1[count_1]) //t1-d>=epision?t>=d 改成这样会不会边界情况好了 但是其他情况错了？ 在边界时候还是没有好
                {
                    count_1+=1;
                // t1=t1-d;
                if(count_1>=9) break;  //9也要变成和tau相对应的
                }                     //这里算错了导致存的位置不对 sim=  0.6   0.8的时候只有一个能对
                if(count_1>9) count_1=9;

                count_0=0;
                t1=old_xiangsidu;
                while(t1>=d1[count_0]) //t1-d>=epision?t>=d 改成这样会不会边界情况好了 但是其他情况错了？ 在边界时候还是没有好
                {
                    count_0+=1;
                // t1=t1-d;
                if(count_0>=9) break;  //9也要变成和tau相对应的
                }                     //这里算错了导致存的位置不对 sim=  0.6   0.8的时候只有一个能对
                if(count_0>9) count_0=9;

                pos=0;
                for(auto y:a.compress[v][9-count_0])
                {
                    //按顺序插入 所以每个的相似度都要计算
                    //xiangsidu=geometricMean_compress_II((a.hg).hyperedges,(a.hg).hyperedges[v].size(),(a.hg).hyperedges[y].size(),v,y);
                    //if(xiangsidu>new_xiangsidu) pos+=1;
                    if(y!=i) pos+=1;
                    else break;
                }
                // if(!a.compress[v][9-count_0].empty()){
                //     if(pos==a.compress[v][9-count_0].size()) pos--;
                //     if(a.compress[v][9-count_0][pos]==i) (a.compress[v][9-count_0]).erase((a.compress[v][9-count_0]).begin()+pos);
                // }
                if(pos<(a.compress[v][9-count_0]).size()) (a.compress[v][9-count_0]).erase((a.compress[v][9-count_0]).begin()+pos);
                //std::cout<<"222"<<'\n';
                pos=0;
                for(auto y:a.compress[v][9-count_1])
                {
                    //按顺序插入 所以每个的相似度都要计算
                    xiangsidu=geometricMean_compress_II((a.hg).hyperedges,(a.hg).hyperedges[v].size(),(a.hg).hyperedges[y].size(),v,y);
                    if(xiangsidu>new_xiangsidu) pos+=1;
                    else break;
                }
                //(a.compress[v][9-count_1]).insert((a.compress[v][9-count_1]).begin()+pos,i);
                if(pos<(a.compress[v][9-count_1]).size()) (a.compress[v][9-count_1]).insert((a.compress[v][9-count_1]).begin()+pos,i);
                else (a.compress[v][9-count_1]).push_back(i);
                if(count_1>count_0)
                {
                    for(size_t k=count_0+1;k<=count_1;k++)
                    {
                        for(size_t h=0;h<a.Cluster_Index[k].size();h++)
                        {
                            if(a.Cluster_Index[k][h].eid==v) {a.Cluster_Index[k][h].count+=1;break;}
                            if(h==a.Cluster_Index[k].size()-1) {cp22.eid=v;cp22.count=1;a.Cluster_Index[k].push_back(cp22);}
                        }
                    }
                }
                if(count_1<count_0)
                {
                    for(size_t k=count_1+1;k<=count_0;k++)
                    {
                        for(size_t h=0;h<a.Cluster_Index[k].size();h++)
                        {
                            if(a.Cluster_Index[k][h].eid==v) {a.Cluster_Index[k][h].count-=1;break;}
                           // if(h==a.Cluster_Index[k].size()-1) {cp22.eid=v;cp22.count=1;a.Cluster_Index[k].push_back(cp22);}
                        }
                    }
                }
                    // std::cout<<"333"<<'\n';
                if(i!=v){
                    // Neighbor.insert(std::make_pair(-xiangsidu,v));  //1-
                    cp3.eid=v;
                    cp3.xiangsidu=-new_xiangsidu;
                    Neighbor.push_back(cp3);  //1-
                }
            }
           
    }
        sort(Neighbor.begin(),Neighbor.end(),comp3);
        
        for(auto p:Neighbor)
        {
            count_1=0;
            t1=-p.xiangsidu;
            while(t1>=d1[count_1]) //t1-d>=epision?t>=d 改成这样会不会边界情况好了 但是其他情况错了？ 在边界时候还是没有好
            {
                count_1+=1;
               // t1=t1-d;
               if(count_1>=9) break;  //9也要变成和tau相对应的
            }                     //这里算错了导致存的位置不对 sim=  0.6   0.8的时候只有一个能对
            if(count_1>9) count_1=9;
            compress[9-count_1].push_back(p.eid);
        }

        // if(is_new_edge) (a.compress).push_back(compress);
        // else 
        // {
            (a.compress)[i].clear();//有无问题
            for(auto y:compress)
            {
                (a.compress[i]).push_back(y);
            }
       // }
        
    // }
      //更新i的ci   以及所有ci排序
      size_t sum=0;
      
      for(size_t k=0;k<10;k++)
      {
            sum+=a.compress[9-k].size();
            if(sum==0) continue;
           // if(!is_new_edge){
                for(size_t h=0;h<a.Cluster_Index[k].size();h++)
                {
                    if(a.Cluster_Index[k][h].eid==i) {a.Cluster_Index[k][h].count=sum;break;}
                    if(h==a.Cluster_Index[k].size()-1) {cp22.eid=i;cp22.count=sum;a.Cluster_Index[k].push_back(cp22);}
                }
           // }
           // else {cp22.eid=i;cp22.count=sum;a.Cluster_Index[k].push_back(cp22);}
      }


    for(unsigned int k=0;k<10;k++) {
        sort(a.Cluster_Index[k].begin(),a.Cluster_Index[k].end(),comp2);
    }


    end=clock();
    a.updata_remove_time=double(end - start) / double(CLOCKS_PER_SEC);
    std::cout<<"remove time  "<<a.updata_remove_time<<'\n';
}


void compress_IIII_cluster_hash(double d,double d_segment_count,std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u)
{
    a.output["algo"] = "optimize_compress_IIII";
    clock_t start, end;

    start = clock();
    /**********处理过程*********/
    intvec edge_to_cluster(e_id_to_edge.size(),0);//判断超边是否已经被加入过当前cluster
    intintvec Cluster(e_id_to_edge.size(),intvec{});//保存cluster 
    //  还是要把cluster里面超边的连接信息保存下来，不能只保存一个超边编号，会导致只知道他们在一个cluster里面而不知道他们的连接关系(实际上是可以知道的 不过要凭超边编号去原图中找连接关系)
    // 和普通图还是区别很大的 普通图中没有一条边和有一条边是不同的子图   超图中超边被加入cluster之后，由于共享节点的关系 不止会和引入它的超边有连接关系
    size_t Cluster_ID=1;
    intvec n_sim;//聚类时临时存储 与i相似超边编号，在每次结束时记得先把里面的内容清空
    intvec core_table(e_id_to_edge.size(),0);//记录是不是core
    //判断一下sim值是不是能够被d_segment_count整除
    double temp_sim=sim;
    unsigned int i=0;
    //std::map<size_t,int> cluster_i_jia_1={};
    double d1[11]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
    while (temp_sim>d1[i])// 这里要不要也改成打表的形式  避免多次与d相减？  只考虑非边界情况就可以了  这里好像阴差阳错都可以
    {
        i+=1;
        if(temp_sim<d1[i]) {temp_sim=-1;break;}
        if(i>=9) break;
    }
    
    for(auto x:a.Cluster_Index[i])
    {
        if(-x.count<u) break;
        else core_table[x.eid]=1;
    }


    for(auto x:a.Cluster_Index[i-1])
    {
        if(-x.count<u) break;
        if(core_table[x.eid]==1) continue;
        if(checkcore_compress_IIIII(d,sim,u,x.eid,a.compress[x.eid],e_id_to_edge)) core_table[x.eid]=1;
    }


    for(auto x:a.Cluster_Index[i]){ //改Cluster index的数据结构和gs index一样
        //if(compress[x.second].empty()) continue;
        if(-x.count<u) break;
        if(edge_to_cluster[x.eid]!=0) continue;
        if(temp_sim>0){//epision而不是0，为了防止边界情况temo_sim是一个大于0的极小值的情况
            if(core_table[x.eid]!=1) continue;   
        }//对于sim大于0.9的时候专门设置的
        edge_to_cluster[x.eid]=Cluster_ID;
        Cluster[Cluster_ID-1].push_back(x.eid); //??????
        std::queue<size_t> Q;
        Q.push(x.eid);
        while(!Q.empty()){
            size_t temp=Q.front();
            Q.pop();
            //计算N_SIM_V
            compute_compress_II_n_sim(d,sim,temp,a.compress[temp],n_sim,e_id_to_edge);//瓶颈在于这里core节点的邻居扩散？
            for ( auto elem1:n_sim) {
                //double xiangsidu=geometricMean(sim,e_id_to_edge,,,,);
                if(edge_to_cluster[elem1]!=Cluster_ID){//还没被加入过当前cluster
                    
                    edge_to_cluster[elem1]=Cluster_ID;
                    Cluster[Cluster_ID-1].push_back(elem1);  //??????    直接在这里定义一个局部变量  最后调用函数将结果保存到csv文件中
                    // if(checkcore_compress_IIIII(d,i,cluster_i_jia_1,d_segment_count,sim,u,elem1,a.compress[elem1],e_id_to_edge)){ //?
                    //     Q.push(elem1);
                    // }
                    if(core_table[elem1]==1){ //?
                        Q.push(elem1);
                    }
                }
       
            }
            n_sim.clear();
        }
        
       // N_SIM.clear();
        Cluster_ID++;
       
    }


    if(temp_sim<0)//
    {
        for(auto x:a.Cluster_Index[i-1]){
            //if(compress[x.second].empty()) continue;
            if(-x.count<u) break;
            if(edge_to_cluster[x.eid]!=0) continue;
            if(core_table[x.eid]!=1) continue;   

            edge_to_cluster[x.eid]=Cluster_ID;
            Cluster[Cluster_ID-1].push_back(x.eid); //??????
            std::queue<size_t> Q;
            Q.push(x.eid);
            while(!Q.empty()){
                size_t temp=Q.front();
                Q.pop();
                //计算N_SIM_V
                compute_compress_II_n_sim(d,sim,temp,a.compress[temp],n_sim,e_id_to_edge);//瓶颈？
                for ( auto elem1:n_sim) {
                    //double xiangsidu=geometricMean(sim,e_id_to_edge,,,,);
                    if(edge_to_cluster[elem1]!=Cluster_ID){//还没被加入过当前cluster
                        
                        edge_to_cluster[elem1]=Cluster_ID;
                        Cluster[Cluster_ID-1].push_back(elem1);  //??????    直接在这里定义一个局部变量  最后调用函数将结果保存到csv文件中
                    
                        if(core_table[elem1]==1){ //?
                            Q.push(elem1);
                        }
                    }
        
                }
                n_sim.clear();
            }
            
            // N_SIM.clear();
            Cluster_ID++;
       
        }
    }

    end = clock();

    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.write_others(Cluster_ID-1);
    a.writecluster("../output/",Cluster,sim,u);
    std::cout<<Cluster_ID-1<<'\n';

}


void compress_IIII_cluster_bianjie_hash(double d,double d_segment_count,std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u)
{
    a.output["algo"] = "optimize_compress_IIII";
    clock_t start, end;

    start = clock();
    /**********处理过程*********/
    intvec edge_to_cluster(e_id_to_edge.size(),0);//判断超边是否已经被加入过当前cluster
    intintvec Cluster(e_id_to_edge.size(),intvec{});//保存cluster 
    //  还是要把cluster里面超边的连接信息保存下来，不能只保存一个超边编号，会导致只知道他们在一个cluster里面而不知道他们的连接关系(实际上是可以知道的 不过要凭超边编号去原图中找连接关系)
    // 和普通图还是区别很大的 普通图中没有一条边和有一条边是不同的子图   超图中超边被加入cluster之后，由于共享节点的关系 不止会和引入它的超边有连接关系
    size_t Cluster_ID=1;
    intvec n_sim;//聚类时临时存储 与i相似超边编号，在每次结束时记得先把里面的内容清空
    intvec core_table(e_id_to_edge.size(),0);
    //判断一下sim值是不是能够被d_segment_count整除
    double temp_sim=sim;
    unsigned int i=0;
    double d1[11]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
    unsigned geshu=(unsigned int)(1.0/d_segment_count);
    //std::map<size_t,int> cluster_i_jia_1={};
    while (temp_sim>d1[i])//  这里要不要也改成打表的形式  避免多次与d相减？    !=0
    {
        i+=1;
        if(i>=9) break;
    }
    
    for(auto x:a.Cluster_Index[i])
    {
        if(-x.count<u) break;
        unsigned int count_s=0;
        for(unsigned int w=0;w<geshu-i;w++)
        {
            // temp2=d1[9-w];
            count_s+=a.compress[x.eid][w].size();
            if(count_s>=u) break; // temp2+d-sim<epision ?
        }
        if(count_s>=u) core_table[x.eid]=1;
    }
    // for(auto x:a.Cluster_Index[i-1])
    // {
    //     cluster_i_jia_1[x.second]=-x.first;
    // }
    // unsigned geshu=(unsigned int)(1.0/d_segment_count);
   //double temp2=0;
    for(auto x:a.Cluster_Index[i]){
        //if(compress[x.second].empty()) continue;
        if(-x.count<u) break;
        if(edge_to_cluster[x.eid]!=0) continue;
        
        edge_to_cluster[x.eid]=Cluster_ID;
        Cluster[Cluster_ID-1].push_back(x.eid); //??????
        std::queue<size_t> Q;
        Q.push(x.eid);
        while(!Q.empty()){
            size_t temp=Q.front();
            Q.pop();
            //计算N_SIM_V
            //compute_compress_II_n_sim(d,sim,temp,a.compress[temp],n_sim,e_id_to_edge);//瓶颈在于这里core节点的邻居扩散？
            for(unsigned k=0;k<geshu-i;k++)
            {
              //  temp2=d1[9-k];// 9->geshu
                //if(temp2<d1[i]) break; // sim最大0.9
                    for ( auto elem1:a.compress[temp][k]) {
                        //double xiangsidu=geometricMean(sim,e_id_to_edge,,,,);
                        if(edge_to_cluster[elem1]!=Cluster_ID){//还没被加入过当前cluster
                            
                            edge_to_cluster[elem1]=Cluster_ID;
                            Cluster[Cluster_ID-1].push_back(elem1);  //??????    直接在这里定义一个局部变量  最后调用函数将结果保存到csv文件中
                            // unsigned int count_sim=0;
                            // for(unsigned int w=0;w<geshu-i;w++)
                            // {
                            //     // temp2=d1[9-w];
                            //     count_sim+=a.compress[elem1][w].size();
                            //     if(count_sim>=u) break; // temp2+d-sim<epision ?
                            // }
                            if(core_table[elem1]==1){ //?
                                Q.push(elem1);
                            }
                        }
            
                    }
                
            }
           // n_sim.clear();
        }
        
       // N_SIM.clear();
        Cluster_ID++;
       
    }


    end = clock();

    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.write_others(Cluster_ID-1);
    a.writecluster("../output/",Cluster,sim,u);
    std::cout<<Cluster_ID-1<<'\n';

}

//另外索引(加一个维度)
bool checkcore_compress_III(double d,std::vector<std::vector<unsigned int>> &segment_count,double d_segment_count,double sim,unsigned int u,size_t i,const intintvec &compress,const intintvec & Hyperedge)
{
    unsigned int sum=0;
    double temp_sim=sim;
    while(temp_sim>epision)
    {
        temp_sim-=d_segment_count;
        sum+=1;
    }
    if(segment_count[i][10-sum]<u) return false;
    if(sum<9){
        if(segment_count[i][8-sum]>=u) return true;
    }
    
    unsigned int count=0;
    for(unsigned k=0;k<(unsigned int)(1.0/d);k++)
    {
        if(count>=u) return true;
        double temp=(9-k)*d;
        if(temp+d-sim<epision) break;
        if(temp>=sim) {count+=compress[k].size();continue;}
        for(auto z:compress[k])
        {
            if(geometricMean2(sim,Hyperedge,Hyperedge[i].size(),Hyperedge[z].size(),i,z)) {
                count+=1;
                if(count>=u) return true;
            }
            else break;
            //if(count>=u) return true;
        }

    }
    if(count>=u) return true;
    else return false;
    
}
void compress_III_construct(double d,double d_segment_count,std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u)
{
    
    clock_t start;
    start = clock();
    /***************得到Nv邻域*******************/  
   
    douintset Neighbor;   //单个neighbor（包含相似度以及编号）
    std::set<size_t> computed_neighbor;//已经计算过了的邻居（不包含相似度）
    //compresstuplevecvec compress(e_id_to_edge.size(),compresstuplevec{});  //存储压缩之后的索引结构
    //size_t N = init_nodes.size();//顶点数
    //intintvec tables(N,intvec{});//索引表 有多少个顶点就建多少个表
   // std::vector<intintvec> compress(e_id_to_edge.size(),intintvec{});
    unsigned int geshu1=(unsigned int)(1.0/d);
    unsigned int geshu=(unsigned int)(1.0/d_segment_count);
    std::vector<std::vector<std::vector<size_t>>> compress(e_id_to_edge.size(),std::vector<std::vector<size_t>>(geshu1));
   // double  d=0.2; //划分间隙  固定区间划分还是根据第一个的相似度  固定长度划分?

   // double d_segment_count=0.1;  //分段统计间隙
    std::vector<std::vector<unsigned int>> segment_count(e_id_to_edge.size(),std::vector<unsigned int>{});//分段统计  这里很多零会导致用的空间变多


    std::cout<<"start!"<<'\n';
    tables_construct(e_id_to_edge,init_nodes,node_index,a);
    // for(size_t i=0;i<e_id_to_edge.size();i++)//建表
    // {
    //     for (const auto& elem:e_id_to_edge[i]) {
    //          auto j = node_index[elem];
    //          tables[j].push_back(i);
    //     }
    // }
   // std::cout<<"2!"<<'\n';
    for(size_t i=0;i<e_id_to_edge.size();i++)
    {

        generate_neighbors_by_tables2(i,e_id_to_edge,node_index,a,computed_neighbor);
        //generate_neighbors_by_enumerate2(i,e_id_to_edge,computed_neighbor);
        for(auto v:computed_neighbor)
        {
            double xiangsidu=geometricMean_compress_II(e_id_to_edge,e_id_to_edge[i].size(),e_id_to_edge[v].size(),i,v);
            if(i!=v){
                Neighbor.insert(std::make_pair(-xiangsidu,v));  //1-
            }
        }
        
        double t1=0;
        unsigned int count_1=0;

        for(auto p:Neighbor)
        {
            count_1=0;
            t1=-p.first;
            while(t1>=d)
            {
                count_1+=1;
                t1=t1-d;
            }
            if(count_1>9) count_1=9;
            compress[i][geshu1-count_1-1].push_back(p.second);
        }
        
        
        for(unsigned int k=0;k<geshu;k++) segment_count[i].push_back(0); 
        for(unsigned int k=0;k<geshu;k++)
        {
            segment_count[i][k]=compress[i][k].size();
        }
        for(unsigned int k=1;k<geshu;k++) segment_count[i][k]+=segment_count[i][k-1];// 总和而不是区间内部个数
        computed_neighbor.clear();
        Neighbor.clear();

    }
    std::cout<<"compute ok!"<<'\n';
    

    clock_t e_tm = clock();
    a.init_time=double(e_tm - start) / double(CLOCKS_PER_SEC);
    a.output["inition time"]= std::to_string(a.init_time);
    a.write_indexIII(compress,segment_count);
// 统计一下内存大小
    size_t  neicun=0;
    for(size_t i=0;i<a.tables.size();i++) {neicun+=a.tables[i].size()*8;}
    for(size_t i=0;i<compress.size();i++) {
        for(auto x:compress[i])
        {
            neicun+=x.size()*8;
        }
    }
    for(size_t i=0;i<segment_count.size();i++) {neicun+=segment_count[i].size()*4;}
    std::cout<<"neicun"<<' '<<neicun<<'\n';
    // std::cout<<segment_count[4714][0]<<' '<<segment_count[4714][1]<<' '<<segment_count[4714][2]<<' '<<segment_count[4714][3]<<' '<<segment_count[4714][4]<<' '<<'\n';
    // std::cout<<segment_count[4714][5]<<' '<<segment_count[4714][6]<<' '<<segment_count[4714][7]<<' '<<segment_count[4714][8]<<' '<<segment_count[4714][9]<<' '<<'\n';
    

}
void compress_III_cluster(double d,double d_segment_count,std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u)
{
    a.output["algo"] = "optimize_compress_III";
    clock_t start, end;
    start = clock();
    /**********处理过程*********/
    intvec edge_to_cluster(e_id_to_edge.size(),0);//判断超边是否已经被加入过当前cluster
    intintvec Cluster(e_id_to_edge.size(),intvec{});//保存cluster 
    //  还是要把cluster里面超边的连接信息保存下来，不能只保存一个超边编号，会导致只知道他们在一个cluster里面而不知道他们的连接关系(实际上是可以知道的 不过要凭超边编号去原图中找连接关系)
    // 和普通图还是区别很大的 普通图中没有一条边和有一条边是不同的子图   超图中超边被加入cluster之后，由于共享节点的关系 不止会和引入它的超边有连接关系
    size_t Cluster_ID=1;

    intvec n_sim;//聚类时临时存储 与i相似超边编号，在每次结束时记得先把里面的内容清空

    for(size_t i=0;i<e_id_to_edge.size();i++){
        if(a.compress[i].empty()) continue;
        if(edge_to_cluster[i]!=0) continue;

        if(!checkcore_compress_III(d,a.segment_count,d_segment_count,sim,u,i,a.compress[i],e_id_to_edge)) continue;   
        
        edge_to_cluster[i]=Cluster_ID;
        Cluster[Cluster_ID-1].push_back(i); //??????
        std::queue<size_t> Q;
        Q.push(i);
        while(!Q.empty()){
            size_t temp=Q.front();
            Q.pop();
            //计算N_SIM_V
            compute_compress_II_n_sim(d,sim,temp,a.compress[temp],n_sim,e_id_to_edge);
            for ( auto elem1:n_sim) {
                //double xiangsidu=geometricMean(sim,e_id_to_edge,,,,);
                if(edge_to_cluster[elem1]!=Cluster_ID){//还没被加入过当前cluster
                    
                    edge_to_cluster[elem1]=Cluster_ID;
                    Cluster[Cluster_ID-1].push_back(elem1);  //??????    直接在这里定义一个局部变量  最后调用函数将结果保存到csv文件中
                   
                    if(checkcore_compress_III(d,a.segment_count,d_segment_count,sim,u,elem1,a.compress[elem1],e_id_to_edge)){ //?
                        Q.push(elem1);
                    }
                }
       
            }
            n_sim.clear();
        }
        
       // N_SIM.clear();
        Cluster_ID++;
       
    }
    end = clock();

    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.write_others(Cluster_ID-1);
    a.writecluster("../output/",Cluster,sim,u);
    std::cout<<Cluster_ID-1<<'\n';

}