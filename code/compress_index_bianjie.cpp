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
#include "compress_index_bianjie.h"
#include "utils.h"
#define epision 1e-9
//--------------------------------------------------------------------- Utility functions ------------------------------------------------------------------------------

void compress_index::writecluster(std::string folder,intintvec Cluster){
    //std::cout << "core: \n";把core保存到csv文件中，core是运行过程中的变量，只存在于运行过程，我们要看到他的话要么print，要么写道csv文件里再看
    std::string file = folder + "cluster_"+output["algo"]+"_"+hg.dataset+".csv";// output是什么？ 好像就是打了个表 用关键词去对应另一个
    std::cout<<"writing to: "<<file<<"\n";
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
    std::string file1="cluster_"+output["algo"]+"_"+hg.dataset+".csv";
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
    double divisionResult = static_cast<double>(num1>num2?num2:num1) / result;

    if(divisionResult<sim) return false;

    if(((double)num1<sim*sim*num2)||((double)num2<sim*sim*num1)) return false;

    unsigned int count=0;
    divisionResult=0;
    for (unsigned int i=0;i<Hyperedge[e1].size();i++) {
             // if((static_cast<double>(count+Hyperedge[e1].size()-i) / result)<sim) break;
              auto it = std::find(Hyperedge[e2].begin(), Hyperedge[e2].end(), Hyperedge[e1][i]);
              if (it != Hyperedge[e2].end()) {
                count++;// 这里有提前终止的可能？
                divisionResult = static_cast<double>(count) / result;
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

    for(unsigned k=0;k<(unsigned int)(1.0/d);k++)
    {
        if(count>=u) return true;
        double temp=(9-k)*d;
        if(temp+d-sim<epision) break;//这里用小于0还是小于epision？
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

    for(unsigned k=0;k<(unsigned int)(1.0/d);k++)
    {
        double temp=(9-k)*d;
        if(temp+d-sim<epision) break;
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

    }
   
}
void compress_II_construct(double d,std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u)
{
    
    clock_t start;
    start = clock();
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
            while(t>=d)//
            {
                count_1+=1;
                t=t-d;
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
    a.writecluster("../output/",Cluster);
    std::cout<<Cluster_ID-1<<'\n';
}

//改进压缩索引，增加cluster_index
bool comp2(core_pair2 &a,core_pair2 &b)
{
    if(a.count!=b.count) return a.count<b.count;
    else return a.eid<b.eid;
}
bool checkcore_compress_IIII(double d,unsigned int k, std::map<size_t,int>&Cluster_i_jia_1,double d_segment_count,double sim,unsigned int u,size_t i,const intintvec &compress,const intintvec & Hyperedge)
{
    //使用map时
    // if(k<9){
    //     auto it=Cluster_Index[k+1].find(i);
    //     if(it!=Cluster_Index[k+1].end()&&it->second>=u) return true;
    // }
    
    //使用set时
    
    auto it=Cluster_i_jia_1.find(i);
    if(it!=Cluster_i_jia_1.end()&&it->second<u) return false;
    
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
bool checkcore_compress_IIIII(double d,double sim,unsigned int u,size_t i,const intintvec &compress,const intintvec & Hyperedge)
{
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

void compress_IIII_construct(double d,double d_segment_count,std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,compress_index& a,double sim,unsigned int u)
{
    
    clock_t start;
    start = clock();
    /***************得到Nv邻域*******************/  
   
    douintset Neighbor;   //单个neighbor（包含相似度以及编号）
    std::set<size_t> computed_neighbor;//已经计算过了的邻居（不包含相似度）
    //compresstuplevecvec compress(e_id_to_edge.size(),compresstuplevec{});  
    //size_t N = init_nodes.size();//顶点数
    //intintvec tables(N,intvec{});//索引表 有多少个顶点就建多少个表   这里该怎么不用N限制大小？同时又能用下标访问？（为了insert考虑）
    //函数里的tables可以限定大小，类自带的tables不限定大小，函数里面的tables计算完成之后再复制给类里面的tables
    //core-order  tables  cluster-index compress  neighbor 这几个都可以这样处理  在函数里面先算，最后再写到类里面的对应变量去存储
    //std::vector<intintvec> compress(e_id_to_edge.size(),intintvec{});//存储压缩之后的索引结构
    unsigned int geshu1=(unsigned int)(1.0/d);
    std::vector<std::vector<std::vector<size_t>>> compress(e_id_to_edge.size(),std::vector<std::vector<size_t>>(geshu1));
    //double  d=0.1; //划分间隙  固定区间划分还是根据第一个的相似度  固定长度划分?

    //double d_segment_count=0.1;  //分段统计间隙
    unsigned int geshu=(unsigned int)(1.0/d_segment_count);//会不会出问题？
    //std::vector<std::map<size_t,int>> Cluster_Index(geshu,std::map<size_t,int>{});//用map嘛这里
    std::vector<std::vector<core_pair2>> Cluster_Index(geshu,std::vector<core_pair2>{});

    //std::vector<int> segment_count;//初始化为空
    std::vector<int> segment_count;
    core_pair2 cp;
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
       
        
        
        for(unsigned int k=0;k<geshu;k++) segment_count.push_back(0); 
        for(unsigned int k=0;k<geshu;k++)
        {
            segment_count[k]=compress[i][k].size();
        }
        for(unsigned int k=1;k<geshu;k++) segment_count[k]+=segment_count[k-1];// 总和而不是区间内部个数
        for(unsigned int k=0;k<geshu;k++) {
            //大的存在后面 小的存在前面
            if(segment_count[k]!=0) {
                cp.count=-segment_count[k];
                cp.eid=i;
                Cluster_Index[geshu-1-k].push_back(cp);
            }
        }
        computed_neighbor.clear();
        segment_count.clear();
        Neighbor.clear();

    }

    for(unsigned int k=0;k<geshu;k++) {
        sort(Cluster_Index[k].begin(),Cluster_Index[k].end(),comp2);
    }

    //使用set时，对pair第二个值进行排序 就可以对第一个值使用find了
    // auto sortBySecond=[](const std::pair<size_t,int> &a,const std::pair<size_t,int> &b) {return a.second>b.second;};
    // for(unsigned int k=0;k<10;k++)
    // {
    //     std::sort(Cluster_Index[k].begin(),Cluster_Index[k].end(),sortBySecond);
    // }


    std::cout<<"compute ok!"<<'\n';
    

    clock_t e_tm = clock();
    a.init_time=double(e_tm - start) / double(CLOCKS_PER_SEC);
    a.output["inition time"]= std::to_string(a.init_time);
    a.write_indexII(compress);
    a.write_indexIIII(Cluster_Index);//这个时间需要不需要统计？
    // 统计一下内存大小
    size_t  neicun=0;
    for(size_t i=0;i<a.tables.size();i++) {neicun+=a.tables[i].size()*8;}
    for(size_t i=0;i<compress.size();i++) {
        for(auto x:compress[i])
        {
            neicun+=x.size()*8;
        }
    }
    for(size_t i=0;i<Cluster_Index.size();i++) {neicun+=Cluster_Index[i].size()*12;}
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
    while (temp_sim>epision)// 
    {
        i+=1;
        temp_sim-=d_segment_count;
        if(i>=9) break;
    }

    // for(auto x:a.Cluster_Index[i-1])
    // {
    //     cluster_i_jia_1[x.second]=-x.first;
    // }
    unsigned geshu=(unsigned int)(1.0/d_segment_count);
    double temp2=0;
    for(auto x:a.Cluster_Index[i]){
        //if(compress[x.second].empty()) continue;
        if(-x.count<u) break;
        if(edge_to_cluster[x.eid]!=0) continue;
        // if(temp_sim>epision){//
        //     if(!checkcore_compress_IIIII(d,sim,u,x.second,a.compress[x.second],e_id_to_edge)) continue;   
        // }
        edge_to_cluster[x.eid]=Cluster_ID;
        Cluster[Cluster_ID-1].push_back(x.eid); //??????
        std::queue<size_t> Q;
        Q.push(x.eid);
        while(!Q.empty()){
            size_t temp=Q.front();
            Q.pop();
            //计算N_SIM_V
            //compute_compress_II_n_sim(d,sim,temp,a.compress[temp],n_sim,e_id_to_edge);//瓶颈在于这里core节点的邻居扩散？
            for(unsigned k=0;k<geshu;k++)
            {
                temp2=(9-k)*d;
                if(temp2+d-sim<epision) break;
                else{ 
                    for ( auto elem1:a.compress[temp][k]) {
                        //double xiangsidu=geometricMean(sim,e_id_to_edge,,,,);
                        if(edge_to_cluster[elem1]!=Cluster_ID){//还没被加入过当前cluster
                            
                            edge_to_cluster[elem1]=Cluster_ID;
                            Cluster[Cluster_ID-1].push_back(elem1);  //??????    直接在这里定义一个局部变量  最后调用函数将结果保存到csv文件中
                            // if(checkcore_compress_IIIII(d,i,cluster_i_jia_1,d_segment_count,sim,u,elem1,a.compress[elem1],e_id_to_edge)){ //?
                            //     Q.push(elem1);
                            // }
                            unsigned int count_sim=0;
                            for(unsigned int w=0;w<geshu;w++)
                            {
                                temp2=(9-w)*d;
                                if((temp2+d-sim<epision)||count_sim>=u) break;
                                count_sim+=a.compress[elem1][w].size();
                            }
                            if(count_sim>=u){ //?
                                Q.push(elem1);
                            }
                        }
            
                    }
                }
            }
           // n_sim.clear();
        }
        
       // N_SIM.clear();
        Cluster_ID++;
       
    }


    // if(temp_sim<0)//
    // {
    //     for(auto x:a.Cluster_Index[i-1]){
    //         //if(compress[x.second].empty()) continue;
    //         if(-x.first<u) break;
    //         if(edge_to_cluster[x.second]!=0) continue;
    //         if(!checkcore_compress_IIIII(d,sim,u,x.second,a.compress[x.second],e_id_to_edge)) continue;   

    //         edge_to_cluster[x.second]=Cluster_ID;
    //         Cluster[Cluster_ID-1].push_back(x.second); //??????
    //         std::queue<size_t> Q;
    //         Q.push(x.second);
    //         while(!Q.empty()){
    //             size_t temp=Q.front();
    //             Q.pop();
    //             //计算N_SIM_V
    //             compute_compress_II_n_sim(d,sim,temp,a.compress[temp],n_sim,e_id_to_edge);//瓶颈？
    //             for ( auto elem1:n_sim) {
    //                 //double xiangsidu=geometricMean(sim,e_id_to_edge,,,,);
    //                 if(edge_to_cluster[elem1]!=Cluster_ID){//还没被加入过当前cluster
                        
    //                     edge_to_cluster[elem1]=Cluster_ID;
    //                     Cluster[Cluster_ID-1].push_back(elem1);  //??????    直接在这里定义一个局部变量  最后调用函数将结果保存到csv文件中
                    
    //                     if(checkcore_compress_IIIII(d,sim,u,elem1,a.compress[elem1],e_id_to_edge)){ //?
    //                         Q.push(elem1);
    //                     }
    //                 }
        
    //             }
    //             n_sim.clear();
    //         }
            
    //         // N_SIM.clear();
    //         Cluster_ID++;
       
    //     }
    // }

    end = clock();

    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.write_others(Cluster_ID-1);
    a.writecluster("../output/",Cluster);
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
    a.writecluster("../output/",Cluster);
    std::cout<<Cluster_ID-1<<'\n';

}