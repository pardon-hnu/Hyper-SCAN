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
//#include "utils.h"

//--------------------------------------------------------------------- Utility functions ------------------------------------------------------------------------------

void Gs_index::writecluster(std::string folder,intintvec Cluster,double sim,unsigned int u){
    //std::cout << "core: \n";把core保存到csv文件中，core是运行过程中的变量，只存在于运行过程，我们要看到他的话要么print，要么写道csv文件里再看
    //std::string file = folder + "cluster_"+output["algo"]+"_"+hg.dataset+".csv";// output是什么？ 好像就是打了个表 用关键词去对应另一个
    
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
        //this->max_neighbor=max_neighbor;
 }

 void Gs_index::write_index(std::vector<std::vector<core_pair>> core_order){

    for(size_t i=0;i<core_order.size();i++)
    {
        (this->core_order).push_back(core_order[i]);
    }

    // for(size_t i=0;i<neighbor_order.size();i++)
    // {
    //     (this->neighbor_order).push_back(neighbor_order[i]);
    // }
 }


void generate_neighbors_by_tables(size_t i,intintvec &e_id_to_edge,intIntMap& node_index,intintvec &tables,std::set<size_t> &neighbors)//这里不是调用类里面的a.tables，而是调用函数里面的tables，因为类里面没有tables，这里用set还是直接用vector？？
{
    //保证每个只记录一次
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

void generate_neighbors_by_enumerate(size_t i,intintvec &e_id_to_edge,std::set<size_t> &neighbors)//传引用，一次记下所有的邻居，再去挨个计算相似度
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



double geometricMean_gs1(std::vector<HashSet> &hash,const intintvec & Hyperedge,unsigned int num1,unsigned int num2,size_t e1,size_t e2) {        /**********得到相似度*/
    
    double result = sqrt(static_cast<double>(num1) * num2);
    double divisionResult = 0;  
    unsigned int count=0;

    size_t temp=0;
    if(num1>num2) {temp=e1;e1=e2;e2=temp;}// e1代表两条超边中规模较小的那条
    
    for (unsigned int j=0;j<Hyperedge[e1].size();j++) {
              //auto j=node_index[x];
             // auto it = hash[e2].find(Hyperedge[e1][j]); //hash[e2][j] 会在没有找到的时候自动创建一个？
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

    // size_t temp=0;
    // if(num1>num2) {temp=e1;e1=e2;e2=temp;}// e1代表两条超边中规模较小的那条
    
    for (unsigned int j=0;j<Hyperedge[e1].size();j++) {
              auto i=node_index[Hyperedge[e1][j]];
             // auto it = hash1[e2].find(Hyperedge[e1][j]); //因为这里变慢？
              if (hash1[i]==1) {
                count++;
              } 
    }
    divisionResult = static_cast<double>(count) / result;
    return divisionResult;
}
double geometricMean_gs(const intintvec & Hyperedge,unsigned int num1,unsigned int num2,size_t e1,size_t e2) {        /**********得到相似度*/
    
    double result = sqrt(static_cast<double>(num1) * num2);
    double divisionResult = 0;  
   // unsigned int count=0;
   std::set<size_t> gongtongdingdian;
    for (unsigned int j=0;j<Hyperedge[e1].size();j++) {
              auto it = std::find(Hyperedge[e2].begin(), Hyperedge[e2].end(), Hyperedge[e1][j]);
              if (it != Hyperedge[e2].end()) {
                //if(node_index[Hyperedge[e1][j]]<i) {return -2;}//已经在前面的表中算过了
                //count++;//如果超边里面有重复顶点 这里直接用count统计的话有被统计多次的可能 同时也有只统计一次的时候 会造成差异？
                gongtongdingdian.insert(Hyperedge[e1][j]);
              }    
    }
    divisionResult = static_cast<double>(gongtongdingdian.size()) / result;
    return divisionResult;
}

void gs_index_construct(std::string dataset, intintvec e_id_to_edge,intvec init_nodes, intIntMap& node_index,Gs_index& a,double sim,unsigned int u)
{
    
   
    clock_t start,e_tm;
    start = clock();
    /***************得到Nv邻域*******************/  
    std::vector<core_pair> Neighbor;
    //douintsetvec Neighbor(e_id_to_edge.size(),douintset{});   //neighbor-order
    size_t N = init_nodes.size();//顶点数
    intintvec tables(N,intvec{});//索引表 有多少个顶点就建多少个表
    std::set<size_t> computed_neighbor;//已经计算过了的邻居（不包含相似度）
    std::vector< HashSet > hash;
    hash.reserve(e_id_to_edge.size());
    //std::vector<std::vector<short>> hash1(e_id_to_edge.size(),std::vector<short>(5000));
    for(size_t i=0;i<e_id_to_edge.size();i++)
    {
        hash.emplace_back(e_id_to_edge[i].size());
    }

    for(size_t i=0;i<e_id_to_edge.size();i++)//建表
    {
        for (const auto& elem:e_id_to_edge[i]) {
             auto j = node_index[elem];
             tables[j].push_back(i);
            // hash[i].insert(elem);  //超边中有相同顶点只被添加一次 按理说不该有相同值
            // hash1[i][j]=1;
        }

    }
    e_tm = clock();
    a.LI_time=double(e_tm - start) / double(CLOCKS_PER_SEC);
    std::cout<<a.LI_time<<'\n';

    double xiangsidu=0;
    core_pair cp;
    for(size_t i=0;i<e_id_to_edge.size();i++)
    {

        generate_neighbors_by_tables(i,e_id_to_edge,node_index,tables,computed_neighbor);
        //generate_neighbors_by_enumerate(i,e_id_to_edge,computed_neighbor);
        for(auto v:computed_neighbor)
        {
           // xiangsidu=geometricMean_gs1(hash,e_id_to_edge,e_id_to_edge[i].size(),e_id_to_edge[v].size(),i,v);  //有问题
           xiangsidu=geometricMean_gs(e_id_to_edge,e_id_to_edge[i].size(),e_id_to_edge[v].size(),i,v);
           // if(e_id_to_edge[i].size()<e_id_to_edge[v].size()) xiangsidu=geometricMean_gs2(node_index,hash1[v],e_id_to_edge,e_id_to_edge[i].size(),e_id_to_edge[v].size(),i,v);
          // else {xiangsidu=geometricMean_gs2(node_index,hash1[i],e_id_to_edge,e_id_to_edge[v].size(),e_id_to_edge[i].size(),v,i);}
            
            if(i!=v){
                cp.eid=v;
                cp.xiangsidu=-xiangsidu;
                Neighbor.push_back(cp);  //1-
            }
        }

        sort(Neighbor.begin(),Neighbor.end(),comp);
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

        //         double xiangsidu=geometricMean_gs(e_id_to_edge,e_id_to_edge[i].size(),e_id_to_edge[tables[j][k]].size(),i,tables[j][k]);
        //         if(i!=tables[j][k]){
        //             Neighbor[i].insert(std::make_pair(-xiangsidu,tables[j][k]));//-1
        //         }
        //     }

        // }
        (a.neighbor_order).push_back(Neighbor);
        Neighbor.clear();
        computed_neighbor.clear();
    }

    std::cout<<"compute ok!"<<'\n';
    //construct core-order
    size_t max_u=0;
    for(size_t i=0;i<(a.neighbor_order).size();i++)
    {
        if(a.neighbor_order[i].size()>max_u) max_u=a.neighbor_order[i].size();
    }
    //douintsetvec core_order(max_u,douintset{});// 虽然能够自动排序，自动判断是否有重复，但是无法按照下标进行访问
    std::vector<std::vector<core_pair>> core_order(max_u,std::vector<core_pair>{});
    //如果要大量访问的话先将其中内容复制到vector中？  否则就用迭代器？  范围访问？
    // douintset mySet;
    // intdouprvec myVector(mySet.begin(), mySet.end());
    
    // 使用范围遍历打印每个元素   在聚类的时候对于set可以这样访问每一个元素
    //   for(auto x:Neighbor[j])
    //         {
    //             std::cout<<x.first<<' '<<x.second<<'\n';
    //         }
    
    //   for(auto x:core_order[j])
    //         {
    //             std::cout<<x.first<<' '<<x.second<<'\n';
    //         }
    std::cout<<"max_u ok!"<<'\n';

    // for(size_t i=1;i<=max_u;i++)
    // {
    //     //u=i时的gs-index存在core-order[i-1]
    //     //neighbor[j]中的第u个节点存在neighbor[j][i-1]处

    //     for(size_t j=0;j<Neighbor.size();j++)
    //     {
    //         if(Neighbor[j].size()>=i) {
    //              auto it = Neighbor[j].begin();    //这里如果不能忍受  core-order就别设置为与set有关的类型
    //              std::advance(it, i-1);
    //              core_order[i-1].insert(std::make_pair(it->first,j));   //no 1-?
    //         }

    //     }
    // }
    for(size_t j=0;j<(a.neighbor_order).size();j++)
    {
        unsigned int i=1;
        for(auto x: a.neighbor_order[j])
        {
          cp.eid=j;
          cp.xiangsidu=x.xiangsidu;
          core_order[i-1].push_back(cp);   //no 1-?
          ++i;
        }
    }
    for(size_t j=0;j<max_u;j++)
    {
        sort(core_order[j].begin(),core_order[j].end(),comp);
    }
    std::cout<<"core-order ok!"<<'\n';
    a.write_index(core_order);//这个时间需要？
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
    std::cout<<"neicun:  "<<neicun<<'\n';
    
}

void gs_index_cluster(intintvec e_id_to_edge,Gs_index& a,double sim,unsigned int u)
{
    a.output["algo"] = "optimize_gs_index";
    // size_t max_u=0;
    // for(size_t i=0;i<a.neighbor_order.size();i++)
    // {
    //     if(a.neighbor_order[i].size()>max_u) max_u=a.neighbor_order[i].size();
    // }

    clock_t start, end;
    start = clock();
    /**********处理过程*********/
    intvec edge_to_cluster(e_id_to_edge.size(),0);//判断超边是否已经被加入过当前cluster
    intintvec Cluster(e_id_to_edge.size(),intvec{});//保存cluster 
    //  还是要把cluster里面超边的连接信息保存下来，不能只保存一个超边编号，会导致只知道他们在一个cluster里面而不知道他们的连接关系(实际上是可以知道的 不过要凭超边编号去原图中找连接关系)
    // 和普通图还是区别很大的 普通图中没有一条边和有一条边是不同的子图   超图中超边被加入cluster之后，由于共享节点的关系 不止会和引入它的超边有连接关系
    size_t Cluster_ID=1;

    for(const auto &elem:a.core_order[u-1]){
        if(-elem.xiangsidu<sim) break;    //1-
        if(edge_to_cluster[elem.eid]!=0) continue;

        edge_to_cluster[elem.eid]=Cluster_ID;
        Cluster[Cluster_ID-1].push_back(elem.eid); //??????
        std::queue<size_t> Q;
        Q.push(elem.eid);
        while(Q.empty()!=true){
            size_t temp=Q.front();
            Q.pop();
            for ( auto elem1:a.neighbor_order[temp]) {
                if(-elem1.xiangsidu<sim) break;   //1-

                if(edge_to_cluster[elem1.eid]!=Cluster_ID){//还没被加入过当前cluster
                            
                    edge_to_cluster[elem1.eid]=Cluster_ID;
                    Cluster[Cluster_ID-1].push_back(elem1.eid);  //??????    直接在这里定义一个局部变量  最后调用函数将结果保存到csv文件中
                    if(a.neighbor_order[elem1.eid].size()>=u){ //?
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
    a.writecluster("../output/",Cluster,sim,u);
    std::cout<<Cluster_ID-1<<'\n';

}

void gs_index_insert( Gs_index &a,intvec insert_hyperedge)
{
    //插入的超边用一个intvec数组存储
    //insert一条超边需要改(a.hg).hyperedges  .init_nodes  .node_index
    //还需要改a.neighbor_order  core_order  max_neighbor(可能)
    clock_t start,end;
    start=clock();
    (a.hg.hyperedges).push_back(insert_hyperedge);

    size_t size_node=(a.hg.init_nodes).size();
    for(size_t i=0;i<insert_hyperedge.size();i++)//gs_index里面在插入删除的时候好像用不上 这几个值，可以不改？
    {
        auto result = (a.hg.node_index).emplace(insert_hyperedge[i],size_node);
        if(result.second) {
            size_node+=1;
            (a.hg.init_nodes).push_back(insert_hyperedge[i]);
        }

    }
    size_t hyperedge_size=(a.hg.hyperedges).size()-1;
    std::set<size_t> neighbors;//存储邻接超边
    std::vector<core_pair> computed_neighbors;
    core_pair cp;
    size_t piot=0;//插入位置
    generate_neighbors_by_enumerate(hyperedge_size,a.hg.hyperedges,neighbors);

    if(neighbors.size()>a.max_neighbor) a.max_neighbor=neighbors.size();//

    for(auto v:neighbors)
    {
        double xiangsidu=geometricMean_gs(a.hg.hyperedges,insert_hyperedge.size(),(a.hg.hyperedges)[v].size(),hyperedge_size,v);
        cp.xiangsidu=-xiangsidu;//这里由于之前的设置 还是使用一个负的相似度，然后从小到大排序
        if(hyperedge_size!=v){
            cp.eid=v;
            for(size_t j=0;j<computed_neighbors.size();j++)//可以采用二分查找
            {
                if(computed_neighbors[j].xiangsidu>=cp.xiangsidu) {piot=j;break;}
            }
            computed_neighbors.insert(computed_neighbors.begin()+piot,cp);  //1-


            cp.eid=hyperedge_size;
            for(size_t j=0;j<(a.neighbor_order[v]).size();j++)//可以采用二分查找
            {
                if(a.neighbor_order[v][j].xiangsidu>=cp.xiangsidu) {piot=j;break;}
            }
            a.neighbor_order[v].insert((a.neighbor_order[v]).begin()+piot,cp);//还要把这条超边加入到neighbors每条超边中
            

            //修改core-order   将neighbors超边 信息在core-order删除(在hyperedge size之后的)
            // size_t count=0;
            // std::pair<double,size_t> temp=std::make_pair(-xiangsidu,hyperedge_size);  cp
            // for(auto x:a.neighbor_order[v])//得到hyperedge_size的插入位置 下标(从0开始) piot
            // {
            //     if(x!=cp) count+=1;
            //     else break;
            // }
            size_t piot2=0;
            auto it=(a.neighbor_order[v]).begin()+piot;
            ++it;
            size_t count2=piot;
            for (;it!= a.neighbor_order[v].end(); ++it){//删
                for(size_t j=0;j<(a.core_order[count2]).size();j++)//可以采用二分查找
                {
                    if(a.core_order[count2][j].xiangsidu>=cp.xiangsidu) {piot2=j;break;}
                }
                a.core_order[count2].erase((a.core_order[count2]).begin()+piot2);
                count2+=1;
            }

            piot2=0;
            it=(a.neighbor_order[v]).begin()+piot;
            count2=piot;
            for (;it!= a.neighbor_order[v].end(); ++it){//插
                for(size_t j=0;j<(a.core_order[count2]).size();j++)//可以采用二分查找
                {
                    if(a.core_order[count2][j].xiangsidu>=cp.xiangsidu) {piot2=j;break;}
                }
                a.core_order[count2].insert((a.core_order[count2]).begin()+piot2,cp);
                count2+=1;
            }

        }
    }
    a.neighbor_order.push_back(computed_neighbors);//修改neighbor_order  
    end=clock();
    a.insert_time=double(end - start) / double(CLOCKS_PER_SEC);
    // size_t max_u=0;
    // for(size_t i=0;i<a.neighbor_order.size();i++)
    // {
    //     if(a.neighbor_order[i].size()>max_u) max_u=a.neighbor_order[i].size();
    // }
    // a.max_neighbor=max_u;
    std::cout<<"insert time  "<<a.insert_time<<'\n';
}

void gs_index_insert_vertex( Gs_index &a,intvec contain_hyperedges)
{


}

// void gs_index_insert( Gs_index &a,intvec insert_hyperedge)
// {
//     //插入的超边用一个intvec数组存储
//     //insert一条超边需要改(a.hg).hyperedges  .init_nodes  .node_index
//     //还需要改a.neighbor_order  core_order  max_neighbor(可能)
//     clock_t start,end;
//     start=clock();
//     (a.hg.hyperedges).push_back(insert_hyperedge);

//     size_t size_node=(a.hg.init_nodes).size();
//     for(size_t i=0;i<insert_hyperedge.size();i++)//gs_index里面在插入删除的时候好像用不上 这几个值，可以不改？
//     {
//         auto result = (a.hg.node_index).emplace(insert_hyperedge[i],size_node);
//         if(result.second) {
//             size_node+=1;
//             (a.hg.init_nodes).push_back(insert_hyperedge[i]);
//         }

//     }
//     size_t hyperedge_size=(a.hg.hyperedges).size()-1;
//     std::set<size_t> neighbors;//存储邻接超边
//     douintset computed_neighbors;
//     generate_neighbors_by_enumerate(hyperedge_size,a.hg.hyperedges,neighbors);

//     //if(neighbors.size()>a.max_neighbor) a.max_neighbor=neighbors.size();//对于core order修改有用,在修改core-order时改这个

//     for(auto v:neighbors)
//     {
//         double xiangsidu=geometricMean_gs(a.hg.hyperedges,insert_hyperedge.size(),(a.hg.hyperedges)[v].size(),hyperedge_size,v);
//         if(hyperedge_size!=v){
//             auto result=a.neighbor_order[v].insert(std::make_pair(-xiangsidu,hyperedge_size));//还要把这条超边加入到neighbors每条超边中
//             computed_neighbors.insert(std::make_pair(-xiangsidu,v));  //1-
//             //修改core-order   将neighbors超边 信息在core-order删除(在hyperedge size之后的)
//             size_t count=0;
//             std::pair<double,size_t> temp=std::make_pair(-xiangsidu,hyperedge_size);
//             for(auto x:a.neighbor_order[v])//得到hyperedge_size的插入位置 下标(从0开始)
//             {
//                 if(x!=temp) count+=1;
//                 else break;
//             }

//             auto it=result.first;
//             ++it;
//             size_t count2=count;
//             for (;it!= a.neighbor_order[v].end(); ++it){//删
//                 a.core_order[count2].erase(std::make_pair(it->first,v));
//                 count2+=1;
//             }

//             it=result.first;
//             count2=count;
//             for (;it!= a.neighbor_order[v].end(); ++it){//插
//                 a.core_order[count2].insert(std::make_pair(it->first,v));
//                 count2+=1;
//             }

//         }
//     }
//     a.neighbor_order.push_back(computed_neighbors);//修改neighbor_order  
//     end=clock();
//     a.insert_time=double(end - start) / double(CLOCKS_PER_SEC);
//     std::cout<<"insert time  "<<a.insert_time<<'\n';
// }

// void gs_index_remove(Gs_index &a,size_t eid)
// {
//     if(eid<0||eid>=(a.hg.hyperedges).size()) return;
//     clock_t start,end;
//     start=clock();

//    // std::vector<size_t> neighbors;//存储邻接超边
//     //douintset computed_neighbors;
//    // generate_neighbors_by_enumerate(eid,a.hg.hyperedges,neighbors);
//     // for(auto x:(a.neighbor_order)[eid])
//     // {
//     //     neighbors.push_back(x.second);
//     // }

//     (a.hg.hyperedges).erase((a.hg.hyperedges).begin()+eid);//在超图中删除eid

//     for(auto v:(a.neighbor_order)[eid])
//     {
//             size_t count=0;
//             std::pair<double,size_t> temp=std::make_pair(v.first,eid);
//             auto result=a.neighbor_order[v.second].begin();

//             for(auto x:a.neighbor_order[v.second])//得到hyperedge_size的插入位置 下标(从0开始)
//             {
//                 if(x!=temp) {count+=1;++result;}
//                 else break;
//             }
            
//             auto it=result;

//             size_t count2=count;
//             for (;it!= a.neighbor_order[v.second].end(); ++it){//删
//                 a.core_order[count2].erase(std::make_pair(it->first,v.second));
//                 count2+=1;
//             }
            
//             it=result;
//             ++it;
//             count2=count;
//             for (;it!= a.neighbor_order[v.second].end(); ++it){//插
//                 a.core_order[count2].insert(std::make_pair(it->first,v.second));
//                 count2+=1;
//             }
            
//             a.neighbor_order[v.second].erase(temp);

//     }

//     (a.neighbor_order).erase((a.neighbor_order).begin()+eid);//删除neighbor-order[i]

//     end=clock();
//     a.remove_time=double(end - start) / double(CLOCKS_PER_SEC);
//     std::cout<<"remove time  "<<a.remove_time<<'\n';
// }
