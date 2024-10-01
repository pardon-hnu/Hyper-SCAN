// #include <bits/stdc++.h>
#include <iostream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <set>
#include <tuple>
#include "hypergraph.h"
#include "readhg.h"
#include "algorithms.h"
#include "gs_index.h"
#include "compress_index.h"

// typedef std::map<std::string, std::string> strstrMap;
typedef std::map<std::string, std::string> strstrMap;
typedef std::tuple<size_t,size_t,size_t> inttriplet;
typedef std::vector<inttriplet> vinttriplet;

double epision=0.2;// 设置的相似度阈值
unsigned int miu=2;// 设置的core阈值   


int main(int argc, char *argv[]) 
{
    
    if (argc >= 2)
    {
        std::istringstream iss( argv[1] );//使用argv[1]来初始化iss对象
        int num_threads;

        if (iss >> num_threads)
        {
            std::cout << num_threads<<"\n";
            Hypergraph h;
            if (argc>=3){
                getHg(argv[2],h);// 无法从.hyp文件中读取数据,则是文件路径没给对
                h.dataset = argv[2];
                std::cout<<"hypergraph ready!"<<'\n';
            }
            // std::string init_type = "nbr"; // or "lub" (local upper bound)
            
            h.initialise();  //初始化 给顶点编号

            std::string alg;

            int iterations;   
            bool log = false;

            if (argc>=4){
                alg = argv[3];
            }
            else{
               alg = "base_line";
            }

            if(argc>=5){
                iterations = atoi(argv[4]); //计算轮次 不给值就默认1次
            }
            else
            iterations = 1;

            if(argc>=6){
                std::string s = argv[5];
                if(s[0]=='1') 
                log = true;           //控制是否给出计算过程中的日志信息
            }

            if(argc>=7){
                miu = atoi(argv[6]);
            }

            if(argc>=8){
                epision = atof(argv[7]);
            }

            std::cout << argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<" "<<argv[6]<<" "<<argv[7]<<"\n";
            std::string filename=argv[2];
            for(int i=1;i<=iterations;i++){
                if(alg == "pSCAN-adp+LI")
                {
                    std::cout <<"pSCAN-adp+LI \n";
                    Algorithm a(h); 
                    //std::cout<<h.hyperedges.size()<<' '<<h.init_nodes.size()<<'\n';
                    optimize_index(h.dataset, h.hyperedges, h.init_nodes, h.node_index,a,epision,miu-1);
                    std::cout<<"Execution time= "<< a.exec_time<<": init_tm= "<<a.output["init_time"]<<"\n";
                }

                if(alg == "GS*-index+LI")
                {
                    std::cout <<"GS*-index+LI \n";
                    Gs_index a(h);  
                    //std::cout<<h.hyperedges.size()<<' '<<h.init_nodes.size()<<'\n';
                    gs_index_construct(h.dataset, h.hyperedges, h.init_nodes, h.node_index,a,epision,miu-1);
                    gs_index_cluster(h.hyperedges,a,epision,miu-1);
                    std::cout<<"Execution time= "<< a.exec_time<<": init_tm= "<<a.init_time<<":  insert time="<<a.insert_time<<"\n";
                }

                if(alg == "LSBI")
                {
                    size_t num_node=0;
                    size_t num_hyperegde=0;
                    std::srand(static_cast<size_t>(std::time(0)));
                    double d1[10]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
                    std::cout <<"LSBI(CI+SI+LI) \n";
                    compress_index a(h);
                    //std::cout<<h.hyperedges.size()<<' '<<h.init_nodes.size()<<'\n';
                    compress_IIII_construct(0.1,0.1,h.dataset, h.hyperedges, h.init_nodes, h.node_index,a,epision,miu-1);//变化tau
                    if(epision==d1[0]||epision==d1[1]||epision==d1[2]||epision==d1[3]||epision==d1[4]||epision==d1[5]||epision==d1[6]||epision==d1[7]||epision==d1[8]||epision==d1[9]) 
                    {        
                        compress_IIII_cluster_bianjie(0.1,0.1,h.dataset, h.hyperedges, h.init_nodes, h.node_index,a,epision,miu-1);
                        std::cout<<"Execution time= "<< a.exec_time<<": init_tm= "<<a.init_time<<"\n";
                         num_node=std::rand() % (h.init_nodes.size())+h.init_nodes.size()/2;
                         if(num_node<h.init_nodes.size())  num_node=h.init_nodes[num_node];
                         num_hyperegde=std::rand() % (2*h.hyperedges.size());
                         updata_insert(a,num_node,num_hyperegde);

                         num_hyperegde=std::rand() % (2*h.hyperedges.size());
                         updata_remove( a,0,num_hyperegde/4);

                    }
                    else 
                    {
                        compress_IIII_cluster(0.1,0.1,h.dataset, h.hyperedges, h.init_nodes, h.node_index,a,epision,miu-1);
                        std::cout<<"Execution time= "<< a.exec_time<<": init_tm= "<<a.init_time<<"\n";
                         num_node=std::rand() % (h.init_nodes.size())+h.init_nodes.size()/2;
                         if(num_node<h.init_nodes.size())  num_node=h.init_nodes[num_node];
                         num_hyperegde=std::rand() % (2*h.hyperedges.size());
                         updata_insert(a,num_node,num_hyperegde);

                         num_hyperegde=std::rand() % (2*h.hyperedges.size());
                         updata_remove( a,0,num_hyperegde/4);
                    }
                }
            }
        }
    }
    
}
