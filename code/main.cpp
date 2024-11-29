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

double epision=0.2;
unsigned int miu=2;

int main(int argc, char *argv[]) 
{
    
    if (argc==5)
    {
        // std::istringstream iss( argv[1] );//
        // int num_threads=1;

        // if (iss >> num_threads)
        // {
        //     std::cout << num_threads<<"\n";
            Hypergraph h;
            // if (argc>=3){
                getHg(argv[1],h);//
                h.dataset = argv[1];
                std::cout<<"hypergraph ready!"<<'\n';
           // }
            // std::string init_type = "nbr"; // or "lub" (local upper bound)
            
            h.initialise();  

            std::string alg;

            int iterations=1;   
            // bool log = false;

            //if (argc>=4){
                alg = argv[2];
           // }
            // else{
            //    alg = "optimize_index";
            // }

            // if(argc>=5){
            //     iterations = atoi(argv[4]); 
            // }
            // else
            // iterations = 1;

            // if(argc>=6){
            //     std::string s = argv[5];
            //     if(s[0]=='1') 
            //     log = true;         
            // }

           // if(argc>=5){
                miu = atoi(argv[3]);
           // }

           // if(argc>=6){
                epision = atof(argv[4]);
           // }

            std::cout << argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
            std::string filename=argv[1];
            for(int i=1;i<=iterations;i++){
                if(alg == "pHSCAN")
                {
                    std::cout <<"pHSCAN \n";
                    Algorithm a(h); 
                    optimize_index(h.dataset, h.hyperedges, h.init_nodes, h.node_index,a,epision,miu-1);
                    std::cout<<"Execution time= "<< a.exec_time<<": init_tm= "<<a.output["init_time"]<<"\n";
                }
                if(alg == "OI-construction")
                {
                    std::cout <<"OI-construction \n";
                    Gs_index a(h);  
                    gs_index_construct(h.dataset, h.hyperedges, h.init_nodes, h.node_index,a);
                }
                if(alg == "LSBI-construction")
                {
                    std::cout <<"LSBI-construction \n";
                    compress_index a(h);
                    compress_IIII_construct(0.1,0.1,h.dataset, h.hyperedges, h.init_nodes, h.node_index,a); 
                }
                // if(alg == "SQuery-OI")
                // {
                //     std::cout <<"optimize_gs_index \n";
                //     Gs_index a(h);  
                //     gs_index_construct(h.dataset, h.hyperedges, h.init_nodes, h.node_index,a);
                //     gs_index_cluster(h.hyperedges,a,epision,miu-1);
                //     std::cout<<"Execution time= "<< a.exec_time<<": init_tm= "<<a.init_time<<":  insert time="<<a.insert_time<<"\n";
                // }
                if(alg == "SQuery")
                {
                    std::cout <<"SQuery \n";
                    compress_index a(h);
                    //std::cout<<h.hyperedges.size()<<' '<<h.init_nodes.size()<<'\n';
                    compress_IIII_construct(0.1,0.1,h.dataset, h.hyperedges, h.init_nodes, h.node_index,a); 
                    compress_IIII_cluster_bianjie(0.1,0.1,h.dataset, h.hyperedges, h.init_nodes, h.node_index,a,epision,miu-1);
                    std::cout<<"Execution time= "<< a.exec_time<<": init_tm= "<<a.init_time<<"\n";
                }

                if(alg == "PQuery")
                {
                    std::cout <<"PQuery \n";
                   // unsigned int workers_number=32;
                    compress_index a(h); 
                    //std::cout<<h.hyperedges.size()<<' '<<h.init_nodes.size()<<'\n';
                    compress_IIII_construct(0.1,0.1,h.dataset, h.hyperedges, h.init_nodes, h.node_index,a);
                    compress_IIII_cluster_parallel(0.1,0.1,h.dataset, h.hyperedges, h.init_nodes, h.node_index,a,epision,miu-1);
                    std::cout<<"Execution time= "<< a.exec_time<<": init_tm= "<<a.init_time<<"\n";
                }
            }
       // }
    }
    
}
