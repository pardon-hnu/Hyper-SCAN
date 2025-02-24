#pragma once


#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <string.h>
#include <math.h>
#include <queue>
#include <algorithm>
#include "global.h"

using namespace std;

// typedef struct
// {
//     int id;
//     double similarity;
// }item_t; 
typedef pair<int,double> item_t; 
typedef struct 
{
    int id;
    int weight;
}weight_t;

static bool cmp(const item_t &a, const item_t &b) {
    if(b.second == a.second)
        return a.first < b.first;
    return a.second > b.second;
}

class WSCAN
{
    private:
        vector<vector<weight_t>> Graph;
        vector<vector<item_t>> CO;
        vector<vector<item_t>> NO;
        vector<int> res;
        string dataset;
    public:
        WSCAN(string dataset)
        {
            this->dataset=dataset;
            string graph_path=dataset_to_filename[dataset];
            load_graph(graph_path);
            this->res.resize(Graph.size());
        }
        void load_graph(string graph_path)
        {
            printf("[LOAD GRAPH] start!!!\n");
            ifstream infile;
            infile.open(graph_path);
            string line;
            getline(infile, line);
            char* temp;
            temp = strtok(const_cast<char*>(line.c_str()), " ");
            int graph_size = atoi(temp);
            printf("graph size: %d\n",graph_size);
            temp = strtok(NULL, " ");
            int edge_num = atoi(temp);
            printf("edge size: %d\n",edge_num);
            Graph.resize(graph_size);
            while (getline(infile, line))
            {
                if (line.size() == 0)
                {
                    break;
                }
                temp = strtok(const_cast<char*>(line.c_str()), " ");
                int u = atoi(temp);
                temp = strtok(NULL, " ");
                int v= atoi(temp);
                temp = strtok(NULL, " ");
                int w= atoi(temp);
                weight_t e1,e2;
                e1.id=v;
                e1.weight=w;
                e2.id=u;
                e2.weight=w;
                Graph[u].push_back(e1);
                Graph[v].push_back(e2);
            }
            infile.close();
            printf("[LOAD GRAPH] end!!!\n");
        }
        double compute_weight_similarity(int v1,int v2)
        {
            double sum1=0; //common
            double sum2=0; //v1
            double sum3=0; //v2

            weight_t temp1,temp2;
            temp1.id=v1;temp1.weight=1;
            temp2.id=v2;temp2.weight=2;
            vector<weight_t> s1=Graph[v1];
            s1.push_back(temp1);
            vector<weight_t> s2=Graph[v2];
            s2.push_back(temp2);

            for(auto nei:Graph[v1])
            {
                sum2+=nei.weight*nei.weight;
                for(auto v:Graph[v2])
                {
                    if(v.id==nei.id)
                    {
                        sum1+=v.weight*nei.weight;
                    }
                }
            }
            for(auto nei:Graph[v2])
            {
                sum3+=nei.weight*nei.weight;
            }
            double score=sum1/(sqrt(sum2)*sqrt(sum3));
            return score;
        }
        void construct_index()
        {
            printf("[CONSTRUCT NO] begin!!!!\n");
            NO.resize(Graph.size());
            for(int i=0;i<Graph.size();i++)
            {
                for(auto nei:Graph[i])
                {
                    if(nei.id>i)
                    {
                        double sim=compute_weight_similarity(i,nei.id);
                        NO[i].push_back(make_pair(nei.id,sim));
                        NO[nei.id].push_back(make_pair(i,sim));
                    }
                }
                sort(NO[i].begin(),NO[i].end(),cmp);
            }
            printf("[CONSTRUCT NO] end!!!!\n");

            printf("[CONSTRUCT CO] begin!!!!\n");
            int Max_deg=0;
            for(int i=0;i<Graph.size();i++)
            {
                if(Graph[i].size()>Max_deg)
                {
                    Max_deg=Graph[i].size();
                }
            }
            CO.resize(Max_deg+1);

            for(int i=0;i<Graph.size();i++)
            {
                for(int j=0;j<NO[i].size();j++)
                {
                    CO[j+1].push_back(make_pair(i,NO[i][j].second));
                }
            }
            for(int i=1;i<Max_deg+1;i++)
            {
                sort(CO[i].begin(),CO[i].end(),cmp);
            }
            printf("[CONSTRUCT CO] end!!!!\n");
        }
        void test_print_no()
        {
            for(int i=0;i<Graph.size();i++)
            {
                printf("id: %d and its neighbor: ",i);
                for(auto nei:NO[i])
                {
                    printf("(%d:%lf),",nei.first,nei.second);
                }
                printf("\n");
            }
        }
        void test_print_co()
        {
            for(int i=1;i<CO.size();i++)
            {
                printf("CO[%d]: ",i);
                for(auto item:CO[i])
                {
                    printf("(%d:%lf),",item.first,item.second);
                }
                printf("\n");
            }
        }
        void batch_pSCAN(vector<int> mus, vector<double> epsilons)
        {
            for(auto mu:mus)
            {
                for(auto epsilon:epsilons)
                {
                    printf("[CLUSTER] mu:%d epsilon %lf \n",mu,epsilon);
                    pSCAN(mu-1,epsilon);
                }
            }
        }
        void pSCAN(int mu,double epsilon)
        {
            
            printf("[INIT] init cluster id of each vertex \n");
            for(int i=0;i<Graph.size();i++)
            {
                res[i]=-1;
            }

            printf("[PROCESS] cluster start !!!! \n");
            int ID=0;
            for(auto item:CO[mu])
            {
                if(res[item.first]!=-1)
                {continue;}
                if(item.second<epsilon)
                {break;}
                // cout<<"now obtain cluster "<<ID<<endl;
                queue<int> BFS;
                BFS.push(item.first);
                // printf("[CLUSTER %d]: ",ID);
                while(BFS.size()!=0)
                {
                    int v=BFS.front();
                    BFS.pop();
                    if(res[v]!=-1) continue;
                    // cout<<"expand core vertex: "<<v<<" , and its neighbor: ";
                    res[v]=ID;
                    // printf("%d,",v);
                    for(auto nei:NO[v])
                    {
                        if(nei.second<epsilon)
                        {
                            break;
                        }
                        int u=nei.first;
                        if(res[u]==ID) continue;
                        // printf("%d,",u);
                        if(NO[u][mu].second>=epsilon)
                        {
                            BFS.push(u);
                            // cout<<u<<"(core),";
                        }
                        else
                        {
                            if(res[u]==-1)
                            {
                                res[u]=ID;
                            }
                            // cout<<u<<"(non-core),";
                        }
                    }
                    
                }
                // printf("\n");
                ID++;
            }
            printf("[OUTPUT] print the result !!!! \n");
            print_res(mu+1,epsilon);
        }
        void print_res(int mu,double epsilon)
        {
            string file="./res/"+dataset+"-"+std::to_string(epsilon)+"-"+std::to_string(mu)+".txt";
            cout<<"write to "<<file<<endl;
            ofstream fwrite;
            fwrite.open(file);
            fwrite<<"c/n vertex_id cluster_id"<<endl;
            for(int i=0;i<res.size();i++)
            {
                if(res[i]==-1) continue;
                if(NO[i][mu].second>=epsilon)
                {
                    fwrite<<"c "<<i<<" "<<res[i]<<endl;
                }
                else
                {
                    fwrite<<"n "<<i<<" "<<res[i]<<endl;
                }
            }
        }
};