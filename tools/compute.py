import math
import os
import sys
sys.path.append("..")
import argparse
import string
import re
import readline


def read_graph(dataset_file):
    hyperedges=[]
    gs=0
    infile=open(dataset_file,'r')
    for line in infile:
        temp_string=line.strip('\n')
        temp=re.split(',',temp_string)
        e=[]
        for str_e in temp:
            if str_e!='':
                e.append(int(str_e))
                if gs<int(str_e):
                    gs=int(str_e)
        hyperedges.append(e)
    infile.close()
    print('load graph successfully, vertex num: '+str(gs+1)+' hyperedge num: '+str(len(hyperedges)))
    return hyperedges,gs+1

def construct_link_index(hyperedges,gs):
    vertex=[]
    for i in range(gs):
        vertex.append(set())
    for i in range(len(hyperedges)):
        for v in hyperedges[i]:
            vertex[v].add(i)
    print('construct link index successfully')
    return vertex

def load_hscan_cluster(hscan_result_file,hyperedges):
    hscan_ecluster=[]
    infile=open(hscan_result_file,'r')
    for line in infile:
        temp_string=line.strip('\n')
        temp=re.split(',',temp_string)
        e=[]
        for str_e in temp:
            if str_e!='':
                e.append(int(str_e))
        hscan_ecluster.append(e)
    infile.close()
    print('load the cluster result of HSCAN successfully, the cluster num: '+str(len(hscan_ecluster)))
    hscan_vcluster=[]
    for i in range(len(hscan_ecluster)):
        vs=set()
        for e in hscan_ecluster[i]:
            for v in hyperedges[e]:
                vs.add(v)
        hscan_vcluster.append(list(vs))
    return hscan_ecluster,hscan_vcluster

def load_scan_cluster(scan_result_file,hyperedges):
    scan_vcluster=[]
    cluster_id={}
    infile=open(scan_result_file,'r')
    new_id=0
    for line in infile:
        temp_string=line.strip('\n')
        if temp_string=='c/n vertex_id cluster_id':
            continue
        temp=re.split(' ',temp_string)
        v=int(temp[1])
        if temp[2] in cluster_id.keys():
            # print('cluster_id: '+temp[2]+' and new id : '+str(cluster_id[temp[2]]))
            scan_vcluster[cluster_id[temp[2]]].append(v)
        else:
            scan_vcluster.append([v])
            cluster_id[temp[2]]=new_id
            new_id+=1
    infile.close()
    print('load the cluster result of SCAN successfully, the cluster num: '+str(len(scan_vcluster)))
    scan_ecluster=[]
    for i in range(len(scan_vcluster)):
        vd={}
        for vertex in scan_vcluster[i]:
            vd[vertex]=1
        edge=set()
        for j in range(len(hyperedges)):
            flag=1
            for v in hyperedges[j]:
                if v not in vd.keys():
                    flag=0
                    break
            if flag==1 :
                edge.add(j)
        scan_ecluster.append(list(edge))
    return scan_ecluster,scan_vcluster

def load_wscan_cluster(wscan_result_file,hyperedges):
    wscan_vcluster=[]
    cluster_id={}
    infile=open(wscan_result_file,'r')
    new_id=0
    for line in infile:
        temp_string=line.strip('\n')
        if temp_string=='c/n vertex_id cluster_id':
            continue
        temp=re.split(' ',temp_string)
        v=int(temp[1])
        if temp[2] in cluster_id.keys():
            # print('cluster_id: '+temp[2]+' and new id : '+str(cluster_id[temp[2]]))
            wscan_vcluster[cluster_id[temp[2]]].append(v)
        else:
            wscan_vcluster.append([v])
            cluster_id[temp[2]]=new_id
            new_id+=1
    infile.close()
    print('load the cluster result of Weighted SCAN successfully, the cluster num: '+str(len(wscan_vcluster)))
    wscan_ecluster=[]
    for i in range(len(wscan_vcluster)):
        vd={}
        for vertex in wscan_vcluster[i]:
            vd[vertex]=1
        edge=set()
        for j in range(len(hyperedges)):
            flag=1
            for v in hyperedges[j]:
                if v not in vd.keys():
                    flag=0
                    break
            if flag==1 :
                edge.add(j)
        wscan_ecluster.append(list(edge))
    return wscan_ecluster,wscan_vcluster

def compute_modularity_2(hyperedges,vertex,ecluster,vcluster):
    modularity=0
    connection_sum=0
    degree_sum=0
    
    connection_sum=len(hyperedges)
    for v_item in vertex:
        degree_sum+=len(v_item)
    
    overlap={}
    for c in vcluster:
        for v in c:
            if v in overlap.keys():
                overlap[v]+=1
            else:
                overlap[v]=1

    # print('dataset has '+str(connection_sum)+' hyperedge')
    # print('dataset the sum of d(v) is '+str(degree_sum))

    for i in range(len(vcluster)):
        cluster_connection=0
        cluster_degree=0
        cluster_connection=len(ecluster[i])
        # print('cluster '+str(i)+' has '+str(cluster_connection)+' hyperedges')
        left=cluster_connection

        for v in vcluster[i]:
            cluster_degree+=len(vertex[v])/overlap[v]

        # print('the sum of d(v) in cluster '+str(i)+' is '+str(cluster_degree))
 
        div=cluster_degree/degree_sum

        dic={}
        for e in ecluster[i]:
            if len(hyperedges[e]) in dic.keys():
                dic[len(hyperedges[e])]+=1
            else:
                dic[len(hyperedges[e])]=1

        right=0
        for s in dic.keys():
            right+=pow(div,s)*dic[s]
        modularity_cluster=left-right
        # print('modularity of cluster '+str(i)+' is '+str(modularity_cluster))
        modularity+=left-right

    return modularity/connection_sum


epsilon_hscans={
    '0.1': '0.100000',
    '0.2': '0.200000',
    '0.3': '0.300000',
    '0.4': '0.400000',
    '0.5': '0.500000',
    '0.6': '0.600000',
    '0.7': '0.700000',
    '0.8': '0.800000',
    '0.9': '0.900000'
}
epsilon_wscans={
    '0.1': '0.100000',
    '0.2': '0.200000',
    '0.3': '0.300000',
    '0.4': '0.400000',
    '0.5': '0.500000',
    '0.6': '0.600000',
    '0.7': '0.700000',
    '0.8': '0.800000',
    '0.9': '0.900000'
}

#set argument
parser=argparse.ArgumentParser()
parser.add_argument('--epsilon','-e',nargs=1,required=True,type=str,dest='epsilon')
parser.add_argument('--mu','-m',nargs=1,required=True,type=str,dest='mu')
args=parser.parse_args()

#get args
epsilon=args.epsilon[0]
mu=args.mu[0]
epsilon_scan=epsilon
epsilon_hscan=epsilon_hscans[epsilon]
epsilon_wscan=epsilon_wscans[epsilon]

dataset_file='/home/hnu/Disk0/ParDon/hypergraph/vldb23NeighborCore/real/DBLP.hyp'
hscan_result_path='/home/hnu/Disk0/ParDon/Code/HSCAN/scan_hypergraph/output_csv/cluster_optimize_gs_index_dblp'
scan_result_path='/home/hnu/Disk0/ParDon/hypergraph/Hyper2General/DBLP/result'
wscan_result_path='/home/hnu/Disk0/ParDon/Github/Hyper-SCAN/weightSCAN/res/dblp'
output_file='dblp-new.txt'

hyperedges,gs=read_graph(dataset_file)
vertex=construct_link_index(hyperedges,gs)

hscan_result_file=hscan_result_path+'_'+mu+'_'+epsilon_hscan+'.csv'
scan_result_file=scan_result_path+'-'+epsilon_scan+'-'+mu+'.txt'
wscan_result_file=wscan_result_path+'-'+epsilon_wscan+'-'+mu+'.txt'

hscan_ecluster,hscan_vcluster=load_hscan_cluster(hscan_result_file,hyperedges)
scan_ecluster,scan_vcluster=load_scan_cluster(scan_result_file,hyperedges)
wscan_ecluster,wscan_vcluster=load_wscan_cluster(wscan_result_file,hyperedges)

modularity_hscan=compute_modularity_2(hyperedges,vertex,hscan_ecluster,hscan_vcluster)
modularity_scan=compute_modularity_2(hyperedges,vertex,scan_ecluster,scan_vcluster)
modularity_wscan=compute_modularity_2(hyperedges,vertex,wscan_ecluster,wscan_vcluster)
print('====================================')
print('dataset: '+output_file+' epsilon: '+epsilon_scan+' mu: '+mu)
print('modularity(hscan): '+str(modularity_hscan))
print('modularity(scan): '+str(modularity_scan))
print('modularity(wscan): '+str(modularity_wscan))
ofile=open(output_file,'a')
ofile.write(' epsilon: '+epsilon_scan+' mu: '+mu+' modularity(hscan): '+str(modularity_hscan)+' modularity(scan): '+str(modularity_scan)+' modularity(wscan): '+str(modularity_wscan))
# ofile.write(' epsilon: '+epsilon_scan[i]+' mu: '+mu[j]+' modularity(wscan): '+str(modularity_wscan))
ofile.write('\n')
ofile.close()