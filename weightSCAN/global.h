#pragma once

#include <iostream>
#include <string>
#include <map>
using namespace std;
std::map <std::string,std::string> dataset_to_filename = {
            {"enron" , "/home/hnu/Disk0/ParDon/hypergraph/hyper2weight/Enron/edge.txt"},
            {"congress" , "/home/hnu/Disk0/ParDon/hypergraph/hyper2weight/congress-bills/edge.txt"},
            {"contact" , "/home/hnu/Disk0/ParDon/hypergraph/hyper2weight/contact-primary-school/edge.txt"},
            {"dblp", "/home/hnu/Disk0/ParDon/hypergraph/hyper2weight/DBLP/edge.txt"},
            {"aminer","/home/hnu/Disk0/ParDon/hypergraph/hyper2weight/aminer/edge.txt"},
            {"drug","/home/hnu/Disk0/ParDon/hypergraph/hyper2weight/NDC-substances/edge.txt"},
            {"ubuntu","/home/hnu/Disk0/ParDon/hypergraph/hyper2weight/threads-ask-ubuntu/edge.txt"},
            {"rpah","/home/hnu/Disk0/ParDon/hypergraph/hyper2weight/rpa_t130000000/edge.txt"},
            {"mag10","/home/hnu/Disk0/ParDon/hypergraph/hyper2weight/mag10/edge.txt"},
            {"coauth-DBLP","/home/hnu/Disk0/ParDon/hypergraph/hyper2weight/coauth-DBLP/edge.txt"}
        };

