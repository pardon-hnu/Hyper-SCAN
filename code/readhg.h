#ifndef READHG_H
#define READHG_H
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <iostream>
#include <iterator>
#include <vector>
#include "hypergraph.h"

std::map <std::string,std::string> dataset_to_filename = {
            {"enron" , "../datasets/Enron.hyp"},
            {"congress" , "../datasets/congress-bills.hyp"},
            {"contact" , "../datasets/contact-primary-school.hyp"},
            {"dblp", "../datasets/DBLP.hyp"},
            {"aminer","../datasets/aminer.hyp"},          
            {"drug","../datasets/NDC-substances.hyp"},
            {"ubuntu","../datasets/threads-ask-ubuntu.hyp"},
            {"rpah","../datasets/rpa_t130000000.hyp"}
        };

template <typename Out>
void split(const std::string &s, char delim, Out result) {
    std::istringstream iss(s);
    std::string item;
    while (std::getline(iss, item, delim)) {
        *result++ = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

void getHg(std::string dataname, Hypergraph & hg){
    std::cout<< dataname<<"\n";
    std::cout << dataset_to_filename[dataname]<<"\n";
    std::ifstream infile(dataset_to_filename[dataname]);
    std::string line;
    // std::map< int, std::vector<std::string>  > Edges;
    size_t i= 0;
    while (std::getline(infile, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        //std::cout<< i<<" "<<line<<"\n";
        std::vector<std::string> x = split(line ,',');   //识别逗号作为分隔符号，从字符串中（比如 1，2，3）提取出一个个超边中的顶点信息（1 和2 以及3）保存在x中
        hg.addEdge(i,x);//x存的是string类型的顶点编号
        i++;
    }
    // hg.print();
}


#endif
