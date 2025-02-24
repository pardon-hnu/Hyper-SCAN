#include <iostream>
#include "weightedSCAN.h"

using namespace std;

int main(int argc, char *argv[]) 
{
    string dataset=argv[1];
    WSCAN worker(dataset);
    worker.construct_index();
    
    // vector<int> mus={2,3,4,5,10,15};
    vector<int> mus={6,7,8,9};
    vector<double> epsilons={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
    worker.batch_pSCAN(mus,epsilons);
    // worker.pSCAN(2,0.6);
    return 0;
}