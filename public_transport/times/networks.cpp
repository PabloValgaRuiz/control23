#include "networks.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

bool RhoMatrix::operator<(const MobilityTransport& M) const{
    if((I < M.I) || ((I == M.I) && (J < M.J)))
        return true;
    else return false;
}

void writeNetwork(const std::vector<RhoMatrix>& n_IJ){

    std::ofstream file("transport_n_IJ.txt");
    for(const auto& link : n_IJ){
        file << link.I << " " << link.J << " " << link.value << std::endl;
    }

}

void getBetweennessCentrality(const std::vector<RhoMatrix>& n_IJ){
    writeNetwork(n_IJ);

    std::string command = "python3 ../betweenness_centrality.py";
    system(command.c_str());
    std::cout << "Betweenness centrality calculated" << std::endl;
}

void getEigenCentrality(const std::vector<RhoMatrix>& n_IJ){
    writeNetwork(n_IJ);

    std::string command = "python3 ../eigen_centrality.py";
    system(command.c_str());
    std::cout << "Eigenvector centrality calculated" << std::endl;
}

std::vector<double> readBetweennessCentrality(const std::vector<RhoMatrix>& n_IJ){
    std::vector<double> centrality;
    std::string filename = "betweenness_centrality.txt";
    std::ifstream file(filename);

    while(!file){
        getBetweennessCentrality(n_IJ);
        file.open(filename);
    }

    double value;
    int i;
    while(file >> i >> value){
        centrality.push_back(value);
    }
    return centrality;
}

std::vector<double> readEigenCentrality(const std::vector<RhoMatrix>& n_IJ){
    std::vector<double> centrality;
    std::string filename = "eigenvector_centrality.txt";
    std::ifstream file(filename);

    while(!file){
        getEigenCentrality(n_IJ);
        file.open(filename);
    }

    double value;
    int i;
    while(file >> i >> value){
        centrality.push_back(value);
    }
    return centrality;
}

std::vector<RhoMatrix> betweennessCentrality(const std::vector<RhoMatrix>& n_IJ){
    auto centrality = readBetweennessCentrality(n_IJ);
    auto centrality_IJ = n_IJ;

    for(int l = 0; l < centrality_IJ.size(); l++){
        //                  VALUE OF CENTRALITY OF THE STATION OF DESTINATION, COULD BE CHANGED TO THE ORIGIN HUB
        centrality_IJ[l].value = centrality[n_IJ[l].J];
    }
    return centrality_IJ;
}

std::vector<RhoMatrix> eigenCentrality(const std::vector<RhoMatrix>& n_IJ){
    auto centrality = readEigenCentrality(n_IJ);
    auto centrality_IJ = n_IJ;

    for(int l = 0; l < centrality_IJ.size(); l++){
        //                  VALUE OF CENTRALITY OF THE STATION OF DESTINATION, COULD BE CHANGED TO THE ORIGIN HUB
        centrality_IJ[l].value = centrality[n_IJ[l].J];
    }
    return centrality_IJ;
}

std::vector<RhoMatrix> transportHubPopulation(const std::vector<RhoMatrix>& n_IJ){
    auto n_Hub_IJ = n_IJ;
    int M = 0;
    for(const auto& link : n_IJ){
        if(link.I >= M) M = link.I + 1;
    }
    
    std::vector<double> hubs_population(M, 0);

    for(const auto& link : n_IJ){
        hubs_population[link.J] += link.value;
    }
    for(auto& link : n_Hub_IJ){
        link.value = hubs_population[link.J];
    }
    return n_Hub_IJ;
}