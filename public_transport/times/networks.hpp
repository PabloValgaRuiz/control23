#pragma once

#include <cstddef>
#include <vector>

struct Link{
    //Number of population that the link contains
    int Pop;
    //Number of population that have come before this link among the ones selected
    int cumulativePop;
};
struct Result{
    double mean{0};
    double mean2{0};

    int population_link{0};
};

struct MobilityTransport;
struct RhoMatrix{
    size_t I{};
    size_t J{};
    double value{};
    RhoMatrix(size_t i, size_t j, double val){
        I = i;
        J = j;
        value = val;
    }
    bool operator<(const RhoMatrix& rho) const{
        if((I < rho.I) || ((I == rho.I) && (J < rho.J)))
            return true;
        else return false;
    }
    bool operator<(const MobilityTransport& M) const;
    //bool operator==(const MobilityTransportMatrix& M);
};
struct CriteriaPopMatrix{
    size_t I{};
    size_t J{};
    double value{};
    double population{};
    CriteriaPopMatrix(RhoMatrix criteria_IJ, RhoMatrix n_IJ){
        I = criteria_IJ.I;
        J = criteria_IJ.J;
        value = criteria_IJ.value;
        population = n_IJ.value;
    }
};
struct MobilityTransport{
    size_t i{};
    size_t j{};
    size_t I{};
    size_t J{};
    size_t population{};
    size_t cumulative{};
    bool operator<(const MobilityTransport& M) const{

        if (i < M.i) return true;
        else if (i > M.i) return false;

        else if (j < M.j) return true;
        else if (j > M.j) return false;

        else if (I < M.I) return true;
        else if (I > M.I) return false;

        else if (J < M.J) return true;
        else if (J > M.J) return false;

        else return false;

    }
    bool operator<(const RhoMatrix& rho) const{
        if((I < rho.I) || ((I == rho.I) && (J < rho.J)))
            return true;
        else return false;
    }
    bool operator==(const MobilityTransport& M) const{
        if((I == M.I) && (J == M.J) && (i == M.i) && (j == M.j)) return true;
        else return false;
    }
};


// bool RhoMatrix::operator==(const MobilityTransportMatrix& M){
//     if((I == M.I) && (J == M.J)) return true;
//     else return false;
// }

typedef MobilityTransport MobTr;
typedef std::vector<MobTr> MobTrMatrix;


std::vector<RhoMatrix> betweennessCentrality(const std::vector<RhoMatrix>& n_IJ);
std::vector<RhoMatrix> eigenCentrality(const std::vector<RhoMatrix>& n_IJ);
std::vector<RhoMatrix> transportHubPopulation(const std::vector<RhoMatrix>& n_IJ);