#include <iostream>
#include "MobMatrix.hpp"
#include "MC_DistDiaNocheF.hpp"
#include "ThreadPool.hpp"
#include <string>
#include <sstream>
#include <random>
#include <cstring>
#include <algorithm>
#include <unordered_map>
#include "Sparse.hpp"
#include "Log.hpp"
#include "benchmark.hpp"

#define MAIN_PROFILE_FUNCTION 1

const std::string path = "../";

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

bool RhoMatrix::operator<(const MobilityTransport& M) const{
    if((I < M.I) || ((I == M.I) && (J < M.J)))
        return true;
    else return false;
}
// bool RhoMatrix::operator==(const MobilityTransportMatrix& M){
//     if((I == M.I) && (J == M.J)) return true;
//     else return false;
// }

typedef MobilityTransport MobTr;
typedef std::vector<MobTr> MobTrMatrix;

Sparse<double> readEigen(const MobMatrix& T, const std::string& name);
Sparse<double> readWeights(const MobMatrix& T, const std::string& name);
std::pair<MobTrMatrix, std::vector<RhoMatrix>> 
    orderLines(const MobMatrix& T, const Sparse<double>& eigenVector, const Sparse<double>& W);
MobTrMatrix chooseLinks(int NlinksChosen, const MobMatrix& T, const MobTrMatrix& n, const std::vector<RhoMatrix>& rho);
void iterations(const MobMatrix& T, const std::vector<MobTrMatrix>& chosenLinks, std::vector<std::vector<Result>>& heatMap, std::mt19937& generator);
std::vector<int> fisherYatesShuffle(int k, std::vector<int> range, std::mt19937& generator);
void countPopulationLinks(const MobMatrix& T, const std::vector<MobTrMatrix>& chosenLinks, std::vector<std::vector<Result>>& heatMap);

#define MUESTRA_MAX 1000

std::string name = "bogota";
                        //beta_c de bogota en p=1: 1/20.6942, 1/1.78102 para areas de ZATs
static double beta = 8.0 / 1.78102;
static double p = 1.0;
static int nPasos = 30;
static int nIterations = 24 * 1; //24 * 32

std::mutex heatMapMutex;

int main(){

    ThreadPool pool{24};

    std::string output = path + "out/" + name + "_transport_30k_beta_8,0.txt";

    MobMatrix T{path + "bogota/mobnetwork.txt", path + "bogota/Poparea.txt"};

    auto eigenVector = readEigen(T, path + "eigenvectors3/bogota_10.txt");

    auto W = readWeights(T, path + "bogota/weights_800m.txt");

    constexpr size_t sizeLinks = 257;
	constexpr size_t sizeTests = 1;

    std::vector<std::vector<Result>> heatMap;
	heatMap.resize(sizeLinks);
	for(auto& v : heatMap)
		v.resize(sizeTests);
    std::vector<MobTrMatrix> vectorChosenLinks;
	vectorChosenLinks.resize(sizeLinks);
 
    //________________________________CHOOSING LINKS_________________________________
    MobTrMatrix n;
    std::vector<RhoMatrix> rho;
    std::tie(n, rho) = orderLines(T, eigenVector, W);

    std::ofstream log("log.txt");
    for(auto r : rho){
        log << "I: " << r.I << ", J: " << r.J << ", value: " << r.value << std::endl;
    }
    log.close();

    std::sort(n.begin(), n.end(), [](MobTr a, MobTr b){
        if((a.I < b.I) || ((a.I == b.I) && (a.J < b.J)))
            return true;
        else return false;
    });

    std::vector<std::future<void>> futures;
    for(size_t i = 0; i < sizeLinks; ++i){
        futures.push_back(std::move(pool.enqueue([&, i]{
            size_t NlcTemp = rho.size() * (i+1) / (1 * sizeLinks + 1); //Care to not choose 0 or all the links (ex: if 1000 links, take from 166 to 866)
            //Choose the Nlc highest component links in the eigenvector
            vectorChosenLinks[i] = chooseLinks(NlcTemp, T, n, rho);
        })));
	}
	for(auto& future : futures){
        future.wait();
    }
    futures.clear();

    Log::debug("Links chosen.");

    std::ofstream log2("log2.txt");
    for(auto v : vectorChosenLinks[0]){
        log2 << "i: " << v.i << ", j: " << v.j << ", I: " << v.I << ", J: " << v.J << std::endl;
    }
    log2.close();
    //return 0;

    //_______________________________COUNTING PEOPLE___________________________________

    countPopulationLinks(T, vectorChosenLinks, heatMap);

    //__________________________________ITERATING______________________________________

    for(int l = 0; l < nIterations; l++){
        futures.push_back(std::move(pool.enqueue([&, l]{
            auto generator = std::unique_ptr<std::mt19937>{new std::mt19937{std::random_device{}()}};
            iterations(T, vectorChosenLinks, heatMap, *generator);
        })));
    }

    for(auto& future : futures)
        future.wait();
    
    std::ofstream f(output);
	for(int i = 0; i < heatMap.size(); i++)//iteracion sobre links
		for(int j = 0; j < heatMap[i].size(); j++)//iteracion sobre tests
                //Cantidad de links: Copia de mas arriba, al elegir los links
			f << heatMap[i][j].population_link << "\t" << ((j+1) * MUESTRA_MAX) * nPasos / heatMap[i].size() << "\t" << heatMap[i][j].mean/nIterations << "\t" <<
            2 * sqrt((heatMap[i][j].mean2 - (heatMap[i][j].mean * heatMap[i][j].mean / (nIterations * nPasos)))/((nIterations * nPasos) * ((nIterations * nPasos)-1))) << "\n";
	f.close();

    return 0;
}


void iterations(const MobMatrix& T, const std::vector<MobTrMatrix>& chosenLinks, std::vector<std::vector<Result>>& heatMap, std::mt19937& generator){
    MC_DistDiaNocheF montecarlo{0, p, T};
    montecarlo.setLambda(beta);
	montecarlo.inicializar(0.0001);
	int PobInf;

    for (int t = 0; t < nPasos; t++){
    {
    mainPROFILE_SCOPE("Iteration");
		montecarlo.iteracion(T);
    }
		Log::info("Iteración: " + std::to_string(t));

        for(int i = 0; i < heatMap.size(); i++){//iteracion sobre links
        mainPROFILE_SCOPE("Iteration over tests");
            //SELECT THE POOL OF PEOPLE THAT HAVE USED PUBLIC TRANSPORT
            size_t totalPopChosen = 0;
            std::vector<int> MCLabels; //Pool
            for(const auto& link : chosenLinks[i])
                if(link.I != link.J)
                    for(int k = 0; k < link.population; k++)
                        if(montecarlo.getDesplazamiento()[link.cumulative + k] != montecarlo.getOrg()[link.cumulative + k]){
                            MCLabels.push_back(link.cumulative + k);
                        }
            //SHUFFLE
            const std::vector<int> sample{fisherYatesShuffle(MUESTRA_MAX, MCLabels, generator)};

            //Take pieces from the sample
            for(int j = 0; j < heatMap[i].size(); j++){//iteracion sobre tests

                int MUESTRA = ((j+1) * MUESTRA_MAX) / heatMap[i].size();
                
                PobInf = 0;
                for(int k = 0; k < MUESTRA && k < sample.size(); k++)
                    if(montecarlo.getEst()[sample[k]] == E)
                        PobInf++;

                std::lock_guard<std::mutex> lock(heatMapMutex);
                heatMap[i][j].mean += static_cast<double>(PobInf);
                heatMap[i][j].mean2 += static_cast<double>(PobInf) * static_cast<double>(PobInf);
            }
            //std::lock_guard<std::mutex> lock(heatMapMutex);
            //Log::debug(std::to_string(chosenLinks[i].Links) + "\t" + std::to_string(MUESTRA_MAX * nPasos) + "\t" + std::to_string(PobInf));
        }
    }
}

void countPopulationLinks(const MobMatrix& T, const std::vector<MobTrMatrix>& chosenLinks, std::vector<std::vector<Result>>& heatMap){
    size_t temp;
    for(int i = 0; i < chosenLinks.size(); i++){
        temp = 0;
        for(const auto& link : chosenLinks[i]){
            temp += link.population;
        }
        for(auto& point : heatMap[i]){
            point.population_link = temp;
        }
    }
}

std::vector<int> fisherYatesShuffle(int k, std::vector<int> range, std::mt19937& generator){

    int n = range.size();

    if(k > n){
        Log::error("Number of sampled people larger than the population in the links: " + std::to_string(k) + " > " + std::to_string(n));
        k = n;
    }
    std::vector<int> reservoir; reservoir.resize(k);

    //Fisher-Yates shuffle is used to shuffle completely a vector.
    //However, we just need the lower k components and after it gets to the k-1th element it won't change the
    //order of the ones before it, so we can stop it there
    for(int i = 0; i < k; i++){
        //Choose a random number from i to n-1
        int j = std::uniform_int_distribution<>{i,n-1}(generator);
        //Swap range[i] with range[j]
        int temp = range[i];
        range[i] = range[j];
        range[j] = temp;
    }

    for(int i = 0; i < k; i++){
        reservoir[i] = range[i];
    }

    return reservoir;
}

Sparse<double> readEigen(const MobMatrix& T, const std::string& name){
mainPROFILE_FUNCTION();

    Sparse<double> eigenVector = T;

    std::ifstream fileEigen{name};

    Log::debug("Eigenvector file opened.");

    int a, b; double c;

    for(int i = 0; i < T.N; i++){
        for(int j = 0; j < T.vecinos[i]; j++){
            while(fileEigen >> a >> b >> c){
                if(a == i && b == T.Mvecinos[i][j]){
                    eigenVector[i][j] = c;
                    break;
                }
            }
        }
    }

    Log::debug("Eigenvector read and created.");
    fileEigen.close();
    return eigenVector;
}
Sparse<double> readWeights(const MobMatrix& T, const std::string& name){
mainPROFILE_FUNCTION();

    std::ifstream fileWeights{name};
    if(fileWeights)
        Log::debug("Weights file opened.");

    int a, b; double c;

    int N = 0, M = 0;
    while(fileWeights >> a >> b >> c){
        if(a >= N){
            N = a + 1; //Not using this
        }
        if(b >= M){
            M = b + 1;
        }
    }
    fileWeights.close();

    //When some people don't use public transport, some nodes may not appear in the weights file, so we use the ones of the mobmatrix
    Sparse<double> Weights(T.N, M);

    fileWeights.open(name);
    while(fileWeights >> a >> b >> c){
        Weights.insert(a, b, c);
    }

    Log::debug("Weights read and created.");
    std::cout << "N: " << T.N << ", M: " << M << std::endl;
    fileWeights.close();
    return Weights;
}

//Order the lines from I to J based on the eigenvector
std::pair<MobTrMatrix, std::vector<RhoMatrix>> 
    orderLines(const MobMatrix& T, const Sparse<double>& eigenVector, const Sparse<double>& W /*Links-Lines (i-I)*/ ){
    //Agents from link i->j, that use line I->J (n_jJ^iI), we need to keep this vector
    MobTrMatrix n{};

    MC_DistDiaNocheF montecarlo{0, p, T};
    montecarlo.setLambda(beta);
	montecarlo.inicializar(0.0001);

    //I_jJ^iI is not necessary to calculate since we only need I_IJ which we can calculate directly
    Sparse<double> I_IJ(W.M, W.M);
    Sparse<double> n_IJ(W.M, W.M);
    
    //If not everyone uses public transport, the 'cumul' variable will not reach the entirety of n_ij, so we have to update it at the
    //end of the i,j loop. The people that are not reached are the ones that will not use public transport
    size_t cumul_link = 0;
    size_t cumul = 0;
    for(size_t i = 0; i < T.N; i++)
    for(size_t j = 0; j < T.vecinos[i]; j++){
        //Loop on the members of n_ij
        //We have interest in knowing its first element (the next to the last of the previous iteration)
        int pos_first = n.end() - n.begin();
        for(size_t I = 0; I < W.vecinos[i]; I++){
        for(size_t J = 0; J < W.vecinos[T.Mvecinos[i][j]]; J++){
            if(W.Mpesos[i][I] != 0 && W.Mpesos[T.Mvecinos[i][j]][J] != 0){//if W_iI and W_jJ are not zero
            
                double temp = W.Mpesos[i][I] * W.Mpesos[T.Mvecinos[i][j]][J] * T.Mpesos[i][j];
                //insert the element n_jJ^iI and its cumulative number
                n.push_back(MobTr{
                i,
                static_cast<size_t>(T.Mvecinos[i][j]),
                static_cast<size_t>(W.Mvecinos[i][I]),
                static_cast<size_t>(W.Mvecinos[T.Mvecinos[i][j]][J]),
                static_cast<size_t>(temp),
                cumul});
                //Da esto la posicion en el montecarlo?? Si todos los demas elementos son cero, debería
                cumul += static_cast<size_t>(temp);
                //Find if the line I->J is in I_IJ
                //                        Row I                                                                           This is just J
                auto iterator = std::find(I_IJ.Mvecinos[W.Mvecinos[i][I]].begin(), I_IJ.Mvecinos[W.Mvecinos[i][I]].end(), W.Mvecinos[T.Mvecinos[i][j]][J]);
                //if it is there, add the infected
                
                if(iterator != I_IJ.Mvecinos[W.Mvecinos[i][I]].end()){
                    int position = iterator - I_IJ.Mvecinos[W.Mvecinos[i][I]].begin(); //position in the vector so we can use it in Mpesos
                    I_IJ.Mpesos[W.Mvecinos[i][I]][position] += temp * eigenVector.Mpesos[i][j];
                    n_IJ.Mpesos[W.Mvecinos[i][I]][position] += temp;
                }
                //if not, add the line and then the infected
                else{
                    int _I = W.Mvecinos[i][I];
                    int _J = W.Mvecinos[T.Mvecinos[i][j]][J];

                    I_IJ.insert(_I, _J, temp * eigenVector.Mpesos[i][j]);
                    n_IJ.insert(_I, _J, temp);

                }
            }
            }
        }

        //If there are 0 people in n_ij, we just move on
        if(T.Mpesos[i][j] == 0) continue;

        auto it_last = n.end() - 1;

        //What if sum n_jJ^iI does not equal n_ij? 
        bool notEveryoneUsesTransport = true;

        if(notEveryoneUsesTransport){
            cumul_link += T.Mpesos[i][j];
            cumul = cumul_link;
        }
        else{

            //We have to add those agents one by one. A single loop should do it, since
            //with every truncation we lose 1 agent at max
            int ntemp_ij = it_last->cumulative + it_last->population - n[pos_first].cumulative;
            int difference = T.Mpesos[i][j] - ntemp_ij; //this is the number of agents they disagree

            //Just add them in order
            auto it = n.begin()+pos_first;
            for(int l = 0; l < difference && it < n.end(); l++){
                (it->population)++;
                if(l != 0){
                    it->cumulative = (it-1)->cumulative + (it-1)->population;
                }
                ++it;
            }
            for(int l = difference; it != n.end(); l++){
                if(l != 0)
                    it->cumulative = (it-1)->cumulative + (it-1)->population;
                ++it;
            }
            //Remember that cumul was counting independently, we have to update it after solving the problem before
            cumul = n.back().cumulative + n.back().population;
            //There should be the same as in n_ij so let's make sure
            if(i == 399 && T.Mvecinos[i][j] == 126){}
            if(it_last->cumulative + it_last->population - n[pos_first].cumulative != T.Mpesos[i][j]){
                std::cout << i << " " << T.Mvecinos[i][j] << std::endl;
                std::cout << "range: " << &(*it_last) - &(n[pos_first]) << std::endl;
                std::cout << "ntemp_ij: " << ntemp_ij << " , difference: " << difference << std::endl;
                std::cout << (it_last->cumulative + it_last->population - n[pos_first].cumulative) << " vs " << T.Mpesos[i][j] << std::endl;

                Log::error("Failed to calculate n_jJ^iI: the population doesn't add up.");
            }
        }

    }

    std::vector<RhoMatrix> rho_IJ{};

    for(int I = 0; I < n_IJ.N; I++){
        for(int J = 0; J < n_IJ.vecinos[I]; J++){
            if(n_IJ[I][J] != 0)
                rho_IJ.emplace_back(static_cast<size_t>(I), static_cast<size_t>(n_IJ.Mvecinos[I][J]), I_IJ[I][J] / n_IJ[I][J]);
            else rho_IJ.emplace_back(static_cast<size_t>(I), static_cast<size_t>(n_IJ.Mvecinos[I][J]), 0);
        }
    }

    //Sort the vector from higher to lower
    std::sort(rho_IJ.begin(), rho_IJ.end(), [](RhoMatrix a, RhoMatrix b){return a.value > b.value;});

    return {n, rho_IJ};
}

MobTrMatrix chooseLinks(int NlinksChosen, const MobMatrix& T, const MobTrMatrix& n, const std::vector<RhoMatrix>& rho){
mainPROFILE_FUNCTION();

    if(NlinksChosen > rho.size()) Log::error("More chosen links than there are in the network");

    std::vector<RhoMatrix> rho_cut;
    std::copy(rho.begin(), rho.begin() + NlinksChosen, std::back_inserter(rho_cut));

    
    std::sort(rho_cut.begin(), rho_cut.end());
    MobTrMatrix n_cut{};

    auto first1 = n.begin();
    auto last1 = n.end();
    auto first2 = rho_cut.begin();
    auto last2 = rho_cut.end();
    auto d_first = std::back_inserter(n_cut);

    //MODIFIED STL intersection algorithm -> it will copy more than just the first appearance of a value on the first vector, even if
    //it only appears once in the second one
    while (first1 != last1 && first2 != last2) {
        if (*first1 < *first2) {
            ++first1;
        } else  {
            if (!(*first2 < *first1)) {
                *d_first++ = *first1++; // *first1 and *first2 are equivalent.
            }
            else{ //This else is the only modification of the intersection algorithm
                ++first2;
            }
        }
    }

    Log::debug("Links chosen");

    return n_cut;
}