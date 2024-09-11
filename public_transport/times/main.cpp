#include <iostream>
#include "networks.hpp"
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


Sparse<double> readEigen(const MobMatrix& T, const std::string& name);
Sparse<double> readWeights(const MobMatrix& T, const std::string& name);
std::tuple<MobTrMatrix, std::vector<RhoMatrix>, std::vector<RhoMatrix>> 
    orderLines(const MobMatrix& T, const Sparse<double>& eigenVector, const Sparse<double>& W);
MobTrMatrix chooseLinks(int NlinksChosen, const MobMatrix& T, const MobTrMatrix& n, const std::vector<RhoMatrix>& rho);
std::vector<double> chooseNumberOfLinks(const std::vector<RhoMatrix>& n_IJ, const MobMatrix& T, const MobTrMatrix& n, const std::vector<RhoMatrix>& CRITERIA_IJ);
void iterations(const MobMatrix& T, const std::vector<MobTrMatrix>& chosenLinks, std::vector<Result>& infected_results, std::vector<Result>& time_results, std::mt19937& generator);
std::vector<int> fisherYatesShuffle(int k, std::vector<int> range, std::mt19937& generator);
void countPopulationLinks(const MobMatrix& T, const std::vector<MobTrMatrix>& chosenLinks, std::vector<Result>& results);

static const int nIterations = 24 * 8; //24 * 8
static const double p = 1.0;

static const std::string name = "bogota"; //beta_c de bogota en p=1: 1/20.6942, 1/1.78102 para areas de ZATs

static const int MUESTRA_MAX = 20000;
static const int STEPS = 10000;
static const int THRESHOLD = 1000;
static const double beta = 3.0 / 1.78102;
static const int nPasos = 30;

std::mutex resultsMutex;

int main(){

    ThreadPool pool{24};

    std::string output = path + "out/eigen_centrality/" + name + "_transport_"+std::to_string(MUESTRA_MAX/1000)+"k_"+std::to_string(nPasos)+"d_beta_3,0_nu_1,00.txt";

    MobMatrix T{path + "bogota/mobnetwork.txt", path + "bogota/Poparea.txt"};

    auto eigenVector = readEigen(T, path + "eigenvectors3/bogota_10.txt");

    auto W = readWeights(T, path + "bogota/weights_800m.txt");
 
    //________________________________CHOOSING LINKS_________________________________
    MobTrMatrix n;
    std::vector<RhoMatrix> n_IJ; //USE FOR THE POPULATION CURVE
    std::vector<RhoMatrix> rho; //USE FOR THE TRANSPORT CURVE
    std::tie(n, n_IJ, rho) = orderLines(T, eigenVector, W);
    auto transportHubPop_IJ = transportHubPopulation(n_IJ);
    auto betweenness_centrality_IJ = betweennessCentrality(n_IJ);
    auto eigen_centrality_IJ = eigenCentrality(n_IJ);

    std::vector<RhoMatrix>& CRITERIA_IJ = eigen_centrality_IJ; //SELECT CRITERIA HERE

                                        //Do this before sorting the criteria vector
    const std::vector<double> NlcVector = chooseNumberOfLinks(n_IJ, T, n, CRITERIA_IJ);
    const size_t sizeLinks = NlcVector.size();

    //Sort the vector from higher to lower
    std::sort(CRITERIA_IJ.begin(), CRITERIA_IJ.end(), [](RhoMatrix a, RhoMatrix b){return a.value > b.value;});
    std::sort(n.begin(), n.end(), [](MobTr a, MobTr b){ return ((a.I < b.I) || ((a.I == b.I) && (a.J < b.J))) ? true : false; } );

    std::vector<Result> infected_results(sizeLinks);
    std::vector<Result> time_results(sizeLinks);
    std::vector<MobTrMatrix> vectorChosenLinks(sizeLinks);

    std::vector<std::future<void>> futures;
    for(size_t i = 0; i < NlcVector.size(); ++i){
        futures.push_back(std::move(pool.enqueue([&, i]{
            //FOR EIGENVECTORS USE rho, FOR POPULATION USE n_IJ
            vectorChosenLinks[i] = chooseLinks(NlcVector[i], T, n, CRITERIA_IJ);

        })));
	} //4996
    
	for(auto& future : futures){
        future.wait();
    }
    futures.clear();

    //_______________________________COUNTING PEOPLE___________________________________

    countPopulationLinks(T, vectorChosenLinks, infected_results);
    countPopulationLinks(T, vectorChosenLinks, time_results);

    //__________________________________ITERATING______________________________________

    for(int l = 0; l < nIterations; l++){
        futures.push_back(std::move(pool.enqueue([&, l]{
            auto generator = std::unique_ptr<std::mt19937>{new std::mt19937{std::random_device{}()}};
            iterations(T, vectorChosenLinks, infected_results, time_results, *generator);
        })));
    }

    for(auto& future : futures)
        future.wait();
    
    std::ofstream f(output);
    f << "population" << "\t" << "links" << "\t" << "tests" << "\t" << "detected" << "\t" << "time" << "\t" << "error" << "\t" << "time_error" << "\n";
	for(int i = 0; i < infected_results.size(); i++){//iteracion sobre links
            //Cantidad de links: Copia de mas arriba, al elegir los links
		// f << heatMap[i][j].population_link << "\t" << ((j+1) * MUESTRA_MAX) * nPasos / heatMap[i].size() << "\t" << heatMap[i][j].mean/nIterations << "\t" <<
        //     2 * sqrt((heatMap[i][j].mean2 - (heatMap[i][j].mean * heatMap[i][j].mean / (nIterations * nPasos)))/((nIterations * nPasos) * ((nIterations * nPasos)-1))) << "\n";
        f << infected_results[i].population_link << "\t"
		  << vectorChosenLinks[i].size() << "\t"
          << MUESTRA_MAX * nPasos << "\t"
          << infected_results[i].mean/nIterations << "\t"
          << time_results[i].mean/nIterations << "\t"
          //Std deviation
          << 1.96 * sqrt((infected_results[i].mean2 - (infected_results[i].mean * infected_results[i].mean / nIterations))/nIterations) << "\t"
          //Std error
          //<< 2 * sqrt((results[i].mean2 - (results[i].mean * results[i].mean / (nIterations)))/((nIterations) * ((nIterations)-1))) << "\n";
          << 1.96 * sqrt((time_results[i].mean2 - (time_results[i].mean * time_results[i].mean / nIterations))/nIterations) << std::endl;

          std::cout << infected_results[i].mean2 << "\t" << infected_results[i].mean << "\t" << nIterations << std::endl;
          std::cout << time_results[i].mean2 << "\t" << time_results[i].mean << "\t" << nIterations << std::endl;

    }
    f.close();

    return 0;
}


void iterations(const MobMatrix& T, const std::vector<MobTrMatrix>& chosenLinks, std::vector<Result>& infected_results, std::vector<Result>& time_results, std::mt19937& generator){
    MC_DistDiaNocheF montecarlo{0, p, T};
    montecarlo.setLambda(beta);
	montecarlo.inicializar(0.0001);

	std::vector<int> infected_means(chosenLinks.size(), 0);
    std::vector<int> time_means(chosenLinks.size(), nPasos);

    for (int t = 0; t < nPasos; t++){
    {
    mainPROFILE_SCOPE("Iteration");
		montecarlo.iteracion(T);
    }
		Log::info("Iteración: " + std::to_string(t));

        for(int l = 0; l < chosenLinks.size(); l++){//iteracion sobre links
        mainPROFILE_SCOPE("Iteration over tests");
            //SELECT THE POOL OF PEOPLE THAT HAVE USED PUBLIC TRANSPORT
            size_t totalPopChosen = 0;
            std::vector<int> MCLabels; //Pool
            for(const auto& link : chosenLinks[l])
                // if(link.I != link.J)
                for(int k = 0; k < link.population; k++)
                    if(montecarlo.getDesplazamiento()[link.cumulative + k] != montecarlo.getOrg()[link.cumulative + k]){
                        MCLabels.push_back(link.cumulative + k);
                    }
            //SHUFFLE
            const std::vector<int> sample{fisherYatesShuffle(MUESTRA_MAX, MCLabels, generator)};

            //Take pieces from the sample
            int PobInf = 0;
            for(int k = 0; k < MUESTRA_MAX && k < sample.size(); k++)
                if(montecarlo.getEst()[sample[k]] == E || montecarlo.getEst()[sample[k]] == I)
                    PobInf++;

            infected_means[l] += PobInf;
            if(infected_means[l] >= THRESHOLD && t < time_means[l])
                time_means[l] = t;
            //std::lock_guard<std::mutex> lock(heatMapMutex);
            //Log::debug(std::to_string(chosenLinks[i].Links) + "\t" + std::to_string(MUESTRA_MAX * nPasos) + "\t" + std::to_string(PobInf));
        }
    }

    std::lock_guard<std::mutex> lock(resultsMutex);
    for(int l = 0; l < chosenLinks.size(); l++){
        time_results[l].mean += static_cast<double>(time_means[l]);
        time_results[l].mean2 += static_cast<double>(time_means[l]) * static_cast<double>(time_means[l]);
        infected_results[l].mean += static_cast<double>(infected_means[l]);
        infected_results[l].mean2 += static_cast<double>(infected_means[l]) * static_cast<double>(infected_means[l]);
    }
}

void countPopulationLinks(const MobMatrix& T, const std::vector<MobTrMatrix>& chosenLinks, std::vector<Result>& results){
    size_t temp;
    for(int i = 0; i < chosenLinks.size(); i++){
        temp = 0;
        for(const auto& link : chosenLinks[i]){
            temp += link.population;
        }
        results[i].population_link = temp;
    }
}

std::vector<double> chooseNumberOfLinks(const std::vector<RhoMatrix>& n_IJ, const MobMatrix& T, const MobTrMatrix& n, const std::vector<RhoMatrix>& CRITERIA_IJ){

    std::vector<double> NlcVector;

    //n_IJ has the population of each link I,J but it might not be ordered as CRITERIA_IJ, which could be n_IJ or rho_IJ.
    //Create the matrix that will be ordered as the criteria but has population information
    std::vector<CriteriaPopMatrix> CRITERIA_pop_IJ;
    for(int i = 0; i < n_IJ.size(); i++){
        CRITERIA_pop_IJ.emplace_back(CRITERIA_IJ[i], n_IJ[i]);
    }
    std::sort(CRITERIA_pop_IJ.begin(), CRITERIA_pop_IJ.end(), [](CriteriaPopMatrix a, CriteriaPopMatrix b){return a.value > b.value;});

    //print to debug
    for(int i = 0; i < 10; i++){
        std::cout << CRITERIA_pop_IJ[i].I << " " << CRITERIA_pop_IJ[i].J << " " << CRITERIA_pop_IJ[i].value << " " << CRITERIA_pop_IJ[i].population << std::endl;
    }

    double totalPop = 0;
    for(int i = 0; i < CRITERIA_pop_IJ.size(); i++){
        const auto& link = CRITERIA_pop_IJ[i];
        if(link.I != link.J){ //If they go from I to I, they are not using public transport and cannot be tested
            totalPop += link.population;
            if(totalPop >= MUESTRA_MAX + 1000){
                NlcVector.push_back(i+1);
                std::cout << i+1 << ", " << totalPop << std::endl;
                break;
            }
        }
    }

    int step = 10000;
    int threshold = totalPop + step;
    for(int i = NlcVector[0]; i < CRITERIA_pop_IJ.size(); i++){
        const auto& link = CRITERIA_pop_IJ[i];
        if(link.I != link.J){
            totalPop += link.population;
            if(totalPop >= threshold){
                NlcVector.push_back(i+1);
                threshold = totalPop + step;
                std::cout << i+1 << ", " << totalPop << std::endl;
            }
        }
    }
    return NlcVector;
    // n_IJ.size() * (i+1 + offset) / (links_zoom * sizeLinks + 1 + offset);
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

//Sort the lines from I to J based on the eigenvector and return the triplet
std::tuple<MobTrMatrix, std::vector<RhoMatrix>, std::vector<RhoMatrix>>
    orderLines(const MobMatrix& T, const Sparse<double>& eigenVector, const Sparse<double>& W /*Links-Lines (i-I)*/ ){
    //Agents from link i->j, that use line I->J (n_jJ^iI), we need to keep this vector
    MobTrMatrix n{};

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
            if(it_last->cumulative + it_last->population - n[pos_first].cumulative != T.Mpesos[i][j]){
                std::cout << i << " " << T.Mvecinos[i][j] << std::endl;
                std::cout << "range: " << &(*it_last) - &(n[pos_first]) << std::endl;
                std::cout << "ntemp_ij: " << ntemp_ij << " , difference: " << difference << std::endl;
                std::cout << (it_last->cumulative + it_last->population - n[pos_first].cumulative) << " vs " << T.Mpesos[i][j] << std::endl;

                Log::error("Failed to calculate n_jJ^iI: the population doesn't add up.");
            }
        }

    }

    //create a list of links for rho_IJ and n_IJ WITHOUT SORTING IT

    std::vector<RhoMatrix> rho_IJ{};
    std::vector<RhoMatrix> n_IJ_linklist{};

    for(int I = 0; I < n_IJ.N; I++){
        for(int J = 0; J < n_IJ.vecinos[I]; J++){
            n_IJ_linklist.emplace_back(static_cast<size_t>(I), static_cast<size_t>(n_IJ.Mvecinos[I][J]), n_IJ.Mpesos[I][J]);
            if(n_IJ[I][J] != 0)
                rho_IJ.emplace_back(static_cast<size_t>(I), static_cast<size_t>(n_IJ.Mvecinos[I][J]), I_IJ[I][J] / n_IJ[I][J]);
            else rho_IJ.emplace_back(static_cast<size_t>(I), static_cast<size_t>(n_IJ.Mvecinos[I][J]), 0);
        }
    }

    return {n, n_IJ_linklist, rho_IJ};
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


    return n_cut;
}
