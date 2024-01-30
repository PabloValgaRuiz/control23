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

struct Link{
    //Number of population that the link contains
    int Pop;
    //Number of population that have come before this link among the ones selected
    int cumulativePop;
};

struct Result{
    double mean;
    double mean2;
};

Sparse<double> readEigen(const MobMatrix& T, const std::string& name);
Sparse<Link> chooseLinks(int NlinksChosen, const MobMatrix& T, const Sparse<double>& eigenVector);
void iterations(const MobMatrix& T, const std::vector<Sparse<Link>>& chosenLinks, std::vector<std::vector<Result>>& results,
    std::vector<Result>& attackRate, double beta);
int contarInfectadosChosen(const MobMatrix& T, const Sparse<Link>& chosenLinks, MC_DistDiaNocheF& montecarlo, int MUESTRA, std::mt19937& generator);
std::vector<int> fisherYatesShuffle(int k, std::vector<int> range, std::mt19937& generator);
std::vector<int> reservoirSampling(int k, int n, std::mt19937& generator);
std::vector<double> countPopulationLinks(const MobMatrix& T, const std::vector<Sparse<Link>>& chosenLinks);

#define MUESTRA_MAX 200000 //200000
//Boston: 200k

const static std::unordered_map<std::string, double> cityBeta{
    {"baltimore", 0.318387},
    {"ny", 0.0466719},
    {"dallas", 0.025512},
    {"dc", 0.839636},
    {"detroit", 0.0106555},
    {"ma", 0.153499},
    {"los angeles", 0.0918512},
    {"miami", 0.168365},
    {"bogota", 1.0 / 1.78102}
};

std::string path = "../";
std::string name = "miami";

//static const double beta = 3.0 * cityBeta.at(name);
static const double p = 1.0;
static const int nPasos = 300;//300
static const int nIterations = 24*1;

std::mutex resultsMutex;
/*
int main(int argc, char* argv[]){

    Instrumentor::Get().BeginSession("Session Name");
{
    InstrumentationTimer timer("Program");

    ThreadPool pool{24};

    std::string output = name;// + "_beta_8,0.txt";

    MobMatrix T{path + "cities3/" + name + "/mobnetwork.txt", path + "cities3/" + name + "/Poparea.txt"};
    std::cout << T.Pob << std::endl;
    auto eigenVector = readEigen(T, name);

    size_t sizeLinks = 33; //33

    std::vector<std::vector<Result>> results(sizeLinks);
    for(auto& i : results)  i.resize(nPasos);
    std::vector<Result> attackRate(sizeLinks);
    std::vector<Sparse<Link>> vectorChosenLinks(sizeLinks);

    //________________________________CHOOSING LINKS_________________________________

    std::vector<std::future<void>> futures;
    for(size_t i = 0; i < sizeLinks; ++i){
        futures.push_back(std::move(pool.enqueue([&, i]{
            constexpr int offset = 3;
            size_t NlcTemp = T.Links * (i+1 + offset) / (2.8 * sizeLinks + offset); //Care to not choose 0
                                                        //Boston: 2.8 * sizeLinks + offset
            //Choose the Nlc highest component links in the eigenvector
            vectorChosenLinks[i] = chooseLinks(NlcTemp, T, eigenVector);
        })));
	}
	for(auto& future : futures){
        future.wait();
    }
    futures.clear();

    auto link_populations = countPopulationLinks(T, vectorChosenLinks);

    Log::debug("Links chosen.");

    //__________________________________ITERATING_____________________________________
    std::ofstream file(path + "out/attackRateMap/" + output + "beta_5.txt");
    file << "beta" << "\t" << "population" << "\t" << "links" << "\t" << "attackRate" << "\t" << "error" << "\n";

    double beta0 = cityBeta.at(name);
    double beta_diff = beta0/4;
    for(double beta = 5 * beta0; beta < 5.001 * beta0; beta += beta_diff){
        
        attackRate.assign(attackRate.size(), Result({0,0}));
        
        for(int l = 0; l < nIterations; l++){
            futures.push_back(std::move(pool.enqueue([&, l, beta]{

                iterations(T, vectorChosenLinks, results, attackRate, beta);

            })));
        }

        for(auto& future : futures)
        future.wait();
        for(int l = 0; l < vectorChosenLinks.size(); l++){
            file << beta / cityBeta.at(name) << "\t" << link_populations[l] << "\t" << vectorChosenLinks[l].Links << "\t" << attackRate[l].mean / nIterations << "\t"
                 << 1.96 * sqrt((attackRate[l].mean2 - attackRate[l].mean*attackRate[l].mean/nIterations)/(nIterations * (nIterations - 1)))
                 << "\n";
        }
    }
    file.close();
    
}
    Instrumentor::Get().EndSession();
    return 0;
}
*/

int main(int argc, char* argv[]){

    ThreadPool pool{24};

    std::string output = name;// + "_beta_8,0.txt";

    MobMatrix T{path + "cities3/" + name + "/mobnetwork.txt", path + "cities3/" + name + "/Poparea.txt"};
    std::cout << T.Pob << std::endl;
    auto eigenVector = readEigen(T, name);


    //________________________________CHOOSING LINKS_________________________________

    const size_t NlcTemp = 100; //Number of links chosen
    Sparse<Link> chosenLinks = chooseLinks(NlcTemp, T, eigenVector);

    std::ofstream outfile(path + "../figures/links_map/chosen_links/" + std::to_string(NlcTemp) + "_" + name + ".txt");
    for(int i = 0; i < chosenLinks.N; i++){
        for(int j = 0; j < chosenLinks.vecinos[i]; j++){
            outfile << i << "\t" << chosenLinks.Mvecinos[i][j] << "\t" << chosenLinks[i][j].Pop << "\n";
        }
    }

    Log::debug("Links chosen.");

    return 0;
}


void iterations(const MobMatrix& T, const std::vector<Sparse<Link>>& chosenLinks, std::vector<std::vector<Result>>& results,
    std::vector<Result>& attackRate, double beta){
    
    std::mt19937 mt(std::random_device{}());

    for(int l = 0; l < chosenLinks.size(); l++){
        MC_DistDiaNocheF montecarlo{0, p, T};
        montecarlo.setLambda(beta);
	    montecarlo.inicializar(0.0004);

        for(int t = 0; t < nPasos; t++){
            montecarlo.iteracion(T);
            int pobInf = contarInfectadosChosen(T, chosenLinks[l], montecarlo, MUESTRA_MAX, mt);
            
            std::lock_guard<std::mutex> lock(resultsMutex);
            results[l][t].mean += montecarlo.totalI;
            results[l][t].mean2 += static_cast<double>(montecarlo.totalI) * montecarlo.totalI;
            
        }
        std::lock_guard<std::mutex> lock(resultsMutex);
        attackRate[l].mean += montecarlo.totalR;
        attackRate[l].mean2 += static_cast<double>(montecarlo.totalR) * montecarlo.totalR;
    }
}

int contarInfectadosChosen(const MobMatrix& T, const Sparse<Link>& chosenLinks, MC_DistDiaNocheF& montecarlo, int MUESTRA, std::mt19937& generator){

    std::vector<int> MCLabels(chosenLinks.Total);

{
mainPROFILE_SCOPE("Copying");
    //Accounts for the number of people included in the reserve(vector MCLabels)
    int statesIndex = 0;

    for(int i = 0; i < chosenLinks.N; i++){
        for(int j = 0; j < chosenLinks.vecinos[i]; j++){
            //Find the first person in the link in the Monte Carlo agent array
            int MCIndex = chosenLinks[i][j].cumulativePop;

            //Copy the adress of every person included in the reserve
            for(int k = 0; k < chosenLinks[i][j].Pop; k++)
                MCLabels[statesIndex + k] = MCIndex + k;

            //Add the number of included people
            statesIndex += chosenLinks[i][j].Pop;
        }
    }
}

{
mainPROFILE_SCOPE("Shuffle");
    //Shuffle the labels, take the first MUESTRA number of them
    MCLabels = fisherYatesShuffle(MUESTRA, MCLabels, generator);
}
    int pobInfChosen = 0;
    //Test the sampled population and isolate the infected, ie send them to the R compartment
    for(auto label : MCLabels){
        if(montecarlo.getEst()[label] == E || montecarlo.getEst()[label] == I){
            pobInfChosen++;
            montecarlo.getEst()[label] = R;
        }
    }

    return pobInfChosen;
}

std::vector<int> fisherYatesShuffle(int k, std::vector<int> range, std::mt19937& generator){

    int n = range.size();

    if(k > n){
        std::cout << "Number of sampled people larger than the population in the links." << std::endl;;
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

    std::ifstream fileEigen{path + "eigenvectors3/" + name + "_10.txt"};

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

    return eigenVector;

}

Sparse<Link> chooseLinks(int NlinksChosen, const MobMatrix& T, const Sparse<double>& eigenVector){
mainPROFILE_FUNCTION();

    if(NlinksChosen > T.Links) Log::error("More chosen links than there are in the network");

    Sparse<Link> chosenLinksPop; chosenLinksPop.N = eigenVector.N;
    chosenLinksPop.vecinos.resize(T.N);
    chosenLinksPop.Mvecinos.resize(T.N);
    chosenLinksPop.Mpesos.resize(T.N);

    auto eigenTemp = eigenVector;
    int I, J; double prov;
    for(int k = 0; k < NlinksChosen; k++){

        //Find maximum component
        
        I = J = prov = 0;
        for(int i = 0; i < eigenTemp.N; i++){
            for(int j = 0; j < eigenTemp.vecinos[i]; j++){
                if(prov < eigenTemp[i][j]){
                    prov = eigenTemp[i][j];
                    I = i; J = j;
                }
            }
        }
        eigenTemp[I][J] = 0;

        Log::debug("Found maximum component.");

        int cumulativePop = 0; //Calcular la poblacion acumulada hasta cada nodo
        for(int i = 0; i < T.N; i++){
            for(int j = 0; j < T.vecinos[i]; j++){
                if(i == I && j == J){
                    goto outOfLoop;
                }
                cumulativePop += T.Mpesos[i][j];
            }
        }
        outOfLoop:
        Log::debug("Saved cumulative population.");

        //Add it to the chosen links matrix
        chosenLinksPop.Links++;
        chosenLinksPop.vecinos[I]++;
        chosenLinksPop.Mvecinos[I].push_back(T.Mvecinos[I][J]);
        chosenLinksPop[I].push_back(Link{static_cast<int>(T.Mpesos[I][J]), cumulativePop});
        chosenLinksPop.Total += T.Mpesos[I][J]; //Total chosen population

    }

    Log::debug("Links chosen");

    if(chosenLinksPop.Total < MUESTRA_MAX)
        Log::error("MUESTRA_MAX bigger than the number of people in the chosen links.");
    return chosenLinksPop;
}

std::vector<double> countPopulationLinks(const MobMatrix& T, const std::vector<Sparse<Link>>& chosenLinks){
    std::vector<double> results(chosenLinks.size());
    size_t temp;
    for(int l = 0; l < chosenLinks.size(); l++){
        const auto& chosenLink = chosenLinks[l];
        temp = 0;
        for(int i = 0; i < chosenLink.N; i++){
            for(int j = 0; j < chosenLink.vecinos[i]; j++){
                temp += chosenLink[i][j].Pop;
            }
        }
        results[l] = temp;
    }
    return results;
}