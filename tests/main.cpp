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

const std::string path = "../";
struct Link{
    //Number of population that the link contains
    int Pop;
    //Number of population that have come before this link among the ones selected
    int cumulativePop;
};

struct Result{
    double mean;
    double mean2;

    int population_link{0};
};

Sparse<double> readEigen(const MobMatrix& T, const std::string& name);
Sparse<Link> chooseLinks(int NlinksChosen, const MobMatrix& T, const Sparse<double>& eigenVector);
void countPopulationLinks(const MobMatrix& T, const std::vector<Sparse<Link>>& chosenLinks, std::vector<std::vector<Result>>& heatMap);
void iterations(const MobMatrix& T, const std::vector<Sparse<Link>>& chosenLinks, std::vector<std::vector<Result>>& heatMap);
int contarInfectadosChosen(const MobMatrix& T, const Sparse<Link>& chosenLinks, MC_DistDiaNocheF& montecarlo, int MUESTRA, std::mt19937& generator);
std::vector<int> fisherYatesShuffle(int k, std::vector<int> range, std::mt19937& generator);
std::vector<int> reservoirSampling(int k, int n, std::mt19937& generator);

#define MUESTRA_MAX 1000 //50000

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

std::string name = "ma";

static const double beta = 3.0 * cityBeta.at(name);
static const double p = 1.0;
static const int nPasos = 300;//300
static const int nIterations = 1;

std::mutex resultsMutex;


int main(int argc, char* argv[]){

    Instrumentor::Get().BeginSession("Session Name");
{
    InstrumentationTimer timer("Program");

    ThreadPool pool{24};

    std::string output = path + "out/" + name + "_beta_3,0.txt";

    MobMatrix T{path + "cities3/" + name + "/mobnetwork.txt", path + "cities3/" + name + "/Poparea.txt"};
    std::cout << T.Pob << std::endl;
    auto eigenVector = readEigen(T, name);

    Log::debug("EigenVector read.");

    size_t sizeLinks = 33; //33

    std::vector<std::vector<Result>> heatMap(sizeLinks);
    for(auto& i : heatMap)  i.resize(nPasos);

    std::vector<Sparse<Link>> vectorChosenLinks(sizeLinks);

    //________________________________CHOOSING LINKS_________________________________

    std::vector<std::future<void>> futures;
    for(size_t i = 0; i < sizeLinks; ++i){
        futures.push_back(std::move(pool.enqueue([&, i]{
            size_t NlcTemp = T.Links * (i+1) / (3 * sizeLinks); //Care to not choose 0
            //Choose the Nlc highest component links in the eigenvector
            vectorChosenLinks[i] = chooseLinks(NlcTemp, T, eigenVector);
        })));
	}
	for(auto& future : futures){
        future.wait();
    }
    futures.clear();

    Log::info("Links chosen.");

    //_______________________________COUNTING PEOPLE___________________________________

    countPopulationLinks(T, vectorChosenLinks, heatMap);

    Log::info("Population counted.");
    //__________________________________ITERATING_____________________________________
        
    for(int l = 0; l < nIterations; l++){
        futures.push_back(std::move(pool.enqueue([&, l]{

            iterations(T, vectorChosenLinks, heatMap);

        })));
    }
    for(auto& future : futures)
        future.wait();


    std::ofstream f(output);
	for(int i = 0; i < heatMap.size(); i++)//iteracion sobre links
		for(int j = 0; j < heatMap[i].size(); j++)//iteracion sobre tests
                //Cantidad de links: Copia de mas arriba, al elegir los links
			f << heatMap[i][j].population_link << "\t" << ((j+1) * 10000) * nPasos / heatMap[i].size() << "\t" << heatMap[i][j].mean/nIterations << "\t" <<
            2 * sqrt((heatMap[i][j].mean2 - (heatMap[i][j].mean * heatMap[i][j].mean / (nIterations * nPasos)))/((nIterations * nPasos) * ((nIterations * nPasos)-1))) << "\n";
	f.close();


    Instrumentor::Get().EndSession();
}
    return 0;
}

void countPopulationLinks(const MobMatrix& T, const std::vector<Sparse<Link>>& chosenLinks, std::vector<std::vector<Result>>& heatMap){
    size_t temp;
    for(int l = 0; l < chosenLinks.size(); l++){
        const auto& chosenLink = chosenLinks[l];
        temp = 0;
        for(int i = 0; i < chosenLink.N; i++){
            for(int j = 0; j < chosenLink.vecinos[i]; j++){
                temp += chosenLink[i][j].Pop;
            }
        }
        for(auto& point : heatMap[l]){
            point.population_link = temp;
        }
    }
}

void iterations(const MobMatrix& T, const std::vector<Sparse<Link>>& chosenLinks, std::vector<std::vector<Result>>& heatMap){
    
    std::mt19937 mt(std::random_device{}());

    MC_DistDiaNocheF montecarlo{0, p, T};
    montecarlo.setLambda(beta);
    montecarlo.inicializar(0.0004);

    for(int t = 0; t < nPasos; t++){
        montecarlo.iteracion(T);

        
        
        for(int l = 0; l < chosenLinks.size(); l++){

            int statesIndex = 0;

            auto chosenLink = chosenLinks[l]; //The last group of links (all of them)
            std::vector<int> MCLabels(chosenLink.Total);

            for(int i = 0; i < chosenLink.N; i++){
                for(int j = 0; j < chosenLink.vecinos[i]; j++){
                    //Find the first person in the link in the Monte Carlo agent array
                    int MCIndex = chosenLink[i][j].cumulativePop;

                    //Copy the adress of every person included in the reserve
                    for(int k = 0; k < chosenLink[i][j].Pop; k++)
                        MCLabels[statesIndex + k] = MCIndex + k;

                    //Add the number of included people
                    statesIndex += chosenLink[i][j].Pop;
                }
            }

            
            int PobInf = contarInfectadosChosen(T, chosenLinks[l], montecarlo, MUESTRA_MAX, mt);
            
            MCLabels = fisherYatesShuffle(MUESTRA_MAX, MCLabels, mt);

            //Take pieces from the sample
            for(int j = 0; j < heatMap[l].size(); j++){//iteracion sobre tests

                int MUESTRA = ((j+1) * MUESTRA_MAX) / heatMap[l].size();

                PobInf = 0;
                for(auto label : MCLabels){
                    if(montecarlo.getEst()[label] == E || montecarlo.getEst()[label] == I){
                        PobInf++;
                        //montecarlo.getEst()[label] = R;
                    }
                }

                std::lock_guard<std::mutex> lock(resultsMutex);
                heatMap[l][j].mean += static_cast<double>(PobInf);
                heatMap[l][j].mean2 += static_cast<double>(PobInf) * static_cast<double>(PobInf);
            }
        }
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
            //montecarlo.getEst()[label] = R;
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
        chosenLinksPop.Total += static_cast<int>(T.Mpesos[I][J]); //Total chosen population
    }

    Log::debug("Links chosen");

    if(chosenLinksPop.Total < MUESTRA_MAX)
        Log::error("MUESTRA_MAX bigger than the number of people in the chosen links.");
    return chosenLinksPop;
}