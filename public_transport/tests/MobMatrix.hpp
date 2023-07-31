#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

class MobMatrix{
private:

    std::string mobility_network;
    std::string pop_area;

public:
    int N = 0, Nc = 0;
    int Pob = 0;
    double Ratio = 1;
    long long Links = 0;
    std::vector<int> population;
    std::vector<double> area;
    std::vector<int> vecinos;
    std::vector<int> vecinosT;
    std::vector<std::vector<int>> Mvecinos;
    std::vector<std::vector<int>> MvecinosT;
    std::vector<std::vector<double>> Mpesos;
    std::vector<std::vector<double>> MpesosT;

    MobMatrix(const std::string& _mobility_network, const std::string& _pop_area){
        this->mobility_network = _mobility_network;
        this->pop_area = _pop_area;
        readPopArea();
        leer_vecinos();
        leer_matrices();
        calculaTraspuesta();
    }

    void insert(int i, int j, double k){
        Pob += k;
        Links++;
        vecinos[i]++;
        vecinosT[j]++;
        Mvecinos[i].push_back(j);
        Mpesos[i].push_back(k);
        MvecinosT[j].push_back(i);
        MpesosT[j].push_back(k);
    }

    void drop(int i, int j){
        auto j_pos = std::find(Mvecinos[i].begin(), Mvecinos[i].end(), j) - Mvecinos[i].begin();
        auto i_pos = std::find(MvecinosT[j].begin(), MvecinosT[j].end(), i) - MvecinosT[j].begin();

        Pob -= Mpesos[i][j_pos];

        Mvecinos[i].erase(Mvecinos[i].begin() + j_pos);
        MvecinosT[j].erase(MvecinosT[j].begin() + i_pos);

        Mpesos[i].erase(Mpesos[i].begin() + j_pos);
        MpesosT[j].erase(MpesosT[j].begin() + i_pos);

        Links--;
        vecinos[i]--;
        vecinosT[j]--;
    }

    void find_insert(int i, int j, int value = 0){
        //If it exists, add the value to the existing one
        //If not, insert new link
        auto j_it = std::find(Mvecinos[i].begin(), Mvecinos[i].end(), j);
        auto i_it = std::find(MvecinosT[j].begin(), MvecinosT[j].end(), i);

        Pob += value;

        if(j_it != Mvecinos[i].end()){
            int j_pos = j_it - Mvecinos[i].begin();
            int i_pos = i_it - MvecinosT[j].begin();

            Mpesos[i][j_pos] += value;
            MpesosT[j][i_pos] += value;
        }
        else{
            insert(i, j, value);
        }
        
    }


    std::vector<double>& operator[](int i){
        return this->Mpesos[i];
    }
    const std::vector<double>& operator[](int i) const{
        return this->Mpesos[i];
    }

private:

    void readPopArea(){
        Pob = 0;
        N = 0;

        std::ifstream inFile;
        inFile.open(this->pop_area);
        if (!inFile) {
            std::cout << "Unable to open file";
            std::exit(1); // terminate with error
        }

        //Contar la cantidad de parches
        int i, trash; double trash2;
        while(inFile >> i >> trash >> trash2){
            if(i >= N)
                N = i+1;
        }

        this->population.resize(N);
        this->area.resize(N);

        int k = 0;
        inFile.close();
        inFile.open(this->pop_area);
        while(inFile >> i >> population[i] >> area[i]){
            population[i] *= Ratio;
            Pob += population[i];
            if(i >= N)
                N = i+1;
            //std::cout << population[i] << std::endl;
            if(i != k++) {
                std::cerr << "Fichero de áreas incompleto" << std::endl;
                }
        }
        inFile.close();
        std::cout << N << " " << std::endl;
    }

    void leer_vecinos()
    {
        int i, j1, j2, trash;
        vecinos.resize(N);
        vecinosT.resize(N);
        for ( i = 0 ; i < N ; i++)
            vecinos[i] = 0;
        
        std::ifstream inFile;
        inFile.open(this->mobility_network);
        if (!inFile) {
            std::cout << "Unable to open file";
            std::exit(1); // terminate with error
        }
        Links = 0;
        while(inFile >> j1 >> j2 >> trash){
            vecinos.at(j1)++;
            vecinosT.at(j2)++;
            Links++;
        }
        inFile.close();
    }

    void leer_matrices()
    {
        Mvecinos.resize(N);
        Mpesos.resize(N);
        for(int i = 0; i < N; i++){
            Mvecinos[i].resize(vecinos[i]);
            Mpesos[i].resize(vecinos[i]);
        }
        std::ifstream inFile;
        inFile.open(this->mobility_network);
        if (!inFile){
            std::cout << "Unable to open file";
            std::exit(1); // terminate with error
        }
        int i,j,trash;
        int temp;

        //Recuento de la poblacion: los datos de pesos y poblaciones no dan una matriz Mpesos de numeros enteros
        //por lo que hay que truncarlos y recalcular las poblaciones
        Pob = 0;
        int population_temp;
        for(i = 0 ; i < N ; i++)
        {   
            temp = 0;
            for (j = 0 ; j < vecinos[i] ; j++)
            {
                inFile >> trash >> Mvecinos[i][j] >> Mpesos[i][j];
                temp += Mpesos[i][j];
            }
            population_temp = 0;
            for(j = 0; j < vecinos[i]; j++){
                Mpesos[i][j] = static_cast<int>(population[i] * Mpesos[i][j] / temp); //Mpesos/temp es equivalente a R_ij -> ahora Mpesos tiene informacion de poblacion total
                population_temp += Mpesos[i][j];
            }
            population[i] = population_temp;
            Pob += population_temp;
        }
        inFile.close();
    }

    void calculaTraspuesta()
    {  
        MvecinosT.resize(N);
        MpesosT.resize(N);
        for(int i = 0; i < N; i++){
            MvecinosT[i].resize(vecinosT[i]);
            MpesosT[i].resize(vecinosT[i]);
        }

        int i = 0, j = 0;
        int B;
        std::vector<int> AUX; AUX.resize(N);

        for (i = 0; i < N; i++){
            for (j = 0; j < vecinos[i]; j++){
                B = Mvecinos[i][j];
                MvecinosT[B][AUX[B]] = i;
                MpesosT[B][AUX[B]] = Mpesos[i][j];
                AUX[B] = AUX[B] + 1;
            }
        }
    }
};