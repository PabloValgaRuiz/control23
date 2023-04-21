#pragma once

#include <stdlib.h>
#include <vector>
#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>

class MobMatrix{
private:

    std::string mobility_network;
    std::string pop_area;

public:
    int N = 0, Nc = 0, Pob = 0;
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

    //Default constructor does nothing, every vector to length zero
    MobMatrix(){}

    //Constructor from data
    MobMatrix(const std::string& _mobility_network, const std::string& _pop_area){
        this->mobility_network = _mobility_network;
        this->pop_area = _pop_area;
        readPopArea();
        leer_vecinos();
        leer_matrices();
        calculaTraspuesta();
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
        int i, j1, j2, trash1, trash2, trash;
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
        int I,i,j,trash;
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
    
public:
    void calculaTraspuesta()
    {  
        MvecinosT.resize(N);
        MpesosT.resize(N);
        for(int i = 0; i < N; i++){
            MvecinosT[i].resize(vecinosT[i]);
            MpesosT[i].resize(vecinosT[i]);
        }
        int k = 0;
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