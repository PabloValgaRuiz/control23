#include "config.hpp"
#include "MC_Dist.hpp"
#include "benchmark.hpp"
#include <random>
#include <iostream>


MC_Dist::MC_Dist(double _lambda, double _p, const MobMatrix& M)
:lambda{_lambda}, p{_p}{
    PROFILE_FUNCTION();
    
    Est.resize(M.Pob); Org.resize(M.Pob); Des.resize(M.Pob); Desplazamiento.resize(M.Pob);
    InfectadosDes.resize(M.N); InfectadosOrg.resize(M.N); probInf.resize(M.N); n_eff.resize(M.N);
    totalDesS.resize(M.N); totalDesE.resize(M.N); totalDesA.resize(M.N); totalDesP.resize(M.N); totalDesI.resize(M.N); totalDesD.resize(M.N); totalDesR.resize(M.N);
    totalOrgS.resize(M.N); totalOrgE.resize(M.N); totalOrgA.resize(M.N); totalOrgP.resize(M.N); totalOrgI.resize(M.N); totalOrgD.resize(M.N); totalOrgR.resize(M.N);

    int k = 0;
    for(int i = 0; i < M.N; i++){
        for(int j = 0; j < M.vecinos[i]; j++){
            for(int l = 0; l < M.Mpesos[i][j]; l++){
                Org[k] = i;
                Des[k] = M.Mvecinos[i][j];
                k++;
            }
        }
    }    
}

void MC_Dist::calcInfectados(){
    PROFILE_FUNCTION();

    for(int i = 0; i < InfectadosDes.size(); i++){
        InfectadosOrg[i] = InfectadosDes[i] = 0;
        totalDesS[i] = totalDesE[i] = totalDesA[i] = totalDesP[i] = totalDesI[i] = totalDesD[i] = totalDesR[i] = 0;
        totalOrgS[i] = totalOrgE[i] = totalOrgA[i] = totalOrgP[i] = totalOrgI[i] = totalOrgD[i] = totalOrgR[i] = 0;
    }
    pobInf = totalS = totalE = totalA = totalP = totalI = totalD = totalR = 0;

    for(int k = 0; k < Est.size(); k++){
        
    #ifdef SIS
        if(Est[k] == E){
            InfectadosOrg[Org[k]]++;
            InfectadosDes[Desplazamiento[k]]++;
            pobInf++;
        }
    #endif
    #ifdef SEIR
        if(Est[k] == S){
            totalDesS[Desplazamiento[k]]++;
            totalOrgS[Org[k]]++;
            totalS++;
        }
        else if(Est[k] == E){
            totalDesE[Desplazamiento[k]]++;
            totalOrgE[Org[k]]++;
            totalE++;
        }
        else if(Est[k] == I){
            totalDesI[Desplazamiento[k]]++;
            totalOrgI[Org[k]]++;
            totalI++;
        }
        else if(Est[k] == R){
            totalDesR[Desplazamiento[k]]++;
            totalOrgR[Org[k]]++;
            totalR++;
        }
    #endif
    #ifdef SIR
        if(Est[k] == E){
            InfectadosOrg[Org[k]]++;
            InfectadosDes[Desplazamiento[k]]++;
            pobInf++;
        }
    #endif
    #ifdef SEAPIDR
        
        if(Est[k] == S){
            totalDesS[Desplazamiento[k]]++;
            totalOrgS[Org[k]]++;
            totalS++;
        }
        else if(Est[k] == E){
            totalDesE[Desplazamiento[k]]++;
            totalOrgE[Org[k]]++;
            totalE++;
        }
        else if(Est[k] == A){
            totalDesA[Desplazamiento[k]]++;
            totalOrgA[Org[k]]++;
            totalA++;
        }
        else if(Est[k] == P){
            totalDesP[Desplazamiento[k]]++;
            totalOrgP[Org[k]]++;
            totalP++;
        }
        else if(Est[k] == I){
            totalDesI[Desplazamiento[k]]++;
            totalOrgI[Org[k]]++;
            totalI++;
        }
        else if(Est[k] == D){
            totalDesD[Desplazamiento[k]]++;
            totalOrgD[Org[k]]++;
            totalD++;
        }
        else if(Est[k] == R){
            totalDesR[Desplazamiento[k]]++;
            totalOrgR[Org[k]]++;
            totalR++;
        }
    #endif

    }


}

void MC_Dist::calc_neff(){
    PROFILE_FUNCTION();

    for(int i = 0; i < n_eff.size(); i++){
        n_eff[i] = 0;
    }
    for(int k = 0; k < Est.size(); k++){
        n_eff[Desplazamiento[k]]++;
    }
}

void MC_Dist::inicializar(double _rhoInicial){
    PROFILE_FUNCTION();

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for(int i = 0; i < Est.size(); i++){
        if(dist(mt) < _rhoInicial)
            Est[i] = I;
        else Est[i] = S;
    }
}

void MC_Dist::inicializar(std::vector<double> _rhoInicial){
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for(int i = 0; i < Est.size(); i++){
        if(dist(mt) < _rhoInicial[Org[i]])
            Est[i] = E;
        else Est[i] = S;
    }
}

void MC_Dist::desplaza(){
    PROFILE_FUNCTION();

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for(int k = 0; k < Des.size(); k++){
        if (dist(mt) < p)		                //SI se desplaza
            Desplazamiento[k] = Des[k];
        else									//NO se desplaza
            Desplazamiento[k] = Org[k];
    }
}

