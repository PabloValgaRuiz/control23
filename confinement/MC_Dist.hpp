#pragma once
#include "config.hpp"
#include "MobMatrix.hpp"
#include <vector>

enum Estado : unsigned char {S, E, A, P, I, D, R};
class MC_Dist{          //Montecarlo para individuos Distinguibles -> base

public:

    int pobInf;
    int totalS, totalE, totalA, totalP, totalI, totalD, totalR;

    MC_Dist(double _lambda, double _p, const MobMatrix& M);

    virtual void calculaLambda0(const MobMatrix& M) = 0;
    virtual void inicializar(double _rhoInicial);
    virtual void inicializar(std::vector<double> _rhoInicial);
    virtual void iteracion(const MobMatrix& T) = 0;
    // virtual void contarInfectados() = 0;

    virtual void calcInfectados();
    virtual void calc_neff();
    virtual void desplaza();

    const std::vector<Estado>& getEst() const {return Est;}
    std::vector<Estado>& getEst() {return Est;} //Modificable
    const std::vector<int>& getOrg() const {return Org;}
    const std::vector<int>& getDes() const {return Des;}
    const std::vector<int>& getDesplazamiento() const {return Desplazamiento;}
    const std::vector<int>& getInfectadosDes() const {return InfectadosDes;}
    const std::vector<int>& getInfectadosOrg() const {return InfectadosOrg;}
    const std::vector<bool>& getIsDisplaced() const {return isDisplaced;}
    const std::vector<double>& getProbInf() const {return probInf;}
    double getLambda0() const {return lambda0;}

    void setLambda(double _lambda){
    #ifdef SIS
        lambda = _lambda * lambda0;
    #endif

    #ifdef SEIR
        lambda = _lambda * lambda0;
    #endif

    #ifdef SIR
        lambda = _lambda * lambda0;
    #endif

    #ifdef SEAPIDR
        lambdaP = _lambda * lambda0;
        lambdaI = _lambda * lambda0;
        lambdaA = _lambda * lambda0 / 2;
    #endif
    }

protected:
    double lambda;  //Probabilidad de infeccion por contacto
    double lambda0, mu = 0.2;
    double lambdaP = 0.07, lambdaI = 0.07, lambdaA = 0.035;
    double nu = 0.50; //   1.0/2.6
    double alpha = 1.0/2.6, delta = 1.0/3.0, gamma = 1.0/14.0, muI = 1.0/4.2, muA = 1.0/6.8, x = 0.35;
    double p;		//Probabilidad de desplazamiento
    std::vector<int> Org, Des, Desplazamiento, InfectadosDes, InfectadosOrg, n_eff;
    std::vector<bool> isDisplaced;
    std::vector<int> totalDesS, totalDesE, totalDesA, totalDesP, totalDesI, totalDesD, totalDesR;
    std::vector<int> totalOrgS, totalOrgE, totalOrgA, totalOrgP, totalOrgI, totalOrgD, totalOrgR;
    std::vector<Estado> Est;
    std::vector<double> probInf;
    
};