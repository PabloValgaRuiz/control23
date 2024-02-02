#include "config.hpp"
#include "MC_DistDiaNocheF.hpp"
#include "benchmark.hpp"
#include <chrono>
#include <random>
#include <math.h>


MC_DistDiaNocheF::MC_DistDiaNocheF(double _lambda, double _p, const MobMatrix& M)
: MC_DistDiaNocheF::MC_Dist{_lambda,_p,M}
{
    PROFILE_FUNCTION();
    
    calculaLambda0(M);
    probInfDes.resize(M.N); probInfOrg.resize(M.N); fvector.resize(M.N);
}

void MC_DistDiaNocheF::update(const MobMatrix& T){
    PROFILE_FUNCTION();

    this->desplaza();
    this->calc_neff();
    generafvector(T);
    this->calculateZs(T);
    this->calcInfectados();
}

void MC_DistDiaNocheF::iteracion(const MobMatrix& T){

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    this->update(T);

{PROFILE_SCOPE("Iteration1");
    for(int i = 0; i < T.N; i++){
        if(n_eff[i] != 0){
            #ifdef SIS
                probInfDes[i] = 1 - lambda * static_cast<double>(InfectadosDes[i])/n_eff[i];
            #endif
            #ifdef SEIR
                probInfDes[i] = 1 - lambda * static_cast<double>(totalDesI[i])/n_eff[i];
            #endif
            #ifdef SIR
                probInfDes[i] = 1 - lambda * static_cast<double>(InfectadosDes[i])/n_eff[i];
            #endif
            #ifdef SEAPIDR
                probInfDes[i] = 1 - static_cast<double>(lambdaA * totalDesA[i] + lambdaP * totalDesP[i] + lambdaI * totalDesI[i])/n_eff[i];
            #endif
        }
        else probInfDes[i] = 1;
        if(T.population[i] != 0){
            #ifdef SIS
                probInfOrg[i] = 1 - lambda * static_cast<double>(InfectadosOrg[i])/T.population[i];
            #endif
            #ifdef SEIR
                probInfOrg[i] = 1 - lambda * static_cast<double>(totalOrgI[i])/T.population[i];
            #endif
            #ifdef SIR
                probInfOrg[i] = 1 - lambda * static_cast<double>(InfectadosOrg[i])/T.population[i];
            #endif
            #ifdef SEAPIDR
                probInfOrg[i] = 1 - static_cast<double>(lambdaA * totalOrgA[i] + lambdaP * totalOrgP[i] + lambdaI * totalOrgI[i])/T.population[i];
            #endif
        }
        else probInfOrg[i] = 1;
    }
}
    //Si zD * fi = 8.2, hay que hacer 8 contactos, y luego una probabilidad de 0.2 de tener un 9no contacto
    int contactosD;
    int contactosN;
    double probContactoD;
    double probContactoN;
    double probInfIndiv;
    
{PROFILE_SCOPE("Iteration2");
    //Contagios por el dia
    for(int k = 0; k < Est.size(); k++){
        switch(Est[k]){
            case S:
                contactosD = static_cast<int>(zD * fvector[Desplazamiento[k]]);
                contactosN = static_cast<int>(zN * sigma);
                probContactoD = zD * fvector[Desplazamiento[k]] - static_cast<double>(contactosD);
                probContactoN = zN * sigma - static_cast<double>(contactosN);

                if(dist(mt) < probContactoD)
                    probInfIndiv = 1 - pow(probInfDes[Desplazamiento[k]], contactosD + 1); //9 contactos
                else probInfIndiv = 1 - pow(probInfDes[Desplazamiento[k]], contactosD);     //8 contactos

                if(dist(mt) < probInfIndiv){
                    Est[k] = E;
                }
                else{   //Contagios por la noche
                    if(dist(mt) < probContactoN)
                        probInfIndiv = 1 - pow(probInfOrg[Org[k]], contactosN + 1);
                    else probInfIndiv = 1 - pow(probInfOrg[Org[k]], contactosN);
                    if(dist(mt) < probInfIndiv){
                        Est[k] = E;
                    }
                }

                break;
            
            case E:
#ifdef SIS
                if(dist(mt) < mu)
                    Est[k] = S;
#endif
#ifdef SEIR
                if(dist(mt) < nu)
                    Est[k] = I;
                break;
            case I:
                if(dist(mt) < mu)
                    Est[k] = R;
                break;
#endif
#ifdef SIR
                if(dist(mt) < mu)
                    Est[k] = R;
#endif
#ifdef SEAPIDR
                double temp = dist(mt);
                if(temp < (1-x)*nu)
                    Est[k] = A;
                else if(temp < nu)
                    Est[k] = P;
#endif
                break;
#ifdef SEAPIDR
            case A:
                if(dist(mt) < muA)
                    Est[k] = R;
                break;
            
            case P:
                if(dist(mt) < alpha)
                    Est[k] = I;
                break;
            
            case I:
                temp = dist(mt);
                if(temp < delta)
                    Est[k] = D;
                else if(temp < delta + muI)
                    Est[k] = R;
                break;
            case D:
                if(dist(mt) < gamma)
                    Est[k] = R;
#endif
        }
    }
}
    
}


void MC_DistDiaNocheF::calculaLambda0(const MobMatrix& T){
    PROFILE_FUNCTION();

    double zD = 0, zN = 0;

    std::vector<double> n_eff_p1; n_eff_p1.resize(T.N);
    for(int i = 0; i < T.N; i++)
        for(int j = 0; j < T.vecinos[i]; j++)
            n_eff_p1[T.Mvecinos[i][j]] += T.Mpesos[i][j];

    double temp1 = 0;
    for(int i = 0; i < T.N; i++){
        temp1 += n_eff_p1[i] * f(i, 1, T);
    }
    zD = T.Pob * kD / temp1;

    temp1 = 0;
    for(int i = 0; i < T.N; i++){
        temp1 += T.population[i] * sigma;
    }
    zN = T.Pob * kN / temp1;
    
    //Calcular el maximo de zD*fi + zN*sigma -> maximo de fi si las sigmas son iguales
    double F = 0, SIGMA = 0;
    double temp2 = 0; temp1 = 0;
    for(int i = 0; i < T.N; i++){
        temp1 = f(i, 0, T);
        temp2 = sigma; //sigma[i] si fuese un vector
        if(zD*F + zN*SIGMA < zD*temp1 + zN*temp2) {F = temp1;  SIGMA = temp2;}
    }
    #ifdef SIS
        lambda0 = mu / (zD * F + zN * SIGMA);
    #endif

    #ifdef SEIR
        lambda0 = mu / (zD * F + zN * SIGMA);
    #endif

    #ifdef SIR
        lambda0 = mu / (zD * F + zN * SIGMA);
    #endif

    #ifdef SEAPIDR
        lambda0 = 1 / ( (zD * F + zN * SIGMA) * ( (1 - x)/(2*muA) + x/(delta + muI) + x/alpha ) );
    #endif
    
}

double MC_DistDiaNocheF::f(int i, double _p, const MobMatrix& T) const{

    double result = 0;  //Calcular la poblacion efectiva neff
    for(int j = 0; j < T.vecinosT[i]; j++){
        result += T.MpesosT[i][j] * _p;
    }
    for(int j = 0; j < T.vecinos[i]; j++){
        result += T.Mpesos[i][j] * (1 - _p);
    }
    if(T.area[i] != 0)
        return result/T.area[i];
    else return 0;
}

double MC_DistDiaNocheF::f(int i, const MobMatrix& T) const{
    return MC_DistDiaNocheF::f(i, this->p, T);
}

void MC_DistDiaNocheF::generafvector(const MobMatrix& T){
    PROFILE_FUNCTION();

    for(int i = 0; i < T.N; i++){
        fvector[i] = f(i, T);
    }
}

void MC_DistDiaNocheF::calculateZs(const MobMatrix&T){
    PROFILE_FUNCTION();

    //Calcular las z en funcion de p=1

    zD = zN = 0;
    double prov1 = 0, prov2 = 0;

    std::vector<double> n_eff_p1; n_eff_p1.resize(T.N);
    for(int i = 0; i < T.N; i++)
        for(int j = 0; j < T.vecinos[i]; j++)
            n_eff_p1[T.Mvecinos[i][j]] += T.Mpesos[i][j];


    for(int i = 0; i < T.N; i++){
        prov1 += n_eff_p1[i] * f(i, 1, T);
        prov2 += T.population[i] * sigma;
    }
    zD = T.Pob * kD / prov1;
    zN = T.Pob * kN / prov2;
}