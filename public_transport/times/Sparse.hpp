#pragma once

#include <vector>

template <typename Type>
class Sparse{
public:
    int N = 0;
    int M = 0;
    unsigned long long Links = 0;
    int Total{};

    std::vector<int> vecinos;
    std::vector<int> vecinosT;
    std::vector<std::vector<int>> Mvecinos;
    std::vector<std::vector<int>> MvecinosT;
    std::vector<std::vector<Type>> Mpesos;
    std::vector<std::vector<Type>> MpesosT;

    Sparse(){}

    Sparse(int _N, int _M){
        N = _N;
        M = _M;
        vecinos.resize(N);
        vecinosT.resize(M);
        Mvecinos.resize(N);
        MvecinosT.resize(M);
        Mpesos.resize(N);
        MpesosT.resize(M);
    }
    Sparse(int _N) : Sparse(_N, _N){}

    Sparse(const Sparse& T){
        N = T.N;
        M = T.M;
        Links = T.Links;
        Total = T.Total;
        vecinos = T.vecinos;
        vecinosT = T.vecinosT;
        Mvecinos = T.Mvecinos;
        MvecinosT = T.MvecinosT;
        Mpesos = T.Mpesos;
        MpesosT = T.MpesosT;
    }

    Sparse(const MobMatrix& T){
        N = T.N;
        M = T.N;
        Links = T.Links;
        Total = T.Pob;
        vecinos = T.vecinos;
        vecinosT = T.vecinosT;
        Mvecinos = T.Mvecinos;
        MvecinosT = T.MvecinosT;

        Mpesos.resize(T.N);
        MpesosT.resize(T.N);
        for(int i = 0; i < T.N; i++){
            Mpesos[i].resize(T.vecinos[i]);
            MpesosT[i].resize(T.vecinosT[i]);
        }
    }

    void insert(int i, int j, Type k){
        Links++;
        vecinos[i]++;
        vecinosT[j]++;
        Mvecinos[i].push_back(j);
        Mpesos[i].push_back(k);
        MvecinosT[j].push_back(i);
        MpesosT[j].push_back(k);
    }

    std::vector<Type>& operator[](int i){
        return this->Mpesos[i];
    }
    const std::vector<Type>& operator[](int i) const{
        return this->Mpesos[i];
    }

};