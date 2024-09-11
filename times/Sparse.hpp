#include "MobMatrix.hpp"

template <typename Type>
class Sparse{

public:
    int N = 0;
    long long Links = 0;
    int Total{0};

    std::vector<int> vecinos;
    std::vector<int> vecinosT;
    std::vector<std::vector<int>> Mvecinos;
    std::vector<std::vector<int>> MvecinosT;
    std::vector<std::vector<Type>> Mpesos;
    std::vector<std::vector<Type>> MpesosT;

    Sparse(){}

    Sparse(int N){
        this->N = N;
        vecinos.resize(N);
        vecinosT.resize(N);
        Mvecinos.resize(N);
        MvecinosT.resize(N);
        Mpesos.resize(N);
        MpesosT.resize(N);
    }

    Sparse(const Sparse& T){
        N = T.N;
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
        Links = T.Links;
        Total = T.Pob;
        vecinos = T.vecinos;
        vecinosT = T.vecinosT;
        Mvecinos = T.Mvecinos;
        MvecinosT = T.MvecinosT;
        Mpesos = T.Mpesos;
        MpesosT = T.MpesosT;
    }

    std::vector<Type>& operator[](int i){
        return this->Mpesos[i];
    }
    const std::vector<Type>& operator[](int i) const{
        return this->Mpesos[i];
    }

};