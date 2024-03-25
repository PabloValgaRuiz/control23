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

        Mpesos.resize(T.N);
        MpesosT.resize(T.N);
        for(int i = 0; i < T.N; i++){
            Mpesos[i].resize(T.vecinos[i]);
            MpesosT[i].resize(T.vecinosT[i]);
        }
    }

    std::vector<Type>& operator[](int i){
        return this->Mpesos[i];
    }
    const std::vector<Type>& operator[](int i) const{
        return this->Mpesos[i];
    }

};