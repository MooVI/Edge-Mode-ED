#if 0

#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<Eigen/MPRealSupport>
#include<NumericalMethods/NumericalMethods/Random.h>
#include<NumericalMethods/NumericalMethods/Statistics.h>
#include<NumericalMethods/NumericalMethods/SamplingForLoops.h>
#include<NumericalMethods/NumericalMethods/Mod.h>
#include<NumericalMethods/NumericalMethods/Useful.h>
#include<Plotting/Plotter/Plotter.h>
#include<iostream>
#include<mpreal.h>
#include<chrono>
#include<functional>
#include<string>

//typedef mpfr::mpreal mpreal;
typedef double mpreal;
typedef unsigned long ulong;

template<class C, class T>
auto contains(const C& v, const T& x)
-> decltype(end(v), true) {
    return end(v) != std::find(begin(v), end(v), x);
}

template <int N>
class powersoftwo {
    int table [N + 1];
public:

    int operator()(int i) {
        return table [i];
    }

    const int * getTable() {
        return table;
    };

    powersoftwo() {
        for (int i = 0, pow = 1; i < N + 1; i++) {
            table [i] = pow;
            pow = pow << 1;
        }
    }
};

int main(int argc, char** argv) {

    constexpr int maxwidth = 16;
    powersoftwo < maxwidth + 1 > pows2;

    typedef Eigen::Matrix<mpreal, Eigen::Dynamic, Eigen::Dynamic> Matrixww;
    typedef Eigen::SparseMatrix<mpreal> SparseMatrixww;
    typedef Eigen::Triplet<mpreal> Triplet;
    typedef Eigen::Matrix<mpreal, Eigen::Dynamic, 1> Vectorw;
    typedef Eigen::Array<mpreal, Eigen::Dynamic, 1> Arrayw;
    typedef Eigen::Array<mpreal, Eigen::Dynamic, Eigen::Dynamic> Arrayww;

    constexpr int width = 3;

    
    std::vector<mpreal> zpos = {2};
    std::vector<mpreal> xpos = {0,1,2};
    
    std::vector<Triplet> result(pows2(width - 1));
    SparseMatrixww mat(pows2(width - 1), pows2(width - 1));

    bool z = contains(zpos, 0);
    bool x = contains(xpos, 0);
    if (x and z) {
        result [0] = Triplet(0, 1, -1);
        result [1] = Triplet(1, 0, 1);
    } else if (z) {
        result [0] = Triplet(0, 0, -1);
        result [1] = Triplet(1, 1, 1);
    } else if (x) {
        result [0] = Triplet(0, 1, 1);
        result [1] = Triplet(1, 0, 1);
    } else {
        result [0] = Triplet(0, 0, 1);
        result [1] = Triplet(1, 1, 1);
    }
    for (int i = 1; i < width - 1; i++) {
        bool z = contains(zpos, i);
        bool x = contains(xpos, i);
        if (x and z) {
           for (int j = 0; j < pows2(i); j++) {
                result[j] = Triplet(result[j].row() + pows2(i), result[j].col(), result[j].value());
                result[j + pows2(i)] = Triplet(result[j].row()-pows2(i), result[j].col() + pows2(i), -result[j].value());
            } 
        }
        else if (z){
            for (int j = 0; j < pows2(i); j++) {
                result[j] = Triplet(result[j].row(), result[j].col(), -result[j].value());
                result[j + pows2(i)] = Triplet(result[j].row() + pows2(i), result[j].col() + pows2(i), -result[j].value());
            }
        }
        else if (x) {
            for (int j = 0; j < pows2(i); j++) {
                result[j] = Triplet(result[j].row() + pows2(i), result[j].col(), result[j].value());
                result[j + pows2(i)] = Triplet(result[j].row()-pows2(i), result[j].col() + pows2(i), result[j].value());
            }
        } else {
            for (int j = 0; j < pows2(i); j++) {
                result[j + pows2(i)] = Triplet(result[j].row() + pows2(i), result[j].col() + pows2(i), result[j].value());
            }
        }
    }
    z = contains(zpos, width-1);
    x = contains(xpos, width-1);
    if (x and z) {
        for (int j = 0; j < result.size(); j++) {
            result[j] = Triplet(pows2(width - 1) - result[j].row()-1, result[j].col(), -result[j].value());
        }
    } else if (z) {
        for (int j = 0; j < result.size(); j++) {
            result[j] = Triplet(result[j].row(), result[j].col(), -result[j].value());
        }
    } else if (x) {
        for (int j = 0; j < result.size(); j++) {
            result[j] = Triplet(pows2(width - 1) - result[j].row()-1, result[j].col(), -result[j].value());
        }
    } 
    
    mat.setFromTriplets(result.begin(), result.end());

    std::cout << mat;

    return 0;
}

#endif