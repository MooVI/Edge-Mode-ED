/* 
 * File:   main.cpp
 * Author: Jack Kemp
 *
 * Created on 13 October 2014, 11:54
 */

#if 1

#include<Eigen/Dense>
#include<Eigen/unsupported/Eigen/MPRealSupport>
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

using NumMethod::posmod;

typedef double mpreal;
typedef unsigned long ulong;

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

template <int N>
inline int sigma_z_j(ulong x, ulong y, ulong j, powersoftwo<N> pows) {
    return (x == y) * (2 * ((x & pows(j)) >> j) - 1);
}

template <int N>
inline int sigma_x_j(ulong x, ulong y, ulong j, powersoftwo<N> pows) {
    return ((x | pows(j)) == (y | pows(j))) and (((x^y) & pows(j)) >> j);
}

template <int N>
inline int sigma_x_j_x_m(ulong x, ulong y, ulong j, ulong m, powersoftwo<N> pows) {
    return ((x | pows(j) | pows(m)) == (y | pows(j) | pows(m))) and (((x^y) & pows(j)) >> j) and (((x^y) & pows(m)) >> m);
}

/*
 * 
 */
int main() {
    //mpreal::set_default_prec(128);
    ScatterPlotter plotter;
    const int width = 10;

    powersoftwo < width + 1 > pows2;

    typedef Eigen::Matrix<mpreal, Eigen::Dynamic, Eigen::Dynamic> Matrixww;
    typedef Eigen::Matrix<mpreal, Eigen::Dynamic, 1> Vectorw;
    typedef Eigen::Array<mpreal, Eigen::Dynamic, 1> Arrayw;
    
    const double r = 0.05;
    //const double J2 = 1.0;
    const double J = 1.0;
    const double f = r * J;
    //const double V = r*J2;
    std::ofstream outfile;
    outfile.open(std::string("diffsJ2_noV_") + std::to_string(width), std::ios::trunc);
    
    //for (double J2 = 0.5; J2 < 1.21; J2 += 0.01) {
        Matrixww H(pows2(width), pows2(width));

        //std::cout<< std::to_string(sigma_z_j(0, 0, 0, pows2)*((0+1) < width));

        for (ulong i = 0; i < pows2(width); i++) {
            for (ulong j = 0; j < pows2(width); j++) {
                for (ulong jsite = 0; jsite < width; jsite++)
                    H(i, j) += -f * sigma_x_j(i, j, jsite, pows2)
                    - J * sigma_z_j(i, j, jsite, pows2) * sigma_z_j(i, j, jsite + 1, pows2)*(((jsite + 1) < width));
                    //- J2 * sigma_z_j(i, j, jsite, pows2) * sigma_z_j(i, j, jsite + 2, pows2)*(((jsite + 2) < width));
                    //- V * sigma_x_j_x_m(i, j, jsite, jsite + 1, pows2)*((jsite + 1) < width);
            }
        }

        
        

        //std::cout << H << std::endl << std::endl;
        Eigen::SelfAdjointEigenSolver<Matrixww> es(H);
        auto eigs = es.eigenvalues();
        auto evecs = es.eigenvectors();
        const int ispec = 32;
        auto spec = evecs.col(ispec).array();
        //outfile << eigs.transpose() << std::endl;
        Arrayw overlap = Vectorw::Zero(pows2(width));
        Arrayw sigz (pows2(width));
        for (int i=0; i<sigz.size();i++){
            sigz[i] = i % 2 == 0 ? -1 : 1; 
        }
        for (int i = 0; i < sigz.size(); i++) {
                overlap[i] += (sigz*spec*evecs.col(i).array()).sum();
        }
        outfile << overlap << std::endl;
        outfile << eigs << std::endl;
        std::cout << eigs[ispec] << std::endl;
        plotter.plot2(eigs, overlap);
        plotter.wait();
    //}
    return 0;
    
}
#endif 

