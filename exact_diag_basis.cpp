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


/* The basis is |downs> + |ups> -> MSB not set, |downs> - |highs> -> MSB set */


template <int N>
inline int sigma_z_j(ulong x, ulong y, ulong j, powersoftwo<N> pows, ulong widthmask) {
    return -((x^widthmask) == y) * (2 * ((x & pows(j)) >> j) - 1);
}

template <int N>
inline int sigma_x_j(ulong x, ulong y, ulong j, powersoftwo<N> pows, ulong widthmask) {
    return (((x | pows(j)) == (y | pows(j))) and (((x^y) & pows(j)) >> j))*(-2*((j & x & widthmask) >> widthmask)+1);
}

template <int N>
inline int sigma_x_j_x_m(ulong x, ulong y, ulong j, ulong m, powersoftwo<N> pows, ulong widthmask) {
    return (j==m and x==y) ? 1 : 
            (((x | pows(j) | pows(m)) == (y | pows(j) | pows(m))) 
            and (((x^y) & pows(j)) >> j) 
            and (((x^y) & pows(m)) >> m))
            *(-2*((j & x & widthmask) >> widthmask)+1)
            *(-2*((m & x & widthmask) >> widthmask)+1);
}

template <int N>
inline int sigma_z_j_z_m(ulong x, ulong y, ulong j, ulong m, powersoftwo<N> pows) {
    return (x == y) * (2 * ((x & pows(j)) >> j) - 1)*(2 * ((x & pows(m)) >> m) - 1);
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
    outfile.open(std::string("split_basis") + std::to_string(width), std::ios::trunc);
    
    //for (double J2 = 0.5; J2 < 1.21; J2 += 0.01) {
        Matrixww HE(pows2(width-1), pows2(width-1));
        Matrixww HO(pows2(width-1), pows2(Width-1));

        //std::cout<< std::to_string(sigma_z_j(0, 0, 0, pows2)*((0+1) < width));

        for (ulong i = 0; i < pows2(width-1); i++) {
            for (ulong j = 0; j < pows2(width-1); j++) {
                for (ulong jsite = 0; jsite < width; jsite++){
                    HE(i, j) += -f * sigma_x_j(i, j, jsite, pows2, pows2(width-1))
                    - J * sigma_z_j_z_m(i, j, jsite, jsite + 1, pows2)*(((jsite + 1) < width));
                    //- J2 * sigma_z_j(i, j, jsite, pows2) * sigma_z_j(i, j, jsite + 2, pows2)*(((jsite + 2) < width));
                    //- V * sigma_x_j_x_m(i, j, jsite, jsite + 1, pows2)*((jsite + 1) < width);
                    HO(i, j) += -f * sigma_x_j(~i, ~j, jsite, pows2, pows2(width-1))
                    - J * sigma_z_j_z_m(~i, ~j, jsite, jsite + 1, pows2)*(((jsite + 1) < width));
                    //- J2 * sigma_z_j(i, j, jsite, pows2) * sigma_z_j(i, j, jsite + 2, pows2)*(((jsite + 2) < width));
                    //- V * sigma_x_j_x_m(i, j, jsite, jsite + 1, pows2)*((jsite + 1) < width
                }    
            }
        }

        
        

        //std::cout << H << std::endl << std::endl;
        Eigen::SelfAdjointEigenSolver<Matrixww> esE(HE);
        Eigen::SelfAdjointEigenSolver<Matrixww> esO(HO);
        auto eigsE = esE.eigenvalues();
        auto eigsO = esO.eigenvalues();
        auto evecsE = esE.eigenvectors();
        auto evecsO = esO.eigenvectors();
        const int ispec = 32;
        auto spec = evecsE.col(ispec).array();
        //outfile << eigs.transpose() << std::endl;
        Arrayw overlap = Vectorw::Zero(pows2(width-1));
        Arrayw sigz (pows2(width-1));
        for (int i=0; i<sigz.size();i++){
            sigz[i] = i % 2 == 0 ? 1 : -1; 
        }
        for (int i = 0; i < sigz.size(); i++) {
                overlap[i] += (sigz*spec*evecsO.col(i).array()).sum();
        }
        outfile << overlap << std::endl;
        outfile << eigsE << std::endl;
        outfile << eigsO << std::endl;
        std::cout << eigsE[ispec] << std::endl;
        plotter.plot2(eigs, overlap);
        plotter.wait();
    //}
    return 0;
    
}
#endif 

