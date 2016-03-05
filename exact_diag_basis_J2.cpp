/* 
 * File:   main.cpp
 * Author: Jack Kemp
 *
 * Created on 13 October 2014, 11:54
 */

#if 0

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
    return -(x == ~y) * (2 * ((y & pows(j)) >> j) - 1);
}

template <int N>
inline int sigma_x_j(ulong x, ulong y, ulong j, powersoftwo<N> pows, ulong width) {
    return (((x | pows(j)) == (y | pows(j))) and (((x^y) & pows(j)) >> j))*(-2 * ((pows(j) & x & pows(width-1)) >> (width-1)) + 1);
}

template <int N>
inline int sigma_x_j_x_m(ulong x, ulong y, ulong j, ulong m, powersoftwo<N> pows, ulong width) {
    return (j == m and x == y) ? 1 :
            (((x | pows(j) | pows(m)) == (y | pows(j) | pows(m)))
            and (((x^y) & pows(j)) >> j)
            and (((x^y) & pows(m)) >> m))
            *(-2 * ((pows(j) & x & pows(width-1)) >> (width-1)) + 1)
            *(-2 * ((pows(m) & x & pows(width-1)) >> (width-1)) + 1);
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
    const int width = 6;

    powersoftwo < width + 1 > pows2;

    typedef Eigen::Matrix<mpreal, Eigen::Dynamic, Eigen::Dynamic> Matrixww;
    typedef Eigen::Matrix<mpreal, Eigen::Dynamic, 1> Vectorw;
    typedef Eigen::Array<mpreal, Eigen::Dynamic, 1> Arrayw;

    const double r = 0.05;
    const double J2 = 0.00;
    const double J = 1.0;
    const double f = r * J;
    const double V = 0.0;
    
    /*Matrixww test (pows2(width-1), pows2(width-1));
    
    for (ulong i = 0; i < pows2(width-1); i++) {
            
                for (ulong jsite = 0; jsite < width; jsite++){
                    test (i, i) += sigma_z_j_z_m(i, i, jsite, jsite + 2, pows2)*(((jsite + 2) < width));
                }
            }
    for (ulong i = 0; i < pows2(width-1); i++) {
        std::cout<< test(i,i)<< ' ';
    }
     */
    
    std::ofstream outfile;
    //outfile.open(std::string("split_basis") + std::to_string(width), std::ios::trunc);
    std::vector<mpreal> maxoverlap, fs, eigdiffs;
    Matrixww HE(pows2(width - 1), pows2(width - 1));
    Matrixww HO(pows2(width - 1), pows2(width - 1));
    
    NumMethod::ForLoopParams<mpreal> fparams;
    fparams.numPoints = 1000; fparams.start = 0.0; fparams.end = 1.01;
    auto body = [&](mpreal J2, int i) {
        HE = Matrixww::Zero(pows2(width - 1), pows2(width - 1));
        HO = Matrixww::Zero(pows2(width - 1), pows2(width - 1));

        //std::cout<< std::to_string(sigma_z_j(0, 0, 0, pows2)*((0+1) < width));

        for (ulong i = 0; i < pows2(width - 1); i++) {
            for (ulong j = 0; j < pows2(width - 1); j++) {
                for (ulong jsite = 0; jsite < width; jsite++) {
                    HE(i, j) += -f * sigma_x_j(i, j, jsite, pows2, width)
                            - J * sigma_z_j_z_m(i, j, jsite, jsite + 1, pows2)*(((jsite + 1) < width))
                            - J2 * sigma_z_j_z_m(i, j, jsite, jsite + 2, pows2)*(((jsite + 2) < width))
                            - V * sigma_x_j_x_m(i, j, jsite, jsite + 1, pows2, width)*((jsite + 1) < width);
                    HO(i, j) += -f * sigma_x_j(~i, ~j, jsite, pows2, width)
                            - J * sigma_z_j_z_m(~i, ~j, jsite, jsite + 1, pows2)*(((jsite + 1) < width))
                            - J2 * sigma_z_j_z_m(~i, ~j, jsite, jsite + 2, pows2)*(((jsite + 2) < width))
                            - V * sigma_x_j_x_m(~i, ~j, jsite, jsite + 1, pows2, width)*((jsite + 1) < width);
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
        //const int ispec = pows2(width - 2);
        const int ispec = 0;
        auto spec = evecsE.col(ispec).array();
        //outfile << eigs.transpose() << std::endl;
        Arrayw overlap = Vectorw::Zero(pows2(width - 1));
        Arrayw sigz(pows2(width - 1));
        for (int i = 0; i < sigz.size(); i++) {
            sigz[i] = i % 2 == 0 ? 1 : -1;
        }
        for (int i = 0; i < sigz.size(); i++) {
            overlap[i] += (sigz * spec * evecsO.col(i).array()).sum();
        }
        fs.push_back(J2);
        int maxind;
        maxoverlap.push_back(overlap.abs().maxCoeff(&maxind));
        eigdiffs.push_back(fabs(eigsE[ispec]-eigsO[maxind]));
        //PlotterData pd;
        //pd.style = "l";
        //plotter.plot2(eigsE, overlap, pd);
        //std::system("pause");
        return false;
    };
    
    NumMethod::EqualSpaceFor forloop; 
    forloop.loop(body, fparams);
    //plotter.writeToFile("split_basis", eigsE, eigsO);
    //std::cout << eigsE[ispec] << ' ' << eigsO[ispec] << std::endl;
    
    PlotterData pd;
    pd.style = "l";
    pd.input = "Ising_J_ground" + std::to_string(width);
    plotter.plot2(fs, maxoverlap, pd);
    pd.input += '_eigdiff'; 
    plotter.plot2(fs, eigdiffs, pd);
    //plotter.plot2withAnalytic(fs, maxoverlap, [](mpreal x){return sqrt(1/(1.0416 + x*x*(1.04 + 3*x*x)));}, 100, pd);
    plotter.wait();
    return 0;

}
#endif 


