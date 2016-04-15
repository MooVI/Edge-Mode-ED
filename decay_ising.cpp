/* 
 * File:   main.cpp
 * Author: Jack Kemp
 *
 * Created on 13 October 2014, 11:54
 */

#if 0

#include<Eigen/Dense>
//#include<Eigen/unsupported/Eigen/MPRealSupport>
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

constexpr ulong width = 10;
constexpr ulong ignoremask = ~((1 << (width - 1)) - 1);

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
inline int sigma_x_j(ulong x, ulong y, ulong j, powersoftwo<N> pows) {
    return (
            (x == (y^pows(j)))
            ^((j == width - 1) and (
            (x | ignoremask) == ((~y) | ignoremask)
            )
            )
            )*(-2 * ((pows(j) & y & pows(width - 1)) >> (width - 1)) + 1);
}

template <int N>
inline int sigma_x_j_x_m(ulong x, ulong y, ulong j, ulong m, powersoftwo<N> pows) {
    return j == m ? (x == y ? 1 : 0) :
            (
            (x == (y^pows(j)^pows(m)))
            or
            ( (j == width - 1) and ((x | ignoremask) == ((~(y^pows(m))) | ignoremask)))
            or
            ( (m == width - 1) and ((x | ignoremask) == ((~(y^pows(j))) | ignoremask)))
            )*(-2 * ((pows(j) & y & pows(width - 1)) >> (width - 1)) + 1)
            *(-2 * ((pows(m) & y & pows(width - 1)) >> (width - 1)) + 1);
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


    powersoftwo < width + 1 > pows2;

    typedef Eigen::Matrix<mpreal, Eigen::Dynamic, Eigen::Dynamic> Matrixww;
    typedef Eigen::Matrix<mpreal, Eigen::Dynamic, 1> Vectorw;
    typedef Eigen::Array<mpreal, Eigen::Dynamic, 1> Arrayw;

    const double r = 0.05;
    const double J2 = 0.00;
    const double J = 1.0;
    const double f = 0.4;
    const double V = 0.0;

    NumMethod::RunningStats<mpreal> statoverlap, stateigdiff;

    Matrixww test(pows2(width - 1), pows2(width - 1));

    for (ulong i = 0; i < pows2(width - 1); i++) {
        for (ulong j = 0; j < pows2(width - 1); j++) {
            // for (ulong jsite = width-1; jsite < width; jsite++){
            test(i, j) = sigma_z_j_z_m(~i, ~j, 0, 1, pows2);
            // }  
        }
    }

    //   std::cout<< test<< std::endl;



    std::ofstream outfile;
    //outfile.open(std::string("split_basis") + std::to_string(width), std::ios::trunc);
    
    Matrixww HE(pows2(width - 1), pows2(width - 1));
    Matrixww HO(pows2(width - 1), pows2(width - 1));

    Arrayw sigz(pows2(width - 1));
    Arrayw sigz2(pows2(width - 1));
    for (int i = 0; i < sigz.size(); i++) {
        sigz[i] = i % 2 == 0 ? 1 : -1;
        sigz2[i] = i % 4 < 2 ? 1 : -1;
    }

    NumMethod::ForLoopParams<mpreal> fparams;
    fparams.numPoints = 4;
    fparams.start = 0.2;
    fparams.end = 0.8;
    std::vector<std::vector<mpreal>> maxoverlap (fparams.numPoints), eigdiffs (fparams.numPoints), times(fparams.numPoints);
    Arrayw moverlaps (pows2(width - 1)), meigdiff(pows2(width - 1));
    std::vector<mpreal> fs;
    auto body = [&](mpreal V, int i) {
        HE = Matrixww::Zero(pows2(width - 1), pows2(width - 1));
        HO = Matrixww::Zero(pows2(width - 1), pows2(width - 1));

        //std::cout<< std::to_string(sigma_z_j(0, 0, 0, pows2)*((0+1) < width));

        for (ulong i = 0; i < pows2(width - 1); i++) {
            for (ulong j = 0; j < pows2(width - 1); j++) {
                for (ulong jsite = 0; jsite < width; jsite++) {
                    HE(i, j) += -f * sigma_x_j(i, j, jsite, pows2)
                            - J * sigma_z_j_z_m(i, j, jsite, jsite + 1, pows2)*(((jsite + 1) < width))
                            - J2 * sigma_z_j_z_m(i, j, jsite, jsite + 2, pows2)*(((jsite + 2) < width))
                            - V * sigma_x_j_x_m(i, j, jsite, jsite + 1, pows2)*((jsite + 1) < width);
                    HO(i, j) += -f * sigma_x_j(~i, ~j, jsite, pows2)
                            - J * sigma_z_j_z_m(~i, ~j, jsite, jsite + 1, pows2)*(((jsite + 1) < width))
                            - J2 * sigma_z_j_z_m(~i, ~j, jsite, jsite + 2, pows2)*(((jsite + 2) < width))
                            - V * sigma_x_j_x_m(~i, ~j, jsite, jsite + 1, pows2)*((jsite + 1) < width);
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
        for (int ispec = 0; ispec < pows2(width - 1); ispec++) {
                //const int ispec = 0;
                auto spec = evecsE.col(ispec).array();
                //outfile << eigs.transpose() << std::endl;
                Arrayw overlap = Vectorw::Zero(pows2(width - 1));
                Arrayw accumulator(pows2(width - 1));
                for (int i = 0; i < sigz.size(); i++) {
                    //for (int j = 0; j < sigz.size(); j++)
                    //    accumulator[j] = spec[j]*(sigz[j] * evecsO.col(i)[j] + f * sigz2[j] * evecsO.col(i)[j^1]);
                    //overlap[i] = accumulator.sum();
                    overlap[i] += (spec * sigz * evecsO.col(i).array()).sum();
                }

                int maxind;
                moverlaps[ispec] = overlap.square().maxCoeff(&maxind);
                meigdiff[ispec] = eigsE[ispec] - eigsO[maxind];
                //stateigdiff.Push();
        }
        NumMethod::LogFor logfor;
        NumMethod::ForLoopParams<mpreal> logparams;
        logparams.start = 0.01; logparams.end = 1e6; logparams.numPoints = 1000;
        //for (mpreal t = 0; t < 1000000000; t += 100000) {
        auto tfunc = [&](mpreal t, int j){
           maxoverlap[i].push_back((moverlaps*meigdiff.unaryExpr([=](mpreal x) {return cos(x*t);})).mean());
           times[i].push_back(t); 
           return false;
        };
        logfor.loop(tfunc,logparams);
        plotter.writeToFile("eigdiff_dist_J2_"+std::to_string(J2), meigdiff);
            //eigdiffs[i].push_back(stateigdiff.Mean());
            statoverlap.Clear();
            //stateigdiff.Clear();
        fs.push_back(V);
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

    PlotterData pd, pd2;
    pd.style = "l";
    pd.input = "V_Ising_test_f_" +std::to_string(f)+"_V_"+ std::to_string(fs[0])+ "_decay_" + std::to_string(width);
    plotter.plot2(times[0], maxoverlap[0], pd);
    //plotter.plot2withAnalytic(times[0], maxoverlap[0], [](mpreal x){return exp(-x*x*1.4e-15);}, 100, pd);
    pd2.input = "V_Ising_test_f_" +std::to_string(f)+"_V_" + std::to_string(fs[1])+ "_decay_" + std::to_string(width);
    plotter.plot2(times[1], maxoverlap[1], pd2);
    //plotter.plot2withAnalytic(times[1], maxoverlap[1], [](mpreal x){return (0.538-0.36)*exp(-x*x*2.8e-6)+0.36;}, 1000, pd2);
    plotter.wait();
    return 0;

}
#endif 




