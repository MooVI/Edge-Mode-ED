/* 
 * File:   main.cpp
 * Author: Jack Kemp
 *
 * Created on 13 October 2014, 11:54
 */

#if 1

#define EIGEN_USE_MKL_ALL 1

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

#

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
inline int sigma_x_j(ulong x, ulong y, ulong j, powersoftwo<N> pows, const int width, const int ignoremask) {
    return (
            (x == (y^pows(j)))
            ^((j == width - 1) and (
            (x | ignoremask) == ((~y) | ignoremask)
            )
            )
            )*(-2 * ((pows(j) & y & pows(width - 1)) >> (width - 1)) + 1);
}

template <int N>
inline int sigma_x_j_x_m(ulong x, ulong y, ulong j, ulong m, powersoftwo<N> pows, const int width, const int ignoremask) {
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

template <int N>
inline int sigma_z_p_x_j_z_m(ulong x, ulong y, ulong p, ulong j, ulong m, powersoftwo<N> pows, const int width, const int ignoremask) {
    return (
            (x == (y^pows(j)))
            ^((j == width - 1) and (
            (x | ignoremask) == ((~y) | ignoremask)
            )
            )
            )*(-2 * ((pows(j) & y & pows(width - 1)) >> (width - 1)) + 1)
            * (2 * ((x & pows(p)) >> p) - 1)
            *(2 * ((x & pows(m)) >> m) - 1);
}


template<typename T>
std::string to_string(T value){
    return std::to_string(value);
}

template<>
std::string to_string(mpfr::mpreal value){
    return std::to_string(value.toDouble());
}


/*
 * 
 */
int main(int argc, char** argv) {
    //mpreal::set_default_prec(128);
    ScatterPlotter plotter;

    constexpr int maxwidth = 16;
    powersoftwo < maxwidth + 1 > pows2;
    
    const bool WRITE_ENERGIES = false;
    const bool WRITE_OVERLAPS = false;
    const bool WRITE_BULK = false;
    const bool WRITE_MEANS = false;
    const bool WRITE_VARS = true;
    const bool WRITE_MAX_OVERLAPS = true;
    const bool WRITE_PAIRED_EDIFFS = false;
    const bool WRITE_ALL_DECAY = false;
    const bool WRITE_PAIRED_DECAY = false;

    const bool CMD_LINE_PARAMS = true;


    typedef Eigen::Matrix<mpreal, Eigen::Dynamic, Eigen::Dynamic> Matrixww;
    typedef Eigen::Matrix<mpreal, Eigen::Dynamic, 1> Vectorw;
    typedef Eigen::Array<mpreal, Eigen::Dynamic, 1> Arrayw;
    typedef Eigen::Array<mpreal, Eigen::Dynamic, Eigen::Dynamic> Arrayww;

    const mpreal r = 0.05;
    const mpreal J2 = 0.00;
    const mpreal J3 = 0.0;
    const mpreal J4 = 0.0;
    const mpreal J = 1.0;
    const mpreal f = 0.0;
    const mpreal V = 0.05;
    const mpreal Js [] = {J3, J4};

    NumMethod::RunningStats<mpreal> statoverlap, stateigdiff;


    const int begin = 8;
    const int end = 13;

    std::string hashlabel = "";
    if (CMD_LINE_PARAMS and argc > 4)
      hashlabel = std::string(argv[4]);

    NumMethod::ForLoopParams<mpreal> fparams;
    NumMethod::GetXFor<mpreal> recordx;
    NumMethod::EqualSpaceFor couplingsfor;
    fparams.start = 0.0;
    fparams.end = 2.0;
    fparams.numPoints = 10;

    if (CMD_LINE_PARAMS)
      fparams = NumMethod::get_for_from_cmd<mpreal>(argv);
    


    std::vector<mpfr::mpreal> maxoverlap, eigdiffs, varoverlaps, vareigdiffs;
    std::vector<mpreal> fs;
    couplingsfor.loop(recordx, fparams);
    fs = recordx.get_x();
    recordx.clear();

    NumMethod::LogFor tfor;
    NumMethod::ForLoopParams<mpreal> tparams;
    std::vector<mpreal> ts;
    if (WRITE_ALL_DECAY | WRITE_PAIRED_DECAY) {
        tparams.start = 0.01;
        tparams.end = 1e14;
        tparams.numPoints = 2000;
        tfor.loop(recordx, tparams);
        ts = recordx.get_x();
        recordx.clear();
    }

    auto body = [&](int width, int i) {



        Arrayw moverlaps(pows2(width - 1)), meigdiff(pows2(width - 1));

        int ignoremask = ~((1 << (width - 1)) - 1);

        Matrixww HE(pows2(width - 1), pows2(width - 1));
        Matrixww HO(pows2(width - 1), pows2(width - 1));

        Arrayw sigz(pows2(width - 1));
        Arrayw sigz2(pows2(width - 1));
        for (int i = 0; i < sigz.size(); i++) {
            sigz[i] = i % 2 == 0 ? 1 : -1;
            sigz2[i] = i % 4 < 2 ? 1 : -1;
        }

        auto couplingsbody = [&](mpreal J2, int j) {
            HE = Matrixww::Zero(pows2(width - 1), pows2(width - 1));
            HO = Matrixww::Zero(pows2(width - 1), pows2(width - 1));

            //std::cout<< std::to_string(sigma_z_j(0, 0, 0, pows2)*((0+1) < width));

            for (ulong i = 0; i < pows2(width - 1); i++) {
                for (ulong j = 0; j < pows2(width - 1); j++) {
                    for (ulong jsite = 0; jsite < width; jsite++) {
                        HE(i, j) += -f * sigma_x_j(i, j, jsite, pows2, width, ignoremask)
                                - J * sigma_z_j_z_m(i, j, jsite, jsite + 1, pows2)*(((jsite + 1) < width))
                                - J2 * sigma_z_j_z_m(i, j, jsite, jsite + 2, pows2)*(((jsite + 2) < width))
                                - V * sigma_x_j_x_m(i, j, jsite, jsite + 1, pows2, width, ignoremask)*((jsite + 1) < width)
                                - Js[jsite % 2] * sigma_z_p_x_j_z_m(i, j, jsite, jsite + 1, jsite + 2, pows2, width, ignoremask)*(((jsite + 2) < width));
                        HO(i, j) += -f * sigma_x_j(~i, ~j, jsite, pows2, width, ignoremask)
                                - J * sigma_z_j_z_m(~i, ~j, jsite, jsite + 1, pows2)*(((jsite + 1) < width))
                                - J2 * sigma_z_j_z_m(~i, ~j, jsite, jsite + 2, pows2)*(((jsite + 2) < width))
                                - V * sigma_x_j_x_m(~i, ~j, jsite, jsite + 1, pows2, width, ignoremask)*((jsite + 1) < width)
                                - Js[jsite % 2] * sigma_z_p_x_j_z_m(~i, ~j, jsite, jsite + 1, jsite + 2, pows2, width, ignoremask)*(((jsite + 2) < width));
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
            std::string label = "Ising_L_" + to_string(width) + "_f_" + to_string(f)
	    + "_J2_"+to_string(J2)+"_V_" + to_string(V);

            std::ofstream outEs;
            if (WRITE_ENERGIES)
                outEs.open((label + "_energies").c_str(), std::ios::trunc);

            if (WRITE_OVERLAPS or WRITE_ALL_DECAY) {
                Arrayww overlap = Matrixww::Zero(pows2(width - 1), pows2(width - 1));
                for (int ispec = 0; ispec < pows2(width - 1); ispec++) {
                    auto spec = evecsE.col(ispec).array();
                    for (int i = 0; i < sigz.size(); i++) {
                        overlap(i, ispec) = (spec * sigz * evecsO.col(i).array()).sum();
                    }
                    
                    int maxind;
                    moverlaps[ispec] = overlap.col(ispec).abs().maxCoeff(&maxind);
                    meigdiff[ispec] = eigsE[ispec] - eigsO[maxind];
		        if (WRITE_VARS or WRITE_MEANS){
		      statoverlap.Push(moverlaps[ispec]);
		      stateigdiff.Push(meigdiff[ispec]);
			}
                    if (WRITE_ENERGIES)
                        outEs << eigsE[ispec] << '\n' << eigsO[ispec] << '\n';
                }
                if (WRITE_OVERLAPS){
                std::ofstream outoverlaps;
                outoverlaps.open((label + "_overlaps_states").c_str(), std::ios::trunc);
                outoverlaps << overlap;
                outoverlaps.flush();
                outoverlaps.close();
                }
                if (WRITE_ALL_DECAY){
                    Arrayw toverlaps (tparams.numPoints);
                    overlap = overlap.array().square();
                    Arrayww alleigdiff(pows2(width-1), pows2(width-1));
                    for (ulong i = 0; i < pows2(width - 1); i++) {
                        for (ulong j = 0; j < pows2(width - 1); j++) {
                            alleigdiff(i, j) = eigsE[j] - eigsO[i];
                        }
                    }
                    auto tfunc = [&](mpreal t, int j) {
                        toverlaps[j] = ((overlap.array()*alleigdiff.unaryExpr([ = ](mpreal x){return cos(x * t);})).sum()/(mpreal) pows2(width-1));
                        return false;
                    };
                    tfor.loop(tfunc, tparams);
                    plotter.writeToFile(label+"_long_decay", ts, toverlaps);
                }
            } else {
                Arrayw overlap = Vectorw::Zero(pows2(width - 1));
                for (int ispec = 0; ispec < pows2(width - 1); ispec++) {
                    auto spec = evecsE.col(ispec).array();
                    for (int i = 0; i < sigz.size(); i++) {
                        overlap(i) = (spec * sigz * evecsO.col(i).array()).sum();
                    }
                    int maxind;
                    moverlaps[ispec] = overlap.abs().maxCoeff(&maxind);
                    meigdiff[ispec] = eigsE[ispec] - eigsO[maxind];
		    if (WRITE_VARS or WRITE_MEANS){
		      statoverlap.Push(moverlaps[ispec]);
		      stateigdiff.Push(meigdiff[ispec]);
		    }
                    if (WRITE_ENERGIES)
                        outEs << eigsE[ispec] << '\n' << eigsO[ispec] << '\n';
                }
            }
            if (WRITE_ENERGIES) {
                outEs.flush();
                outEs.close();
            }
            
            if (WRITE_MAX_OVERLAPS)
                plotter.writeToFile(label + "_maxoverlaps", moverlaps);
            if (WRITE_PAIRED_EDIFFS)
                plotter.writeToFile(label + "_pairedeigdiffs", meigdiff);
            if (WRITE_PAIRED_DECAY) {
                Arrayw toverlaps(tparams.numPoints);
                moverlaps = moverlaps.square();
                auto tfunc = [&](mpreal t, int j) {
                    toverlaps[j] = (moverlaps * meigdiff.unaryExpr([ = ](mpreal x){return cos(x * t);})).mean();
                    return false;
                };
                tfor.loop(tfunc, tparams);
                plotter.writeToFile(label + "_paired_decay", ts, toverlaps);
            }
            if (WRITE_BULK) {
                Arrayww overlap = Matrixww::Zero(pows2(width - 1), width);
                const int bulkspec = pows2(width - 2);
                for (int jsite = 0; jsite < width; jsite++) {
                    Arrayw sigzn(pows2(width - 1));
                    for (int i = 0; i < sigzn.size(); i++) {
                        sigzn[i] = (2 * ((i & pows2(jsite)) >> jsite) - 1);
                    }
                    auto spec = evecsE.col(bulkspec).array();
                    for (int i = 0; i < sigz.size(); i++) {
                        overlap(i, jsite) += (spec * sigzn * evecsO.col(i).array()).sum();
                    }
                    std::ofstream outoverlaps;
                    outoverlaps.open((label + "_overlaps_bulk").c_str(), std::ios::trunc);
                    outoverlaps << overlap;
                    outoverlaps.flush();
                    outoverlaps.close();
                }
            }
            if (WRITE_MEANS) {
                maxoverlap.push_back(statoverlap.Mean());
                eigdiffs.push_back(stateigdiff.Mean());
            }
	    if (WRITE_VARS){
	      varoverlaps.push_back(statoverlap.Variance());
	      vareigdiffs.push_back(stateigdiff.Variance());
	    }
	    if (WRITE_VARS or WRITE_MEANS){
	      statoverlap.Clear();
	      stateigdiff.Clear();
	    }
            return false;
        };
        couplingsfor.loop(couplingsbody, fparams);
	std::string label = hashlabel + "Ising_L_" + std::to_string(width) + "_f_" + std::to_string(f) + "_V_" + std::to_string(V);
        if (WRITE_MEANS) {
            plotter.writeToFile(label + "_meanoverlap", fs, maxoverlap);
            plotter.writeToFile(label + "_meanediff", fs, eigdiffs);
	    maxoverlap.clear();
	    eigdiffs.clear();
        }
	if (WRITE_VARS){
	  plotter.writeToFile(label + "_varoverlap", fs, varoverlaps);
	  plotter.writeToFile(label + "_varediff", fs, vareigdiffs);
	  varoverlaps.clear();
	  vareigdiffs.clear();
	}
        return false;
    };

    NumMethod::Range forloop;
    forloop.loop(body, begin, end, 2);
    std::vector<int> widths = forloop.get_x(begin, end);

    return 0;

}
#endif 






