/*
 * File:   main.cpp
 * Author: Jack Kemp
 *
 * Created on 13 October 2014, 11:54
 */

#if 1

#include<Eigen/Dense>
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

using NumMethod::posmod;

//typedef mpfr::mpreal mpreal;
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
std::string to_string(T value) {
    return std::to_string(value);
}

template<>
std::string to_string(mpfr::mpreal value) {
    return std::to_string(value.toDouble());
}

/*
 *
 */
int main(int argc, char** argv) {
    mpfr::mpreal::set_default_prec(128);
    ScatterPlotter plotter;

    constexpr int maxwidth = 16;
    powersoftwo < maxwidth + 1 > pows2;

    const bool WRITE_ENERGIES = false;
    const bool WRITE_OVERLAPS = false;
    const bool WRITE_BULK = false;
    const bool WRITE_MEANS = false;
    const bool WRITE_VARS = false;
    const bool WRITE_MAX_OVERLAPS = false;
    const bool WRITE_PAIRED_EDIFFS = false;
    const bool WRITE_ALL_DECAY = false;
    const bool WRITE_PAIRED_DECAY = false;

    const bool FINITE_TEMPERATURE = false;
    const bool FINITE_TEMPERATURE_LOOP = false;
    const bool ENERGY_WINDOW_LOOP = true;


    const bool CMD_LINE_PARAMS = false;


    typedef Eigen::Matrix<mpreal, Eigen::Dynamic, Eigen::Dynamic> Matrixww;
    typedef Eigen::Matrix<mpreal, Eigen::Dynamic, 1> Vectorw;
    typedef Eigen::Array<mpreal, Eigen::Dynamic, 1> Arrayw;
    typedef Eigen::Array<mpreal, Eigen::Dynamic, Eigen::Dynamic> Arrayww;

    const mpreal r = 1.0;
    const int staggered = 0;
    const mpreal J2 = 0.0;
    const mpreal J3 = 0.0;
    const mpreal J4 = 0.0;
    const mpreal J = 1.0;
    const mpreal f = 0.05;
    const mpreal V = 0.0;
    const mpreal Js [] = {J3, J4};

    const mpreal T = 0.1;
    const mpreal Ew = 0.0;

    NumMethod::RunningStats<mpfr::mpreal> statoverlap, stateigdiff;
    NumMethod::WeightedRunningStats<mpfr::mpreal> wstatoverlap, wstateigdiff;

    const int begin = 8;
    const int end = 11;

    std::string hashlabel = "";
    if (CMD_LINE_PARAMS and argc > 4)
        hashlabel = std::string(argv[4]);

    NumMethod::ForLoopParams<mpreal> fparams;
    NumMethod::GetXFor<mpreal> recordx;
    NumMethod::EqualSpaceFor couplingsfor;
    fparams.start = 0.0;
    fparams.end = 3.5;
    fparams.numPoints = 8;

    if (CMD_LINE_PARAMS)
        fparams = NumMethod::get_for_from_cmd<mpreal>(argv);


    std::vector<mpfr::mpreal> maxoverlap, eigdiffs, varoverlaps, vareigdiffs;
    std::vector<mpreal> fs;
    couplingsfor.loop(recordx, fparams);
    fs = recordx.get_x();
    recordx.clear();

    NumMethod::EqualSpaceFor tfor;
    NumMethod::ForLoopParams<mpreal> tparams;
    std::vector<mpreal> ts;
    if (WRITE_ALL_DECAY | WRITE_PAIRED_DECAY) {
        tparams.start = 0.01;
        tparams.end = 1e6;
        tparams.numPoints = 1000;
        tfor.loop(recordx, tparams);
        ts = recordx.get_x();
        recordx.clear();
    }
    
    NumMethod::EqualSpaceFor Tfor;
    NumMethod::ForLoopParams<mpreal> Tparams;
    std::vector<mpreal> Ts;
    if (FINITE_TEMPERATURE_LOOP) {
        Tparams.start = 0.05;
        Tparams.end = 10;
        Tparams.numPoints = 500;
        Tfor.loop(recordx, Tparams);
        Ts = recordx.get_x();
        recordx.clear();
    }
    
    NumMethod::EqualSpaceFor Efor;
    NumMethod::ForLoopParams<mpreal> Eparams;
    std::vector<mpreal> Es;
    if (ENERGY_WINDOW_LOOP) {
        Eparams.start = 0.0;
        Eparams.end = 1;
        Eparams.numPoints = 500;
        Efor.loop(recordx, Eparams);
        Es = recordx.get_x();
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
            mpreal beta = 1.0/T;
            const mpreal Js [] = {J3, J4};
            std::vector<mpreal> J2s(width);

            for (int i = 0; i < width; ++i)
                J2s[i] = J2;
            J2s[staggered] *= r;
            HE = Matrixww::Zero(pows2(width - 1), pows2(width - 1));
            HO = Matrixww::Zero(pows2(width - 1), pows2(width - 1));

            //std::cout<< to_string(sigma_z_j(0, 0, 0, pows2)*((0+1) < width));

            for (ulong i = 0; i < pows2(width - 1); i++) {
                for (ulong j = 0; j < pows2(width - 1); j++) {
                    for (ulong jsite = 0; jsite < width; jsite++) {
                        HE(i, j) += -f * sigma_x_j(i, j, jsite, pows2, width, ignoremask)
                                - J * sigma_z_j_z_m(i, j, jsite, jsite + 1, pows2)*(((jsite + 1) < width))
                                - J2s[jsite] * sigma_z_j_z_m(i, j, jsite, jsite + 2, pows2)*(((jsite + 2) < width))
                                - V * sigma_x_j_x_m(i, j, jsite, jsite + 1, pows2, width, ignoremask)*((jsite + 1) < width)
                                - Js[jsite % 2] * sigma_z_p_x_j_z_m(i, j, jsite, jsite + 1, jsite + 2, pows2, width, ignoremask)*(((jsite + 2) < width));
                        HO(i, j) += -f * sigma_x_j(~i, ~j, jsite, pows2, width, ignoremask)
                                - J * sigma_z_j_z_m(~i, ~j, jsite, jsite + 1, pows2)*(((jsite + 1) < width))
                                - J2s[jsite] * sigma_z_j_z_m(~i, ~j, jsite, jsite + 2, pows2)*(((jsite + 2) < width))
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
            mpreal partE;
            std::string label = "Ising_L_" + to_string(width) + "_f_" + to_string(f)
                    + "_J2_" + to_string(J2);

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
                    if (WRITE_VARS or WRITE_MEANS) {
                        if (FINITE_TEMPERATURE) {
                            mpreal weight = exp(-beta * eigsE[ispec]);
                            wstatoverlap.Push(moverlaps[ispec], weight);
                            wstateigdiff.Push(meigdiff[ispec], weight);
                        } else {
                            statoverlap.Push(moverlaps[ispec]);
                            stateigdiff.Push(meigdiff[ispec]);
                        }
                    }
                    if (WRITE_ENERGIES)
                        outEs << eigsE[ispec] << '\n' << eigsO[ispec] << '\n';
                }
                if (WRITE_OVERLAPS) {
                    std::ofstream outoverlaps;
                    outoverlaps.open((label + "_overlaps_states").c_str(), std::ios::trunc);
                    outoverlaps << overlap;
                    outoverlaps.flush();
                    outoverlaps.close();
                }
                if (WRITE_ALL_DECAY) {
                    Arrayw toverlaps(tparams.numPoints);
                    overlap = overlap.array().square();
                    Arrayww alleigdiff(pows2(width - 1), pows2(width - 1));
                    for (ulong i = 0; i < pows2(width - 1); i++) {
                        for (ulong j = 0; j < pows2(width - 1); j++) {
                            alleigdiff(i, j) = eigsE[j] - eigsO[i];
                        }
                    }
                    auto tfunc = [&](mpreal t, int j) {
                        toverlaps[j] = ((overlap.array() * alleigdiff.unaryExpr([ = ](mpreal x){return cos(x * t);})).sum() / (mpreal) pows2(width - 1));
                        return false;
                    };
                    tfor.loop(tfunc, tparams);
                    plotter.writeToFile(label + "_long_decay", ts, toverlaps);
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
                    if (WRITE_VARS or WRITE_MEANS) {
                        if (FINITE_TEMPERATURE) {
                            mpreal weight = exp(-beta * eigsE[ispec]);
                            wstatoverlap.Push(moverlaps[ispec], weight);
                            wstateigdiff.Push(meigdiff[ispec], weight);
                        } else {
                            statoverlap.Push(moverlaps[ispec]);
                            stateigdiff.Push(meigdiff[ispec]);
                        }
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
            if (FINITE_TEMPERATURE_LOOP){
                Arrayw weights (pows2(width-1));
                Arrayw mean_overlaps(Tparams.numPoints);
                Arrayw variance_overlaps(Tparams.numPoints);
                Arrayw mean_ediff(Tparams.numPoints);
                Arrayw variance_ediff(Tparams.numPoints);
                auto Tfunc = [&](mpreal T, int j) {
                    weights = eigsE.unaryExpr([ = ](mpreal x){return exp(-x/T);});
                    NumMethod::WeightedRunningStats<mpreal> Tstatoverlap(moverlaps, weights), Tstatediff(meigdiff, weights);
                    mean_overlaps[j] = Tstatoverlap.Mean();
                    variance_overlaps[j]= Tstatoverlap.Variance();
                    mean_ediff[j] = Tstatediff.Mean();
                    variance_ediff[j] = Tstatediff.Variance();
                    return false;
                };
                Tfor.loop(Tfunc, Tparams);
                plotter.writeToFile(label + "_Tmeanoverlap", Ts, mean_overlaps);
                plotter.writeToFile(label + "_Tvaroverlap", Ts, variance_overlaps);
                plotter.writeToFile(label + "_Tmeanediff", Ts, mean_ediff);
                plotter.writeToFile(label + "_Tvarediff", Ts, variance_ediff);
            }
            if (ENERGY_WINDOW_LOOP){
                Arrayw weights (pows2(width-1));
                Arrayw mean_overlaps(Eparams.numPoints);
                Arrayw variance_overlaps(Eparams.numPoints);
                Arrayw mean_ediff(Eparams.numPoints);
                Arrayw variance_ediff(Eparams.numPoints);
                auto Efunc = [&](mpreal Efraction, int j) {
                    mpreal E = Efraction*(eigsE[pows2(width-1)-1]-eigsE[0])+eigsE[0];
                    weights = eigsE.unaryExpr([ = ](mpreal x){return x < E ? 1.0 : 0.0;});
                    NumMethod::WeightedRunningStats<mpreal> Estatoverlap(moverlaps, weights), Estatediff(meigdiff, weights);
                    mean_overlaps[j] = Estatoverlap.Mean();
                    variance_overlaps[j]= Estatoverlap.Variance();
                    mean_ediff[j] = Estatediff.Mean();
                    variance_ediff[j] = Estatediff.Variance();
                    return false;
                };
                Efor.loop(Efunc, Eparams);
                plotter.writeToFile(label + "_Emeanoverlap", Es, mean_overlaps);
                plotter.writeToFile(label + "_Evaroverlap", Es, variance_overlaps);
                plotter.writeToFile(label + "_Emeanediff", Es, mean_ediff);
                plotter.writeToFile(label + "_Evarediff", Es, variance_ediff);
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
                if (FINITE_TEMPERATURE){
                    maxoverlap.push_back(wstatoverlap.Mean());
                    eigdiffs.push_back(wstateigdiff.Mean());
                } else {
                    maxoverlap.push_back(statoverlap.Mean());
                    eigdiffs.push_back(stateigdiff.Mean());
                }
            }
            if (WRITE_VARS) {
                if (FINITE_TEMPERATURE){
                    varoverlaps.push_back(wstatoverlap.Variance());
                    vareigdiffs.push_back(wstateigdiff.Variance());
                } else {
                    varoverlaps.push_back(statoverlap.Variance());
                    vareigdiffs.push_back(stateigdiff.Variance());
                }
            }
            if (WRITE_VARS or WRITE_MEANS) {
                statoverlap.Clear();
                stateigdiff.Clear();
                wstatoverlap.Clear();
                wstateigdiff.Clear(); 
            }
            return false;
        };
        couplingsfor.loop(couplingsbody, fparams);

        std::string label = hashlabel + "Ising_windowtest_L_" + to_string(width)
                + "_V_" + to_string(V)+"_f_"+to_string(f);
        if (WRITE_MEANS) {
            plotter.writeToFile(label + "_meanoverlap", fs, maxoverlap);
            plotter.writeToFile(label + "_meanediff", fs, eigdiffs);
            maxoverlap.clear();
            eigdiffs.clear();
        }
        if (WRITE_VARS) {
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
