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
#include<complex>
#include <bitset>

using NumMethod::posmod;

//typedef mpfr::mpreal mpreal;
typedef double mpreal;
typedef std::complex<mpreal> complex;
typedef unsigned long ulong;

const int ipow(int base, int exp)
{
    int result = 1;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }

    return result;
}
class Spin3State {
    static constexpr ulong mask = 3;
    ulong state;
    
public:
    static ulong width;
    ulong operator [] (const ulong& i) const{
        return (this->state >> (2*i)) & this->mask;
    }
    ulong get_state() const{
        return this->state;
    }
    Spin3State increment() const{
        if (!(this->state & 2))
            return Spin3State(this->state+1);
        ulong incr = 2;
        for (int i=2; (this->state >> i) & 2; i+=2){
            incr += (1 << i);
        }
        return Spin3State(this->state + incr); 
    }
    
    Spin3State cycle(const ulong& i) const{
        return (*this)[i] == 2 ? Spin3State(this->state - (2 << (2*i))) : Spin3State(this->state + (1 << (2*i)));
    }
    
    Spin3State cycle_all() const{
        ulong new_state = this->state;
        for (ulong i =0; i < width; i++)
          new_state += ((*this)[i] == 2 ? - (2 << (2*i)) : (1 << (2*i)));
        return Spin3State(new_state);
    }
    
    Spin3State bcycle_all() const{
        ulong new_state = this->state;
        for (ulong i =0; i < width; i++)
          new_state += ((*this)[i] == 0 ? (2 << (2*i)) : - (1 << (2*i)));
        return Spin3State(new_state);
    }
    
    Spin3State bcycle(const ulong& i) const{
        return (*this)[i] == 0 ? Spin3State(this->state +(2 << (2*i))) : Spin3State(this->state - (1 << (2*i)));
    }
    
    ulong cycle_get(const ulong& i) const{
        ulong shifted = (this->state >> (2*i))& this->mask;
        return shifted & 2 ? 0 : (shifted+1);
    }
    
    ulong bcycle_get(const ulong& i) const{
        ulong shifted = (this->state >> (2*i)) & this->mask;
        return shifted == 0 ? 2 : (shifted-1);
    }
    
    Spin3State(ulong istate): state(istate){
    }
    
    bool operator == (const Spin3State& other) const{
        return this->state == other.state;
    }
    
    Spin3State operator = (const ulong& x) {
        this->state = x;
        return *this;
    }
    
};

ulong Spin3State::width;

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

const double pi = 3.141592653589793238;

template <int N, typename mpreal>
class RootOfUnity {
    std::complex<mpreal> table [2*N];
public:

     complex operator()(int i) const{
        return table [i];
    }
     
    complex operator[](int i) const{
        return table [i];
    }


    const int * getTable() const{
        return table;
    };

    RootOfUnity() {
        for (int i = 0; i < N; i++) {
            table [i] = std::polar(1., (2*pi*i)/N);
            table [i+N] = table[i]; //allows wrap around
        } 
    }
};

RootOfUnity<3, double> omega;


inline complex sigma_z_j(const Spin3State& x, const Spin3State& y, const ulong& j) {
    return ((mpreal)(x == y))*omega[y[j]];
}


inline complex sigma_x_j(const Spin3State& x, const Spin3State& y, const ulong& j, const ulong& width) {
    return j== width -1 ? ((mpreal) (x==y.bcycle_all().cycle(width-1)))*omega[y[j]^3]: (mpreal)(x == y.cycle(j));
}


inline complex sigma_x_j_conj(const Spin3State& x, const Spin3State& y, const ulong& j, const ulong& width) {
    return j== width -1 ? ((mpreal)(x==y.cycle_all().bcycle(width-1)))*omega[y[j]]: (mpreal) (x == y.bcycle(j));
}


inline complex sigma_z_j_conj_z_m(const Spin3State& x, const Spin3State& y, ulong j, ulong m) {
    return ((mpreal) (x == y))*omega[(y[j]^3) + y[m]];
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
    const bool WRITE_MEANS = true;
    const bool WRITE_VARS = true;
    const bool WRITE_MAX_OVERLAPS = false;
    const bool WRITE_PAIRED_EDIFFS = false;
    const bool WRITE_ALL_DECAY = false;
    const bool WRITE_PAIRED_DECAY = false;
    const bool WRITE_ALL_LEVEL_SPACINGS = false;
    const bool WRITE_MEAN_LEVEL_SPACINGS = false;

    const bool FINITE_TEMPERATURE = false;
    const bool FINITE_TEMPERATURE_LOOP = false;
    const bool ENERGY_WINDOW_LOOP = false;


    const bool CMD_LINE_PARAMS = false;


    typedef Eigen::Matrix<complex, Eigen::Dynamic, Eigen::Dynamic> Matrixww;
    typedef Eigen::Matrix<complex, Eigen::Dynamic, 1> Vectorw;
    typedef Eigen::Array<complex, Eigen::Dynamic, 1> ArraywC;
    typedef Eigen::Array<mpreal, Eigen::Dynamic, 1> Arrayw;
    typedef Eigen::Array<mpreal, Eigen::Dynamic, Eigen::Dynamic> Arrayww;

    const double phi = pi / 5;
    const double theta = pi / 6;
    complex ephi = std::polar(1., phi);
    complex mephi = std::polar(1., -phi);
    complex eth = std::polar(1., theta);
    complex meth = std::polar(1., -theta);
    const double J = 1.0;
    const double f = 0.01;
    
    const double T =0;

    NumMethod::RunningStats<mpfr::mpreal> statoverlap, stateigdiff;
    NumMethod::WeightedRunningStats<mpfr::mpreal> wstatoverlap, wstateigdiff;

    const int begin = 6;
    const int end = 8;

    std::string hashlabel = "";
    if (CMD_LINE_PARAMS and argc > 4)
        hashlabel = std::string(argv[4]);

    NumMethod::ForLoopParams<mpreal> fparams;
    NumMethod::GetXFor<mpreal> recordx;
    NumMethod::EqualSpaceFor couplingsfor;
    fparams.start = 0.0;
    fparams.end = 3.15/2;
    fparams.numPoints = 250;

    if (CMD_LINE_PARAMS)
        fparams = NumMethod::get_for_from_cmd<mpreal>(argv);


    std::vector<mpfr::mpreal> maxoverlap, eigdiffs, varoverlaps, vareigdiffs, meanspacings;
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
        Spin3State::width = width;
        const int nstates = ipow(3, width);
        const int n3states = ipow(3, width-1);

        Arrayw moverlaps(n3states), meigdiff(n3states);

        int ignoremask = ~((1 << (width - 1)) - 1);

        Matrixww HA(n3states, n3states);
        Matrixww HB(n3states, n3states);

        
        Arrayw sigz2(n3states);
        Arrayw sigz(n3states);
        for (int i = 0; i < sigz.size(); i++) {
            //sigz[i] = omega[i % 3];
            sigz[i] = 2*cos(2*pi*(i%3)/3.);
        }

        auto couplingsbody = [&](mpreal theta, int j) {
            mpreal beta = 1.0 / T;
            complex ephi = std::polar(1., phi);
            complex mephi = std::polar(1., -phi);
            complex eth = std::polar(1., theta);
            complex meth = std::polar(1., -theta);
            

            //std::cout<< to_string(sigma_z_j(0, 0, 0, pows2)*((0+1) < width));
            HA = Matrixww::Zero(n3states, n3states);
            HB = Matrixww::Zero(n3states, n3states);

            Spin3State istate(0);
            Spin3State jstate(0);
            for (ulong i = 0; i < n3states; i++) {
                jstate = 0;
                for (ulong j = 0; j < n3states; j++) {
                    for (ulong jsite = 0; jsite < width; jsite++) {
                        HA(i, j) += -f * (ephi * sigma_x_j(istate, jstate, jsite, width) + mephi * sigma_x_j_conj(istate, jstate, jsite, width))
                                - J * (meth * sigma_z_j_conj_z_m(istate, jstate, jsite, jsite + 1)
                                + eth * sigma_z_j_conj_z_m(istate, jstate, jsite + 1, jsite))*((mpreal) ((jsite + 1) < width));
                        Spin3State istateB = istate.cycle_all(), jstateB = jstate.cycle_all();
                        HB(i, j) += -f * (ephi * sigma_x_j(istateB, jstateB, jsite, width) + mephi * sigma_x_j_conj(istateB, jstateB, jsite, width))
                                - J * (meth * sigma_z_j_conj_z_m(istateB, jstateB, jsite, jsite + 1)
                                + eth * sigma_z_j_conj_z_m(istateB, jstateB, jsite + 1, jsite))*((mpreal) ((jsite + 1) < width));
                        //Spin3State istateC = istate.bcycle_all(), jstateC = jstate.bcycle_all();  
                        //HC(i, j) += -f *(ephi*sigma_x_j(istateC, jstateC, jsite) + mephi*sigma_x_j_conj(istateC, jstateC, jsite))
                        //        - J *(meth*sigma_z_j_conj_z_m(istateC, jstateC, jsite, jsite + 1)
                        //             +eth*sigma_z_j_conj_z_m(istateC, jstateC, jsite+1, jsite))*((mpreal) ((jsite + 1) < width));
                    }
                    jstate = jstate.increment();
                    //std::cout << jstate.get_state();
                }
                istate = istate.increment();
            }

            //std::cout << H << std::endl << std::endl;
            Eigen::SelfAdjointEigenSolver<Matrixww> esE(HA);
            Eigen::SelfAdjointEigenSolver<Matrixww> esO(HB);
            auto eigsE = esE.eigenvalues();
            auto eigsO = esO.eigenvalues();
            auto evecsE = esE.eigenvectors();
            auto evecsO = esO.eigenvectors();
            mpreal partE;
            std::string label = "parafermion_L_" + to_string(width) + "_f_" + to_string(f)
                    + "_J_" + to_string(J);

            std::ofstream outEs;
            if (WRITE_ENERGIES)
                outEs.open((label + "_energies").c_str(), std::ios::trunc);

            if (WRITE_OVERLAPS or WRITE_ALL_DECAY) {
                Arrayww overlap = Arrayww::Zero(n3states, n3states);
                for (int ispec = 0; ispec < n3states; ispec++) {
                    auto spec = evecsE.col(ispec).array();
                    for (int i = 0; i < sigz.size(); i++) {
                        overlap(i, ispec) = fabs( (spec * sigz * evecsO.col(i).array()).sum() );
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
                    Arrayww alleigdiff(n3states, n3states);
                    for (ulong i = 0; i < n3states; i++) {
                        for (ulong j = 0; j < n3states; j++) {
                            alleigdiff(i, j) = eigsE[j] - eigsO[i];
                        }
                    }
                    auto tfunc = [&](mpreal t, int j) {
                        toverlaps[j] = ((overlap.array() * alleigdiff.unaryExpr([ = ](mpreal x){return cos(x * t);})).sum() / (mpreal) n3states);
                        return false;
                    };
                    tfor.loop(tfunc, tparams);
                    plotter.writeToFile(label + "_long_decay", ts, toverlaps);
                }
            } else {
                Arrayw overlap = Arrayw::Zero(n3states);
                for (int ispec = 0; ispec < n3states; ispec++) {
                    auto spec = evecsE.col(ispec).array();
                    for (int i = 0; i < sigz.size(); i++) {
                        overlap(i) = fabs( (spec * sigz * evecsO.col(i).array()).sum() );
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
            if (WRITE_ALL_LEVEL_SPACINGS or WRITE_MEAN_LEVEL_SPACINGS){
                //Arrayw energies(nstates-1);
                //for (int i=0;i < n3states; i++){
                    //energies[2*i] = eigsE[i];
                    //energies[2*i+1] = eigsO[i];
                //}
                //std::sort(energies.data(), energies.data()+energies.size());
                //Arrayw spacings = (energies.head(nstates-1)-energies.tail(nstates-1)).abs();
                Arrayw spacings = (eigsE.array().head(n3states-1)-eigsE.array().tail(n3states-1)).abs();
                Arrayw normspacings (n3states-2);
                for (int i=0; i< n3states-2; i++){
                    normspacings[i] = NumMethod::min(spacings[i], spacings[i+1])/NumMethod::max(spacings[i],spacings[i+1]);
                   // normspacings[i] = spacings[i];
                }
                if (WRITE_ALL_LEVEL_SPACINGS){
                    plotter.writeToFile(label + "_level_spacings", normspacings);
                }
                if (WRITE_MEAN_LEVEL_SPACINGS){
                    meanspacings.push_back(normspacings.mean());
                }
            }
            if (FINITE_TEMPERATURE_LOOP){
                Arrayw weights (n3states);
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
                Arrayw weights (n3states);
                Arrayw mean_overlaps(Eparams.numPoints);
                Arrayw variance_overlaps(Eparams.numPoints);
                Arrayw mean_ediff(Eparams.numPoints);
                Arrayw variance_ediff(Eparams.numPoints);
                auto Efunc = [&](mpreal Efraction, int j) {
                    mpreal E = Efraction*(eigsE[n3states-1]-eigsE[0])+eigsE[0];
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
                std::cerr << "NOT IMPLEMENTED FOR PARAFERMIONS";
                /*Arrayww overlap = Matrixww::Zero(n3states, width);
                const int bulkspec = pows2(width - 2);
                for (int jsite = 0; jsite < width; jsite++) {
                    Arrayw sigzn(n3states);
                    for (int i = 0; i < sigzn.size(); i++) {
                        sigzn[i] = (2 * ((i & pows2(jsite)) >> jsite) - 1);
                    }
                    auto spec = evecsE.col(bulkspec).array();
                    for (int i = 0; i < sigz.size(); i++) {
                        overlap(i, jsite) += fabs((spec * sigzn * evecsO.col(i).array()).sum());
                    }
                    std::ofstream outoverlaps;
                    outoverlaps.open((label + "_overlaps_bulk").c_str(), std::ios::trunc);
                    outoverlaps << overlap;
                    outoverlaps.flush();
                    outoverlaps.close();
                 */
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

        std::string label = hashlabel + "parafermion_L_" + to_string(width)
                + "_J_" + to_string(J)+"_f_"+to_string(f) + "_phi_" + to_string(phi);
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
        if (WRITE_MEAN_LEVEL_SPACINGS){
            plotter.writeToFile(label + "_meanspacings", fs, meanspacings);
            meanspacings.clear();
        }
        return false;
    };

    NumMethod::Range forloop;
    forloop.loop(body, begin, end, 1);
    std::vector<int> widths = forloop.get_x(begin, end);

    return 0;

}
#endif
