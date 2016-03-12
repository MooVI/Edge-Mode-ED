/* 
 * File:   main.cpp
 * Author: Jack Kemp
 *
 * Created on 13 October 2014, 11:54
 */

#if 1

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
#include<complex>
#include <bitset>
#include<thread>

using NumMethod::posmod;

typedef double mpreal;
typedef std::complex<mpreal> complex;
typedef unsigned long ulong;

typedef Eigen::Matrix<complex, Eigen::Dynamic, Eigen::Dynamic> Matrixww;
typedef Eigen::Matrix<complex, Eigen::Dynamic, 1> Vectorw;
typedef Eigen::Array<complex, Eigen::Dynamic, 1> ArraywC;
typedef Eigen::Array<mpreal, Eigen::Dynamic, 1> Arrayw;

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

const int width = 7;
const int nstates = ipow(3, width);
const int n3states = ipow(3, width-1);

class Spin3State {
    static constexpr ulong mask = 3;
    ulong state;
public:
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






inline complex sigma_z_j(const Spin3State& x, const Spin3State& y, ulong j) {
    return ((mpreal)(x == y))*omega[y[j]];
}


inline complex sigma_x_j(const Spin3State& x, const Spin3State& y, const ulong& j) {
    return j== width -1 ? ((mpreal) (x==y.bcycle_all().cycle(width-1)))*omega[y[j]^3]: (mpreal)(x == y.cycle(j));
}


inline complex sigma_x_j_conj(const Spin3State& x, const Spin3State& y, const ulong& j) {
    return j== width -1 ? ((mpreal)(x==y.cycle_all().bcycle(width-1)))*omega[y[j]]: (mpreal) (x == y.bcycle(j));
}


inline complex sigma_z_j_conj_z_m(const Spin3State& x, const Spin3State& y, ulong j, ulong m) {
    return ((mpreal) (x == y))*omega[(y[j]^3) + y[m]];
}

template<typename F>
void solve_eigensystem(Eigen::SelfAdjointEigenSolver<Matrixww>& es, const Eigen::MatrixBase<F> & M){
    es.compute(M);
}

/*
 * 
 */
int main() {
    //mpreal::set_default_prec(128);
    ScatterPlotter plotter;
    //Eigen::initParallel();
    

    const double phi = pi/6;
    const double theta = pi/3;
    complex ephi = std::polar(1., phi);
    complex mephi = std::polar(1., -phi);
    complex eth = std::polar(1., theta);
    complex meth = std::polar(1., -theta);
    const double J = 1.0;
    const double f = 0.01;
    
    NumMethod::RunningStats<mpreal> statoverlapA, stateigdiffA, statoverlapB, stateigdiffB ;
    
    std::vector<mpreal> maxoverlapA, fs, eigdiffsA, maxoverlapB, eigdiffsB;
    
    Arrayw sigz(n3states);
        for (int i = 0; i < sigz.size(); i++) {
            //sigz[i] = omega[i % 3];
            sigz[i] = 2*cos(2*pi*(i%3)/3.);
        }
    
    Spin3State tstate(9);
    std::cout<< tstate[1] << " " << omega[1^3] << " " << omega[2^3] << std::endl;
    std::cout<< tstate.cycle(1).get_state() << " " << tstate.cycle(0).get_state() << " " << omega[2] << std::endl;
    

    Eigen::Matrix<complex, Eigen::Dynamic, Eigen::Dynamic> test (n3states, n3states);

//    Spin3State istate(0);
//    Spin3State jstate(0);
//    for (ulong i = 0; i < n3states; i++) {
//        jstate = 0;
//        for (ulong j = 0; j < n3states; j++) {
//            test(i, j) = sigma_x_j_conj(istate.bcycle_all(), jstate.bcycle_all(), 2);
//            jstate = jstate.increment();
//            //std::cout << jstate.get_state();
//        }
//        istate = istate.increment();
//    }
//
//    std::cout << test << std::endl;
    
    Matrixww HA(n3states, n3states);
    Matrixww HB(n3states, n3states);
    Matrixww HC(n3states, n3states);
    NumMethod::ForLoopParams<mpreal> fparams;
    fparams.numPoints = 100; fparams.start = 0; fparams.end = 1.1;
    auto body = [&](mpreal f, int i) {
        HA = Matrixww::Zero(n3states, n3states);
        HB = Matrixww::Zero(n3states, n3states);
        HC = Matrixww::Zero(n3states, n3states);
        Spin3State istate(0);
        Spin3State jstate(0);
        for (ulong i = 0; i < n3states; i++) {
            jstate = 0;
            for (ulong j = 0; j < n3states; j++) {
                for (ulong jsite = 0; jsite < width; jsite++) {
                    HA(i, j) += -f *(ephi*sigma_x_j(istate, jstate, jsite) + mephi*sigma_x_j_conj(istate, jstate, jsite))
                            - J *(meth*sigma_z_j_conj_z_m(istate, jstate, jsite, jsite + 1)
                                 +eth*sigma_z_j_conj_z_m(istate, jstate, jsite+1, jsite))*((mpreal) ((jsite + 1) < width));
                    Spin3State istateB = istate.cycle_all(), jstateB = jstate.cycle_all();  
                    HB(i, j) += -f *(ephi*sigma_x_j(istateB, jstateB, jsite) + mephi*sigma_x_j_conj(istateB, jstateB, jsite))
                            - J *(meth*sigma_z_j_conj_z_m(istateB, jstateB, jsite, jsite + 1)
                                 +eth*sigma_z_j_conj_z_m(istateB, jstateB, jsite+1, jsite))*((mpreal) ((jsite + 1) < width));
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

        
        Eigen::SelfAdjointEigenSolver<Matrixww> esA(HA.rows());
        std::thread t1(solve_eigensystem<Matrixww>,std::ref(esA), std::ref(HA));
        //Eigen::SelfAdjointEigenSolver<Matrixww> esB(HB.rows());
        //std::thread t2(solve_eigensystem<Matrixww>,std::ref(esB), std::ref(HB));
        Eigen::SelfAdjointEigenSolver<Matrixww> esB(HB);
        t1.join();
        //t2.join();
        Arrayw overlapA = Arrayw::Zero(n3states), overlapB = Arrayw::Zero(n3states);
        auto evecsA = esA.eigenvectors();
        auto evecsB = esB.eigenvectors();
        //auto evecsC = esC.eigenvectors();
        auto eigsA = esA.eigenvalues();
        auto eigsB = esB.eigenvalues();
        //auto eigsC = esC.eigenvalues();
        for (int ispec =0; ispec< n3states; ispec++){
        //const int ispec = n3states/3;
        //const int ispec = 0;
        auto specA = evecsA.col(ispec).array();
        //auto specB = evecsB.col(ispec).array();
        for (int i = 0; i < sigz.size(); i++) {
            overlapA[i] =fabs((specA*sigz*evecsB.col(i).array()).sum());
            //overlapB[i] =fabs((specB*sigz*evecsC.col(i).array()).sum()); 
        }
        
        int maxindA, maxindB;
        statoverlapA.Push(overlapA.abs().maxCoeff(&maxindA));
        stateigdiffA.Push(fabs(eigsA[ispec]-eigsB[maxindA]));
        //statoverlapB.Push(overlapB.abs().maxCoeff(&maxindB));
        //stateigdiffB.Push(fabs(eigsB[ispec]-eigsC[maxindB]));
        }
        fs.push_back(f);
        maxoverlapA.push_back(statoverlapA.Mean());
                eigdiffsA.push_back(stateigdiffA.Mean());
        //        maxoverlapB.push_back(statoverlapB.Mean());
        //        eigdiffsB.push_back(stateigdiffB.Mean());
        //plotter.plot2(esB.eigenvalues(), overlapA);
        //plotter.plot2(esC.eigenvalues(), overlapB);
        //plotter.writeToFile("parafermions_"+std::to_string(f), esA.eigenvalues(), esB.eigenvalues(), esC.eigenvalues());
        return false;
    };
    NumMethod::EqualSpaceFor forloop; 
    forloop.loop(body, fparams);
    PlotterData pdA, pdB, pdAe, pdBe;
    pdA.style = "l";
    pdA.input = "parfermion_mean_A_" + std::to_string(width);
    pdB.input = "parfermion_mean_B_" + std::to_string(width);
    //plotter.plot2(fs, maxoverlap, pd);
    pdAe.input = pdA.input + std::string("_eigdiff");
    pdBe.input = pdB.input + std::string("_eigdiff");
    plotter.plot2(fs, eigdiffsA, pdAe);
    plotter.plot2(fs, maxoverlapA, pdA);
    //plotter.plot2(fs, eigdiffsB, pdBe);
    //plotter.plot2(fs, maxoverlapB, pdB);
    plotter.wait();
    return 0;

}
#endif 




