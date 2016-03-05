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
#include<complex>
#include <bitset>

using NumMethod::posmod;

typedef double mpreal;
typedef std::complex<mpreal> complex;
typedef unsigned long ulong;


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
    
    Spin3State operator = (ulong x) {
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



inline complex sigma_z_j(const Spin3State& x, const Spin3State& y, ulong j) {
    return ((mpreal)(x == y))*omega[y[j]];
}


inline mpreal sigma_x_j(const Spin3State& x, const Spin3State& y, ulong j) {
    return (mpreal) (x == y.cycle(j));
}


inline mpreal sigma_x_j_conj(const Spin3State& x, const Spin3State& y, ulong j) {
    return (mpreal) (x == y.bcycle(j));
}


inline complex sigma_z_j_conj_z_m(const Spin3State& x, const Spin3State& y, ulong j, ulong m) {
    return ((mpreal) (x == y))*omega[(y[j]^3) + y[m]];
}



/*
 * 
 */
int main() {
    //mpreal::set_default_prec(128);
    ScatterPlotter plotter;
    const int width = 3;
    const int nstates = ipow(3, width);

    typedef Eigen::Matrix<complex, Eigen::Dynamic, Eigen::Dynamic> Matrixww;
    typedef Eigen::Matrix<complex, Eigen::Dynamic, 1> Vectorw;
    typedef Eigen::Array<mpreal, Eigen::Dynamic, 1> Arrayw;

    const double phi = 0.0;
    const double theta = 0.0;
    complex ephi = std::polar(1., phi);
    complex mephi = std::polar(1., -phi);
    complex eth = std::polar(1., theta);
    complex meth = std::polar(1., -theta);
    const double J = 1.0;
    const double f = 0.0;
    
    std::vector<mpreal> maxoverlap, fs, eigdiffs;
    
    
    
    Spin3State tstate(9);
    std::cout<< tstate[1] << " " << omega[1^3] << " " << omega[2^3] << std::endl;
    std::cout<< tstate.cycle(1).get_state() << " " << tstate.cycle(0).get_state() << " " << omega[2] << std::endl;
    

    Eigen::Matrix<mpreal, Eigen::Dynamic, Eigen::Dynamic> test (nstates, nstates);

    Spin3State istate(0);
    Spin3State jstate(0);
    for (ulong i = 0; i < nstates; i++) {
        jstate = 0;
        for (ulong j = 0; j < nstates; j++) {
            test(i, j) = sigma_x_j_conj(istate, jstate, 2);
            jstate = jstate.increment();
            //std::cout << jstate.get_state();
        }
        istate = istate.increment();
    }

    std::cout << test << std::endl;
    
    Matrixww H(nstates, nstates);
   
    NumMethod::ForLoopParams<mpreal> fparams;
    fparams.numPoints = 2; fparams.start = 0.0; fparams.end = 0.6;
    auto body = [&](mpreal f, int i) {
        H = Matrixww::Zero(nstates, nstates);
        Spin3State istate(0);
        Spin3State jstate(0);
        for (ulong i = 0; i < nstates; i++) {
            jstate = 0;
            for (ulong j = 0; j < nstates; j++) {
                for (ulong jsite = 0; jsite < width; jsite++) {
                    H(i, j) += -f *(ephi*sigma_x_j(istate, jstate, jsite) + mephi*sigma_x_j_conj(istate, jstate, jsite))
                            - J *(meth*sigma_z_j_conj_z_m(istate, jstate, jsite, jsite + 1)
                                 +eth*sigma_z_j_conj_z_m(istate, jstate, jsite+1, jsite))*((mpreal) ((jsite + 1) < width));   
                }
                jstate = jstate.increment();
                //std::cout << jstate.get_state();
            }
            istate = istate.increment();
        }

        Eigen::SelfAdjointEigenSolver<Matrixww> es(H);
        auto eigs = es.eigenvalues();
        plotter.writeToFile("parafermionsall_"+std::to_string(f),eigs);
        return false;
    };
    
    NumMethod::EqualSpaceFor forloop; 
    forloop.loop(body, fparams);
    return 0;

}
#endif 


