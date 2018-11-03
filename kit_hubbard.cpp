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

template <int N>
inline int sigma_y_j_y_m(ulong x, ulong y, ulong j, ulong m, powersoftwo<N> pows, const int width, const int ignoremask) {
    return j == m ? (x==y ? 1 : 0)  :
            ( 
            (x == (y^pows(j)^pows(m)))
            or
            ( (j == width-1) and ((x | ignoremask) == ((~(y^pows(m))) | ignoremask)))
            or
            ( (m == width-1) and ((x | ignoremask) == ((~(y^pows(j))) | ignoremask)))
            )*(-2 * ((pows(j) & y & pows(width-1)) >> (width-1)) + 1)
             *(-2 * ((pows(m) & y & pows(width-1)) >> (width-1)) + 1)
             *(2*(((y & pows(j)) >> j) ^ ((y & pows(m)) >> m))-1);
}

template <int N>
inline int sigma_x_j_x_k_y_l_y_m(ulong x, ulong y, ulong j, ulong k, ulong l, ulong m,  powersoftwo<N> pows, const int width, const int ignoremask) {
    //assumes j,k,l,m different
    return 
            ( 
            (x == (y^pows(j)^pows(m)^pows(l)^pows(k)))
            or
            ( (j == width-1) and ((x | ignoremask) == ((~(y^pows(m)^pows(l)^pows(k))) | ignoremask)))
            or
            ( (m == width-1) and ((x | ignoremask) == ((~(y^pows(j)^pows(l)^pows(k))) | ignoremask)))
            or
            ( (k == width-1) and ((x | ignoremask) == ((~(y^pows(j)^pows(l)^pows(m))) | ignoremask)))
            or
            ( (l == width-1) and ((x | ignoremask) == ((~(y^pows(j)^pows(m)^pows(k))) | ignoremask)))
            )*(-2 * ((pows(l) & y & pows(width-1)) >> (width-1)) + 1)
             *(-2 * ((pows(m) & y & pows(width-1)) >> (width-1)) + 1)
             *(2*(((y & pows(l)) >> l) ^ ((y & pows(m)) >> m))-1)
            *(-2 * ((pows(j) & y & pows(width - 1)) >> (width - 1)) + 1)
            *(-2 * ((pows(k) & y & pows(width - 1)) >> (width - 1)) + 1);
            
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
    const bool PERIODIC = false;
    
    const bool CMD_LINE_PARAMS = false;


    typedef Eigen::Matrix<mpreal, Eigen::Dynamic, Eigen::Dynamic> Matrixww;
    typedef Eigen::Matrix<mpreal, Eigen::Dynamic, 1> Vectorw;
    typedef Eigen::Array<mpreal, Eigen::Dynamic, 1> Arrayw;
    typedef Eigen::Array<mpreal, Eigen::Dynamic, Eigen::Dynamic> Arrayww;

    const mpreal vf = 0.0;
    const mpreal J = -1.0;
    const mpreal Jx = 0.1;

    const int begin = 6;
    const int end = 7;

    std::string hashlabel = "";
    if (CMD_LINE_PARAMS and argc > 4)
        hashlabel = std::string(argv[4]);

    NumMethod::ForLoopParams<mpreal> fparams;
    NumMethod::GetXFor<mpreal> recordx;
    NumMethod::EqualSpaceFor couplingsfor;
    fparams.start = 0.0;
    fparams.end = 1;
    fparams.numPoints = 10;
    
    std::vector<mpreal> gaps;
    std::vector<mpreal> gaps2;
    std::vector<mpreal> genergy;
    std::vector<mpreal> corrsZ, corrsX, corrsXY;

    if (CMD_LINE_PARAMS)
        fparams = NumMethod::get_for_from_cmd<mpreal>(argv);

    std::vector<mpreal> fs;
    couplingsfor.loop(recordx, fparams);
    fs = recordx.get_x();
    recordx.clear();

   
    auto body = [&](int width, int i) {

        int ignoremask = ~((1 << (2*width - 1)) - 1);

        Matrixww HE(pows2(2*width-1), pows2(2*width-1));
        Matrixww HO(pows2(2*width-1), pows2(2*width-1));

        Arrayw sigz(pows2(2*width-1));
        Arrayw sigz2(pows2(2*width-1));
        Arrayw sigzn(pows2(2*width-1));
        const int nbulk = 2*width-1;
        for (int i = 0; i < sigz.size(); i++) {
            sigz[i] = i % 2 == 0 ? 1 : -1;
            sigz2[i] = i % 4 < 2 ? 1 : -1;
            sigzn[i] = (2 * ((i & pows2(nbulk)) >> nbulk) - 1);
        }
        Arrayw sigz2zn = sigzn*sigz2;
        auto couplingsbody = [&](mpreal vf, int j) {
   
            HE = Matrixww::Zero(pows2(nbulk), pows2(nbulk));
            HO = Matrixww::Zero(pows2(nbulk), pows2(nbulk));

            //std::cout<< to_string(sigma_z_j(0, 0, 0, pows2)*((0+1) < width));

            for (ulong i = 0; i < pows2(nbulk); i++) {
                for (ulong j = 0; j < pows2(nbulk); j++) {
                    //Assumes width is multiple of 2
                    for (ulong jsite = 0; jsite < 2*width; jsite+=4) { 
                        HE(i, j) += 
                                - J * sigma_z_j_z_m(i, j, jsite, jsite + 1, pows2)*(((jsite + 1) < 2*width))
                                - J * sigma_z_j_z_m(i, j, jsite+2, jsite + 3, pows2)*(((jsite + 3) < 2*width))
                                - Jx * sigma_x_j_x_m(i, j, jsite, jsite + 2, pows2, 2*width, ignoremask)*(((jsite + 2) < 2*width))
                                - Jx * vf * sigma_y_j_y_m(i, j, jsite, jsite + 2, pows2, 2*width, ignoremask)*((jsite + 2) < 2*width)
                                - Jx * vf * sigma_x_j_x_m(i, j, jsite+1, jsite + 3, pows2, 2*width, ignoremask)*(((jsite + 3) < 2*width))
                                - Jx * sigma_y_j_y_m(i, j, jsite+1, jsite + 3, pows2, 2*width, ignoremask)*((jsite + 3) < 2*width)
                                - Jx * vf * sigma_x_j_x_m(i, j, jsite+2, jsite + 4, pows2, 2*width, ignoremask)*(((jsite + 4) < 2*width))
                                - Jx * sigma_y_j_y_m(i, j, jsite+2, jsite + 4, pows2, 2*width, ignoremask)*((jsite + 4) < 2*width)
                                - Jx * sigma_x_j_x_m(i, j, jsite+3, jsite + 5, pows2, 2*width, ignoremask)*(((jsite + 5) < 2*width))
                                - Jx * vf * sigma_y_j_y_m(i, j, jsite+3, jsite + 5, pows2, 2*width, ignoremask)*((jsite + 5) < 2*width);
                               
                         HO(i, j) += 
                                - J * sigma_z_j_z_m(~i, ~j, jsite, jsite + 1, pows2)*(((jsite + 1) < 2*width))
                                - J * sigma_z_j_z_m(~i, ~j, jsite+2, jsite + 3, pows2)*(((jsite + 3) < 2*width))
                                - Jx * sigma_x_j_x_m(~i, ~j, jsite, jsite + 2, pows2, 2*width, ignoremask)*(((jsite + 2) < 2*width))
                                - Jx * vf * sigma_y_j_y_m(~i, ~j, jsite, jsite + 2, pows2, 2*width, ignoremask)*((jsite + 2) < 2*width)
                                - Jx * vf * sigma_x_j_x_m(~i, ~j, jsite+1, jsite + 3, pows2, 2*width, ignoremask)*(((jsite + 3) < 2*width))
                                - Jx * sigma_y_j_y_m(~i, ~j, jsite+1, jsite + 3, pows2, 2*width, ignoremask)*((jsite + 3) < 2*width)
                                - Jx * vf * sigma_x_j_x_m(~i, ~j, jsite+2, jsite + 4, pows2, 2*width, ignoremask)*(((jsite + 4) < 2*width))
                                - Jx * sigma_y_j_y_m(~i, ~j, jsite+2, jsite + 4, pows2, 2*width, ignoremask)*((jsite + 4) < 2*width)
                                - Jx * sigma_x_j_x_m(~i, ~j, jsite+3, jsite + 5, pows2, 2*width, ignoremask)*(((jsite + 5) < 2*width))
                                - Jx * vf *sigma_y_j_y_m(~i, ~j, jsite+3, jsite + 5, pows2, 2*width, ignoremask)*((jsite + 5) < 2*width);
                    }
                    if (PERIODIC) { 
                        HE(i, j) += 
                                - Jx * vf * sigma_x_j_x_m(i, j, 2*width-2, (ulong) 0, pows2, 2*width, ignoremask)
                                - Jx  * sigma_y_j_y_m(i, j, 2*width-2, (ulong) 0, pows2, 2*width, ignoremask)
                                - Jx  * sigma_x_j_x_m(i, j, 2*width-1, (ulong) 1, pows2, 2*width, ignoremask)
                                - Jx * vf * sigma_y_j_y_m(i, j, 2*width-1, (ulong) 1, pows2, 2*width, ignoremask);
                               
                         HO(i, j) += 
                                - Jx * vf * sigma_x_j_x_m(~i, ~j, 2*width-2, 0, pows2, 2*width, ignoremask)
                                - Jx * sigma_y_j_y_m(~i, ~j, 2*width-2, 0, pows2, 2*width, ignoremask)
                                - Jx * sigma_x_j_x_m(~i, ~j, 2*width-1, 1, pows2, 2*width, ignoremask)
                                - Jx * vf * sigma_y_j_y_m(~i, ~j, 2*width-1, 1, pows2, 2*width, ignoremask);
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
            
            std::string label = "Kitaev_";
            if (PERIODIC){
                label+= "periodic_";
            }
            label+= "L_" + to_string(width) + "_vf_" + to_string(vf) + "_Jx_" + to_string(Jx)+ "_Jz_"+to_string(J);

            
            std::ofstream outEs;
            if (WRITE_ENERGIES){
                outEs.open((label + "_energies").c_str(), std::ios::trunc);
                for (int ispec = 0; ispec < eigsE.size(); ispec++){
                    outEs << eigsE[ispec] << '\n' << eigsO[ispec] << '\n';
                }
                outEs.flush();
                outEs.close();
            }
            Matrixww zz(pows2(2*width-1), pows2(2*width-1));
            Matrixww xx(pows2(2*width-1), pows2(2*width-1));
            Matrixww xy_rung(pows2(2*width-1), pows2(2*width-1));
            
            zz = Matrixww::Zero(pows2(nbulk), pows2(nbulk));
            xx = Matrixww::Zero(pows2(nbulk), pows2(nbulk));
            xy_rung = Matrixww::Zero(pows2(nbulk), pows2(nbulk));
            
            //Construct operators
            for (ulong i = 0; i < pows2(nbulk); i++) {
                for (ulong j = 0; j < pows2(nbulk); j++) {
                    zz(i, j) += sigma_z_j_z_m(i, j, 0, 2*width-1, pows2);
                    xx(i, j) += sigma_x_j_x_m(i, j,0, 2*width-1, pows2, 2*width, ignoremask);
                    xy_rung(i, j) += sigma_x_j_x_k_y_l_y_m(i, j, 0, 2*width-2,1, 2*width-1, pows2, 2*width, ignoremask);
                }
            }
            genergy.push_back(eigsE[0]);
            gaps.push_back((eigsE[1]-eigsE[0]));
            corrsZ.push_back((evecsE.col(0).adjoint()*zz*evecsE.col(0)).value());
            corrsX.push_back((evecsE.col(0).adjoint()*xx*evecsE.col(0)).value());
            corrsXY.push_back((evecsE.col(0).adjoint()*xy_rung*evecsE.col(0)).value());
            return false;
        };
        couplingsfor.loop(couplingsbody, fparams);
        std::string label = "Kitaev_";
        if (PERIODIC){
               label+= "periodic_";
        }
            label+= "L_" + to_string(width)
                    + "_Jx_" + to_string(Jx)+ "_Jz_"+to_string(J);
         
    plotter.writeToFile(label + "_genergy", fs, genergy);
    plotter.writeToFile(label + "_gaps", fs, gaps);
    plotter.writeToFile(label + "_corrsZ", fs, corrsZ);
    plotter.writeToFile(label + "_corrsX", fs, corrsX);
    plotter.writeToFile(label + "_corrsXY_rung", fs, corrsXY);
        return false;
    };
    
    NumMethod::Range forloop;
    forloop.loop(body, begin, end, 2);
   
    std::vector<int> widths = forloop.get_x(begin, end);

    return 0;

}
#endif
