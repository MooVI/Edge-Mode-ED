/*
 * File:   main.cpp
 * Author: Jack Kemp
 *
 * Created on 13 October 2014, 11:54
 */

#if 1

#include<Eigen/Dense>
//#include<Eigen/MPRealSupport>
#include<NumericalMethods/NumericalMethods/Random.h>
#include<NumericalMethods/NumericalMethods/Statistics.h>
#include<NumericalMethods/NumericalMethods/SamplingForLoops.h>
#include<NumericalMethods/NumericalMethods/Mod.h>
#include<NumericalMethods/NumericalMethods/Useful.h>
#include<NumericalMethods/NumericalMethods/1DRootFinding.h>
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

  const bool FAKE_BCS = false;
  const bool FAKE_END_GAMMA = false;
  const bool LONG_RANGE = false;

  const bool X1_Z2_EDGE = false;
  const bool Y1_Y2_Z3_EDGE = false;
  const bool X1_EDGE = false;
  const bool Y1_EDGE = false;
  const bool Y1_Z2_EDGE = false;

  const bool WRITE_ENERGIES = false;
  const bool WRITE_OVERLAPS = false;
  const bool WRITE_BULK = false;
  const bool WRITE_MEANS = false;
  const bool WRITE_VARS = false;
  const bool WRITE_MAX_OVERLAPS = false;
  const bool WRITE_PAIRED_EDIFFS = false;
  const bool WRITE_ALL_DECAY = true;
  const bool WRITE_DECAY_TIME = false;
  const bool WRITE_PAIRED_DECAY = false;
  const bool WRITE_ALL_LEVEL_SPACINGS = false;
  const bool WRITE_MEAN_LEVEL_SPACINGS = false;
  const bool WRITE_CORR = false;
  const bool WRITE_FINITE_T_CORR = false;
  const bool WRITE_BULK_DECAY = false;
  const bool WRITE_FINITE_T_BULK_DECAY = false;
  const bool WRITE_BULK_STATE_DECAY = false;
  const bool WRITE_ALL_SPIN_DECAYS = false;

  const bool WRITE_FINITE_T_DECAY = false;
  const bool WRITE_STATE_DECAY = false;

  const bool FINITE_TEMPERATURE = false;
  const bool FINITE_TEMPERATURE_LOOP = false;
  const bool ENERGY_WINDOW_LOOP = false;


  const bool CMD_LINE_PARAMS = true;


  typedef Eigen::Matrix<mpreal, Eigen::Dynamic, Eigen::Dynamic> Matrixww;
  typedef Eigen::Matrix<mpreal, Eigen::Dynamic, 1> Vectorw;
  typedef Eigen::Array<mpreal, Eigen::Dynamic, 1> Arrayw;
  typedef Eigen::Array<mpreal, Eigen::Dynamic, Eigen::Dynamic> Arrayww;

  const mpreal r = 1.0;
  const mpreal theta = 0.0;
  const int staggered = 0;
  const mpreal mag = 0.0;
  const mpreal J2 = 0.00;
  const mpreal J3 = 0.0;
  const mpreal J4 = 0.0;
  const mpreal Jy = 0.0;
  const mpreal J = 1.0;
  const mpreal f = 0.25;
  const mpreal f_fake_end = 10;
  const mpreal V = 0.0;
  const mpreal alpha = 0;
  const mpreal Js [] = {J3, J4};

  const mpreal T = 0.1;
  const mpreal Ew = 0.0;

  NumMethod::RunningStats<mpfr::mpreal> statoverlap, stateigdiff;
  NumMethod::WeightedRunningStats<mpfr::mpreal> wstatoverlap, wstateigdiff;

  const int begin = 8;
  const int end = 15;

  std::string hashlabel = "";
  if (CMD_LINE_PARAMS and argc > 4)
    hashlabel = std::string(argv[4]);

  std::string qtype = "";
      if (X1_Z2_EDGE){
        qtype = "_x1z2";
      }
      else if (Y1_Z2_EDGE){
          qtype = "_y1z2";
      }
      else if (Y1_EDGE){
          qtype = "_y1";
      }
      else if (X1_EDGE){
        qtype = "_x1";
      }
      else if (Y1_Y2_Z3_EDGE){
        qtype = "_y1y2z3";
      }

  NumMethod::ForLoopParams<mpreal> fparams;
  NumMethod::GetXFor<mpreal> recordx;
  NumMethod::EqualSpaceFor couplingsfor;
  fparams.start = 0.0;
  fparams.end = 1.0;
  fparams.numPoints = 200;

  if (CMD_LINE_PARAMS)
    fparams = NumMethod::get_for_from_cmd<mpreal>(argv);


  std::vector<mpfr::mpreal> maxoverlap, eigdiffs, varoverlaps;
  std::vector<mpfr::mpreal> vareigdiffs, meanspacings, correlations;
  std::vector<mpfr::mpreal> decay_times;
  
  std::vector<mpreal> fs;
  couplingsfor.loop(recordx, fparams);
  fs = recordx.get_x();
  recordx.clear();

  NumMethod::Brent rootfinder;
  mpreal root_min_time = 1;
  mpreal root_max_time = 1e9;
  mpreal root_error = 10;

  const int tstate = 0;
  NumMethod::EqualSpaceFor tfor;
  NumMethod::ForLoopParams<mpreal> tparams;
  std::vector<mpreal> ts;
  std::vector<mpreal> tends = {1e6,1e7,1e8,1e9};
  if (WRITE_ALL_DECAY | WRITE_PAIRED_DECAY | WRITE_FINITE_T_DECAY | WRITE_STATE_DECAY | WRITE_BULK_DECAY | WRITE_ALL_SPIN_DECAYS) {
    tparams.start = 0.0;
    //tparams.end = 1e6;
    tparams.numPoints = 5001;
    tfor.loop(recordx, tparams);
    ts = recordx.get_x();
    recordx.clear();
  }
    
  NumMethod::EqualSpaceFor Tfor;
  NumMethod::ForLoopParams<mpreal> Tparams;
  std::vector<mpreal> Ts;
  if (FINITE_TEMPERATURE_LOOP | WRITE_FINITE_T_DECAY | WRITE_FINITE_T_CORR) {
    Tparams.start = 0.25;
    Tparams.end = 7.0;
    Tparams.numPoints = 28;
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

    tparams.end = tends[i];
    int tstate = pows2(width-1)-4; //NOTICE THIS!!!!!

    Arrayw moverlaps(pows2(width - 1)), meigdiff(pows2(width - 1));

    int ignoremask = ~((1 << (width - 1)) - 1);

    Matrixww HE(pows2(width - 1), pows2(width - 1));
    Matrixww HO(pows2(width - 1), pows2(width - 1));

    Arrayw sigz(pows2(width - 1));
    Arrayw sigz2(pows2(width - 1));
    Arrayw sigz3(pows2(width - 1));
    for (int i = 0; i < sigz.size(); i++) {
      sigz[i] = i % 2 == 0 ? 1 : -1;
      sigz2[i] = i % 4 < 2 ? 1 : -1;
      sigz3[i] = i % 8 < 4 ? 1 : -1;
    }

    auto couplingsbody = [&](mpreal f, int j) {
      // mpreal f = mag*cos(theta);
      //mpreal V = mag*sin(theta);
      //mpreal f = 1*V; //XXXXXX notice this!
      //                     mpreal f = J2;
      mpreal beta = 1.0/T;
      const mpreal Js [] = {J3, J4};
      std::vector<mpreal> J2s(width);
      tparams.end = 10*3.14/((1-f*f)*pow(f,width));// NOTICE THIS!
      if (WRITE_ALL_DECAY | WRITE_PAIRED_DECAY | WRITE_FINITE_T_DECAY | WRITE_STATE_DECAY | WRITE_BULK_DECAY | WRITE_ALL_SPIN_DECAYS) {
        tfor.loop(recordx, tparams);
        ts = recordx.get_x();
        recordx.clear();
      }

      for (int i = 0; i < width; ++i)
	J2s[i] = J2;
      J2s[staggered] *= r;
      HE = Matrixww::Zero(pows2(width - 1), pows2(width - 1));
      HO = Matrixww::Zero(pows2(width - 1), pows2(width - 1));

      //std::cout<< to_string(sigma_z_j(0, 0, 0, pows2)*((0+1) < width));
      if (FAKE_BCS){
	for (ulong i = 0; i < pows2(width - 1); i++) {
	  for (ulong j = 0; j < pows2(width - 1); j++) {
	    for (ulong jsite = 0; jsite < width; jsite++) {
	      HE(i, j) += -f * sigma_x_j(i, j, jsite, pows2, width, ignoremask)*(((jsite + 1) < width))
		- J * sigma_z_j_z_m(i, j, jsite, jsite + 1, pows2)*(((jsite + 1) < width))
		- J2s[jsite] * sigma_z_j_z_m(i, j, jsite, jsite + 2, pows2)*(((jsite + 2) < width))
		- V * sigma_x_j_x_m(i, j, jsite, jsite + 1, pows2, width, ignoremask)*((jsite + 2) < width)
		- Jy * sigma_y_j_y_m(i, j, jsite, jsite + 1, pows2, width, ignoremask)*((jsite + 2) < width)
		- Js[jsite % 2] * sigma_z_p_x_j_z_m(i, j, jsite, jsite + 1, jsite + 2, pows2, width, ignoremask)*(((jsite + 2) < width));
	      HO(i, j) += -f * sigma_x_j(~i, ~j, jsite, pows2, width, ignoremask)*(((jsite + 1) < width))
		- J * sigma_z_j_z_m(~i, ~j, jsite, jsite + 1, pows2)*(((jsite + 1) < width))
		- J2s[jsite] * sigma_z_j_z_m(~i, ~j, jsite, jsite + 2, pows2)*(((jsite + 2) < width))
		- V * sigma_x_j_x_m(~i, ~j, jsite, jsite + 1, pows2, width, ignoremask)*((jsite + 2) < width)
		- Jy * sigma_y_j_y_m(~i, ~j, jsite, jsite + 1, pows2, width, ignoremask)*((jsite + 2) < width)
		- Js[jsite % 2] * sigma_z_p_x_j_z_m(~i, ~j, jsite, jsite + 1, jsite + 2, pows2, width, ignoremask)*(((jsite + 2) < width));
	    }
	  }
	}
      }
      else if(FAKE_END_GAMMA) {
	for (ulong i = 0; i < pows2(width - 1); i++) {
	  for (ulong j = 0; j < pows2(width - 1); j++) {
	    for (ulong jsite = 0; jsite < width; jsite++) {
	      HE(i, j) += -f * sigma_x_j(i, j, jsite, pows2, width, ignoremask)
	      - J * sigma_z_j_z_m(i, j, jsite, jsite + 1, pows2)*(((jsite + 1) < width))
	      - J2s[jsite] * sigma_z_j_z_m(i, j, jsite, jsite + 2, pows2)*(((jsite + 2) < width))
	      - V * sigma_x_j_x_m(i, j, jsite, jsite + 1, pows2, width, ignoremask)*((jsite + 1) < width)
	      - Jy * sigma_y_j_y_m(i, j, jsite, jsite + 1, pows2, width, ignoremask)*((jsite + 2) < width)
	      - Js[jsite % 2] * sigma_z_p_x_j_z_m(i, j, jsite, jsite + 1, jsite + 2, pows2, width, ignoremask)*(((jsite + 2) < width));
	      
	      HO(i, j) += -f * sigma_x_j(~i, ~j, jsite, pows2, width, ignoremask)
	      - J * sigma_z_j_z_m(~i, ~j, jsite, jsite + 1, pows2)*(((jsite + 1) < width))
	      - J2s[jsite] * sigma_z_j_z_m(~i, ~j, jsite, jsite + 2, pows2)*(((jsite + 2) < width))
	      - V * sigma_x_j_x_m(~i, ~j, jsite, jsite + 1, pows2, width, ignoremask)*((jsite + 1) < width)
	      - Jy * sigma_y_j_y_m(~i, ~j, jsite, jsite + 1, pows2, width, ignoremask)*((jsite + 2) < width)
	      - Js[jsite % 2] * sigma_z_p_x_j_z_m(~i, ~j, jsite, jsite + 1, jsite + 2, pows2, width, ignoremask)*(((jsite + 2) < width));
	     
	    }
	    HO(i, j) += -f_fake_end * sigma_x_j(~i, ~j, width-1, pows2, width, ignoremask);
	    HE(i, j) += -f_fake_end * sigma_x_j(i, j, width-1, pows2, width, ignoremask);
	  }
	}
      }
      else if (LONG_RANGE){
	for (ulong i = 0; i < pows2(width - 1); i++) {
	  for (ulong j = 0; j < pows2(width - 1); j++) {
	    for (ulong jsite = 0; jsite < width; jsite++) {
	      HE(i, j) += -f * sigma_x_j(i, j, jsite, pows2, width, ignoremask)
	      - V * sigma_x_j_x_m(i, j, jsite, jsite + 1, pows2, width, ignoremask)*((jsite + 1) < width);
	      HO(i, j) += -f * sigma_x_j(~i, ~j, jsite, pows2, width, ignoremask)
	      - V * sigma_x_j_x_m(~i, ~j, jsite, jsite + 1, pows2, width, ignoremask)*((jsite + 1) < width);
	      for (ulong r = width-1-jsite; r > 0; r--){
		HE(i, j) += - J * pow((mpreal) r, -alpha)*sigma_z_j_z_m(i, j, jsite, jsite + r, pows2);
		HO(i, j) += - J * pow((mpreal) r, -alpha)*sigma_z_j_z_m(~i, ~j, jsite, jsite + r, pows2);
	      }
	    }
	  }
	}
      }
      else {
	for (ulong i = 0; i < pows2(width - 1); i++) {
	  for (ulong j = 0; j < pows2(width - 1); j++) {
	    for (ulong jsite = 0; jsite < width; jsite++) {
	      HE(i, j) += -f * sigma_x_j(i, j, jsite, pows2, width, ignoremask)
	      - J * sigma_z_j_z_m(i, j, jsite, jsite + 1, pows2)*(((jsite + 1) < width))
	      - J2s[jsite] * sigma_z_j_z_m(i, j, jsite, jsite + 2, pows2)*(((jsite + 2) < width))
	      - V * sigma_x_j_x_m(i, j, jsite, jsite + 1, pows2, width, ignoremask)*((jsite + 1) < width)
	      - Jy * sigma_y_j_y_m(i, j, jsite, jsite + 1, pows2, width, ignoremask)*((jsite + 2) < width)
	      - Js[jsite % 2] * sigma_z_p_x_j_z_m(i, j, jsite, jsite + 1, jsite + 2, pows2, width, ignoremask)*(((jsite + 2) < width));
	      HO(i, j) += -f * sigma_x_j(~i, ~j, jsite, pows2, width, ignoremask)
	      - J * sigma_z_j_z_m(~i, ~j, jsite, jsite + 1, pows2)*(((jsite + 1) < width))
	      - J2s[jsite] * sigma_z_j_z_m(~i, ~j, jsite, jsite + 2, pows2)*(((jsite + 2) < width))
	      - V * sigma_x_j_x_m(~i, ~j, jsite, jsite + 1, pows2, width, ignoremask)*((jsite + 1) < width)
	      - Jy * sigma_y_j_y_m(~i, ~j, jsite, jsite + 1, pows2, width, ignoremask)*((jsite + 2) < width)
	      - Js[jsite % 2] * sigma_z_p_x_j_z_m(~i, ~j, jsite, jsite + 1, jsite + 2, pows2, width, ignoremask)*(((jsite + 2) < width));
	    }
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
      std::string label = "Ising_equalspacefreq"+qtype+"_L_" + to_string(width) + "_f_" + to_string(f)
      + "_V_" + to_string(V)+ "_J2_" + to_string(J2) + "_J_" + to_string(J); 
      std::ofstream outEs;
      if (WRITE_ENERGIES)
	outEs.open((label + "_energies").c_str(), std::ios::trunc);
      if (X1_Z2_EDGE or Y1_Z2_EDGE or X1_EDGE or Y1_EDGE or Y1_Y2_Z3_EDGE){
	for (int i = 0; i < pows2(width-2); i++) {
	  evecsO.row(2*i).swap(evecsO.row(2*i+1));
	}
      }
      if (Y1_Y2_Z3_EDGE){
	for (int i = 0; (4*i+3) < pows2(width-1); i++) {
	  evecsO.row(4*i).swap(evecsO.row(4*i+2));
          evecsO.row(4*i+1).swap(evecsO.row(4*i+3));
	}
      }
      if (WRITE_OVERLAPS or WRITE_ALL_DECAY or WRITE_FINITE_T_DECAY
	  or WRITE_STATE_DECAY or WRITE_DECAY_TIME) {
	Arrayww overlap = Matrixww::Zero(pows2(width - 1), pows2(width - 1));
	for (int ispec = 0; ispec < pows2(width - 1); ispec++) {
	  auto spec = evecsE.col(ispec).array();
	  if (not (X1_Z2_EDGE or Y1_Z2_EDGE or X1_EDGE or Y1_EDGE or Y1_Y2_Z3_EDGE)){
	    for (int i = 0; i <  pows2(width-1); i++) {
	      overlap(i, ispec) = (spec * sigz * evecsO.col(i).array()).sum();
	    }
	  }
	  else if (X1_Z2_EDGE){
	    for (int i = 0; i < pows2(width-1); i++) {
	      overlap(i, ispec) = (spec * sigz2 * evecsO.col(i).array()).sum();
	    }
	  }
          else if (X1_EDGE){
            for (int i = 0; i < pows2(width-1); i++) {
	      overlap(i, ispec) = (spec * evecsO.col(i).array()).sum();
            }
          }
          else if(Y1_Z2_EDGE) {
	    for (int i = 0; i < pows2(width-1); i++) {
	      overlap(i, ispec) = (spec * sigz2 * sigz * evecsO.col(i).array()).sum(); //sigz1 gives the (-1¸ 1) part of y1
	    }
	  }
           else if(Y1_EDGE) {
	    for (int i = 0; i < pows2(width-1); i++) {
	      overlap(i, ispec) = (spec * sigz * evecsO.col(i).array()).sum(); //sigz1 gives the (-1¸ 1) part of y1
	    }
	  }
           else if(Y1_Y2_Z3_EDGE) {
	    for (int i = 0; i < pows2(width-1); i++) {
	      overlap(i, ispec) = (-spec * sigz3 * sigz2 * sigz * evecsO.col(i).array()).sum(); //sigz1 gives the (-1¸ 1) part of y1
	    }
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
	if (WRITE_STATE_DECAY){
	  Arrayw toverlaps (tparams.numPoints);
	  Arrayw stateoverlap = overlap.row(tstate).array().square();
	  Arrayw alleigdiff(pows2(width-1));
	  for (ulong i = 0; i < pows2(width - 1); i++) {
	    alleigdiff[i] = eigsE[tstate] - eigsO[i];
	  }
	  auto tfunc = [&](mpreal t, int j) {
	    toverlaps[j] = ((stateoverlap*alleigdiff.unaryExpr([ = ](mpreal x){return cos(x * t);}))).sum();
	    return false;
	  };
	  tfor.loop(tfunc, tparams);
	  plotter.writeToFile(label+"_state_"+to_string(tstate)+"_decay", ts, toverlaps);
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
        if (WRITE_DECAY_TIME) {
	  if (not WRITE_ALL_DECAY){
            overlap = overlap.array().square();
          }
	  Arrayww alleigdiff(pows2(width - 1), pows2(width - 1));
	  for (ulong i = 0; i < pows2(width - 1); i++) {
	    for (ulong j = 0; j < pows2(width - 1); j++) {
	      alleigdiff(i, j) = eigsE[j] - eigsO[i];
	    }
	  }
	  auto tfunc = [&](mpreal t) {
                         return (((overlap.array()
                                   * alleigdiff.unaryExpr(
                                                          [ = ](mpreal x){return cos(x * t);})).sum() / (mpreal) pows2(width - 1)))-exp(-1);
	  };
          decay_times.push_back(rootfinder(tfunc, root_min_time, root_max_time, root_error));
          std::cout<<"TIME:"<<std::endl<<decay_times.back()<<std::endl;
	}
	if (WRITE_FINITE_T_DECAY) {
	  Arrayw toverlaps(tparams.numPoints);
	  if (not WRITE_ALL_DECAY){
	    overlap = overlap.array().square();
	  }
	  auto Tfunc = [&](mpreal T, int j) {
	    Arrayww alleigdiff(pows2(width - 1), pows2(width - 1));
	    Arrayw eweights(pows2(width - 1)); 
	    for (ulong i = 0; i < pows2(width - 1); i++) {
	      eweights(i) = exp(-eigsE[i]/T) + exp(-eigsO[i]/T);
	      for (ulong j = 0; j < pows2(width - 1); j++) {
		alleigdiff(i, j) = eigsE[j] - eigsO[i];
	      }
	    }
	    mpreal zdenom = eweights.sum();
	    auto tfunc = [&](mpreal t, int j) {
	      toverlaps[j] = (((overlap.array()
				* alleigdiff.unaryExpr([ = ](mpreal x){return cos(x * t);})).colwise().sum().transpose()
			       *eweights
			       ).sum())
	      / (mpreal) zdenom;
	      return false;
	    };
	    tfor.loop(tfunc, tparams);
	    plotter.writeToFile(label +"_T_"+to_string(T)+ "_long_decay", ts, toverlaps);
	    return false;
	  };
	  Tfor.loop(Tfunc, Tparams);
		  
	}
      } else {
	Arrayw overlap = Vectorw::Zero(pows2(width - 1));
	for (int ispec = 0; ispec < pows2(width - 1); ispec++) {
	  auto spec = evecsE.col(ispec).array();
	  if (not (X1_Z2_EDGE or Y1_Z2_EDGE)){
	    for (int i = 0; i <  pows2(width-1); i++) {
	      overlap(i) = (spec * sigz * evecsO.col(i).array()).sum();
	    }
	  }
	  else {
	    for (int i = 0; i < pows2(width-1); i++) {
	      overlap(i) = (spec * sigz2 * evecsO.col(i).array()).sum();
	    }
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
	//Arrayw energies(pows2(width)-1);
	//for (int i=0;i < pows2(width-1); i++){
	//energies[2*i] = eigsE[i];
	//energies[2*i+1] = eigsO[i];
	//}
	//std::sort(energies.data(), energies.data()+energies.size());
	//Arrayw spacings = (energies.head(pows2(width)-1)-energies.tail(pows2(width)-1)).abs();
	Arrayw spacings = (eigsE.array().head(pows2(width-1)-1)-eigsE.array().tail(pows2(width-1)-1)).abs();
	Arrayw normspacings (pows2(width-1)-2);
	for (int i=0; i< pows2(width-1)-2; i++){
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
      Arrayww overlap;
      if (WRITE_BULK_DECAY or WRITE_FINITE_T_BULK_DECAY or WRITE_BULK_STATE_DECAY){
	overlap = Matrixww::Zero(pows2(width - 1), pows2(width - 1));
	Arrayw sigzmid(pows2(width - 1));
	for (int i = 0; i < sigzmid.size(); i++) {
	  sigzmid[i] = (2 * ((i & pows2(width/2)) >> (width/2)) - 1);
	}
	for (int ispec = 0; ispec < pows2(width - 1); ispec++) {
	  auto spec = evecsE.col(ispec).array();
	  for (int i = 0; i < pows2(width-1); i++) {
	    overlap(i, ispec) = (spec * sigzmid * evecsO.col(i).array()).sum();
	  }
	}
      }
      if (WRITE_BULK_DECAY){
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
	plotter.writeToFile(label + "_bulk_decay", ts, toverlaps);
      }
      if (WRITE_BULK_STATE_DECAY){
	  Arrayw toverlaps (tparams.numPoints);
	  Arrayw stateoverlap = overlap.row(tstate).array().square();
	  Arrayw alleigdiff(pows2(width-1));
	  for (ulong i = 0; i < pows2(width - 1); i++) {
	    alleigdiff[i] = eigsE[tstate] - eigsO[i];
	  }
	  auto tfunc = [&](mpreal t, int j) {
	    toverlaps[j] = ((stateoverlap*alleigdiff.unaryExpr([ = ](mpreal x){return cos(x * t);}))).sum();
	    return false;
	  };
	  tfor.loop(tfunc, tparams);
	  plotter.writeToFile(label+"_state_"+to_string(tstate)+"_bulk_decay", ts, toverlaps);
      }
      if (WRITE_FINITE_T_BULK_DECAY) {
	//Not the fastest way of doing it (could coalesce with other finite t decay)
	Arrayw toverlaps(tparams.numPoints);
	if (not WRITE_BULK_DECAY){
	  overlap = overlap.array().square();
	}
	auto Tfunc = [&](mpreal T, int j) {
	  Arrayww alleigdiff(pows2(width - 1), pows2(width - 1));
	  Arrayw eweights(pows2(width - 1)); 
	  for (ulong i = 0; i < pows2(width - 1); i++) {
	    eweights(i) = exp(-eigsE[i]/T) + exp(-eigsO[i]/T);
	    for (ulong j = 0; j < pows2(width - 1); j++) {
	      alleigdiff(i, j) = eigsE[j] - eigsO[i];
	    }
	  }
	  mpreal zdenom = eweights.sum();
	  auto tfunc = [&](mpreal t, int j) {
	    toverlaps[j] = (((overlap.array()
			      * alleigdiff.unaryExpr([ = ](mpreal x){return cos(x * t);})).colwise().sum().transpose()
			     *eweights
			     ).sum())
	    / (mpreal) zdenom;
	    return false;
	  };
	  tfor.loop(tfunc, tparams);
	  plotter.writeToFile(label +"_T_"+to_string(T)+ "_bulk_decay", ts, toverlaps);
	  return false;
	};
	Tfor.loop(Tfunc, Tparams);		  
      }
      if (WRITE_ALL_SPIN_DECAYS){
	Arrayww alleigdiff(pows2(width - 1), pows2(width - 1));
	for (ulong i = 0; i < pows2(width - 1); i++) {
	  for (ulong j = 0; j < pows2(width - 1); j++) {
	    alleigdiff(i, j) = eigsE[j] - eigsO[i];
	  }
	}
	Arrayww overlap = Matrixww::Zero(pows2(width - 1), pows2(width - 1));
	Arrayw sigzmid(pows2(width - 1));
	Arrayw toverlaps(tparams.numPoints);
	for (int jsite=0; jsite < (width/2)+1; jsite++){
	  for (int i = 0; i < sigzmid.size(); i++) {
	    sigzmid[i] = (2 * ((i & pows2(jsite)) >> jsite) - 1);
	  }
	  for (int ispec = 0; ispec < pows2(width - 1); ispec++) {
	    auto spec = evecsE.col(ispec).array();
	    for (int i = 0; i < pows2(width-1); i++) {
	      overlap(i, ispec) = (spec * sigzmid * evecsO.col(i).array()).sum();
	    }
	  }
	  overlap = overlap.array().square();
	  auto tfunc = [&](mpreal t, int j) {
	    toverlaps[j] = ((overlap.array() * alleigdiff.unaryExpr([ = ](mpreal x){return cos(x * t);})).sum() / (mpreal) pows2(width - 1));
	    return false;
	  };
	  tfor.loop(tfunc, tparams);
	  plotter.writeToFile(label + "_spin_"+to_string(jsite)+"_decay", ts, toverlaps);
	}
      }
      if (WRITE_CORR){
	Arrayww overlap = Matrixww::Zero(pows2(width - 1), pows2(width - 1));
	Arrayw sigz1xL(pows2(width - 1));
	for (int i = 0; i < sigz1xL.size(); i++) {
	  sigz1xL[i] = (2 * ((i & pows2(width-1)) >> (width-1)) - 1)*sigz[i];
	  }
	for (int ispec = 0; ispec < pows2(width - 1); ispec++) {
	  auto spec = evecsE.col(ispec).array();
	  for (int i = 0; i < pows2(width-1); i++) {
	      overlap(i, ispec) = (spec * sigz1xL * evecsE.col(i).array()).sum();
	    }
	}
	correlations.push_back(overlap.mean());
      }
      if (WRITE_FINITE_T_CORR){
	Arrayww overlap = Matrixww::Zero(pows2(width - 1), pows2(width - 1));
	Arrayw Tcorr(Tparams.numPoints);
	Arrayw sigz1xL(pows2(width - 1));
	for (int i = 0; i < sigz1xL.size(); i++) {
	  sigz1xL[i] = (2 * ((i & pows2(width-1)) >> (width-1)) - 1)*sigz[i];
	}
	auto Tfunc = [&](mpreal T, int j) {
	  Arrayw eweights(pows2(width - 1)); 
	  for (ulong i = 0; i < pows2(width - 1); i++) {
	    eweights(i) = exp(-eigsE[i]/T) + exp(-eigsO[i]/T);
	  }
	  mpreal zdenom = eweights.sum();
	  for (int ispec = 0; ispec < pows2(width - 1); ispec++) {
	    auto spec = evecsE.col(ispec).array();
	    for (int i = 0; i < pows2(width-1); i++) {
	      overlap(i, ispec) = -(spec * sigz1xL * evecsE.col(i).array()).sum()*eweights(i);
	    }
	  }
	  Tcorr[j] = overlap.sum()/zdenom;
	  return false;
	};
	Tfor.loop(Tfunc, Tparams);
	plotter.writeToFile(label +"_finite_T_corr", Ts, Tcorr);
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

    std::string label = hashlabel + "Ising_reson"+qtype+"_L_" + to_string(width) + "_f_" + to_string(f)
      + "_V_" + to_string(V)+ "_J_" + "vary" + "_J2_" + to_string(J); 
    //+ "_f_" + to_string(f);
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
    if (WRITE_CORR){
      plotter.writeToFile(label + "_corr", fs, correlations);
      correlations.clear();
    }
    if (WRITE_DECAY_TIME){
      plotter.writeToFile(label + "_decay_time", fs, decay_times);
      decay_times.clear();
    }
    return false;
  };

  NumMethod::Range forloop;
  forloop.loop(body, begin, end, 2);
  std::vector<int> widths = forloop.get_x(begin, end);

  return 0;

}
#endif

