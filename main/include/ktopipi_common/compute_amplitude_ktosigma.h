#ifndef _FIT_KTOPIPI_GPARITY_COMPUTE_AMPLITUDE_KTOSIGMA_H
#define _FIT_KTOPIPI_GPARITY_COMPUTE_AMPLITUDE_KTOSIGMA_H

#include<config.h>
#include<utils/macros.h>

#include "compute_amplitude.h"

CPSFIT_START_NAMESPACE

template<typename Controls>
typename Controls::outputType computeKtoSigmaAmplitudeType1_2(const typename Controls::inputType &in){
  Controls controls(in);

  static std::map<int,int> idx_map {  {1,0}, {6,1}, {8,2}, {11,3}, {19,4}  };

  for(int t=0;t<controls.size();t++){
#define A(I) controls.A(I-1,t)
#define D(IDX,I,G1,J,G2) controls.C(idx_map[IDX],I,G1,J,G2,t)

    A(1) = 0.5*D(6, 0, VmA, 1, VpA) -0.5*D(1, 0, VmA, 1, VpA);
    A(2) = 0.5*D(11, 0, VmA, 1, VpA) -0.5*D(8, 0, VmA, 1, VpA);
    A(3) = 0.5*D(6, 0, VmA, 1, VpA) + 0.5*D(6, 0, VmA, 0, VmA) -0.5*D(1, 0, VmA, 1, VpA) -0.5*D(1, 0, VmA, 0, VmA);
    A(4) = D(11, 0, VmA, 1, VpA) -0.5*D(8, 0, VmA, 1, VpA) -0.5*D(19, 0, VmA, 0, VmA);
    A(5) = D(6, 0, VmA, 1, VmA) -0.5*D(1, 0, VmA, 1, VmA) -0.5*D(1, 0, VmA, 0, VpA);
    A(6) = D(11, 0, VmA, 1, VmA) -0.5*D(8, 0, VmA, 1, VmA) -0.5*D(19, 0, VmA, 0, VpA);
    A(7) = 0.25*D(6, 0, VmA, 1, VmA) -0.5*D(1, 0, VmA, 1, VmA) + 0.25*D(1, 0, VmA, 0, VpA);
    A(8) = 0.25*D(11, 0, VmA, 1, VmA) -0.5*D(8, 0, VmA, 1, VmA) + 0.25*D(19, 0, VmA, 0, VpA);
    A(9) = 0.25*D(6, 0, VmA, 1, VpA) -0.5*D(1, 0, VmA, 1, VpA)+ 0.25*D(1, 0, VmA, 0, VmA);
    A(10) = 0.25*D(11, 0, VmA, 1, VpA) -0.5*D(8, 0, VmA, 1, VpA) + 0.25*D(19, 0, VmA, 0, VmA);
                            
#undef D
#undef A

  }
  return controls.out;
}


template<typename Controls>
typename Controls::outputType computeKtoSigmaAmplitudeType3(const typename Controls::inputType &in){
  Controls controls(in);

  static std::map<int,int> idx_map {  {2,0}, {3,1}, {7,2}, {10,3}, {14,4}, {16,5}, {18,6}, {21,7}, {23,8} };

  for(int t=0;t<controls.size();t++){
#define A(I) controls.A(I-1,t)
#define D(IDX,I,G1,J,G2) controls.C(idx_map[IDX],I,G1,J,G2,t)

    A(1) = 0.5*D(2, 0, VmA, 1, VpA) -0.5*D(3, 0, VmA, 1, VpA);
    A(2) = 0.5*D(10, 0, VmA, 1, VpA) -0.5*D(7, 0, VmA, 1, VpA);
    A(3) = 0.5*D(2, 0, VmA, 1, VpA) + 0.5*D(2, 0, VmA, 0, VmA) -0.5*D(3, 0, VmA, 1, VpA) 
      -0.5*D(3, 0, VmA, 0, VmA) + 0.5*D(14, 0, VmA, 0, VmA) -0.5*D(16, 0, VmA, 0, VmA);
    A(4) = D(10, 0, VmA, 1, VpA) -0.5*D(7, 0, VmA, 1, VpA) -0.5*D(18, 0, VmA, 0, VmA) 
      + 0.5*D(21, 0, VmA, 0, VmA) -0.5*D(23, 0, VmA, 0, VmA);
    A(5) = D(2, 0, VmA, 1, VmA) -0.5*D(3, 0, VmA, 1, VmA) -0.5*D(3, 0, VmA, 0, VpA) 
      + 0.5*D(14, 0, VmA, 0, VpA) -0.5*D(16, 0, VmA, 0, VpA);
    A(6) = D(10, 0, VmA, 1, VmA) -0.5*D(7, 0, VmA, 1, VmA) -0.5*D(18, 0, VmA, 0, VpA)
      + 0.5*D(21, 0, VmA, 0, VpA) -0.5*D(23, 0, VmA, 0, VpA);
    A(7) = 0.25*D(2, 0, VmA, 1, VmA) -0.5*D(3, 0, VmA, 1, VmA) + 0.25*D(3, 0, VmA, 0, VpA) 
      -0.25*D(14, 0, VmA, 0, VpA) + 0.25*D(16, 0, VmA, 0, VpA);
    A(8) = 0.25*D(10, 0, VmA, 1, VmA) -0.5*D(7, 0, VmA, 1, VmA) + 0.25*D(18, 0, VmA, 0, VpA)
      -0.25*D(21, 0, VmA, 0, VpA) + 0.25*D(23, 0, VmA, 0, VpA);
    A(9) = 0.25*D(2, 0, VmA, 1, VpA) -0.5*D(3, 0, VmA, 1, VpA) + 0.25*D(3, 0, VmA, 0, VmA)
      -0.25*D(14, 0, VmA, 0, VmA) + 0.25*D(16, 0, VmA, 0, VmA);
    A(10) = 0.25*D(10, 0, VmA, 1, VpA) -0.5*D(7, 0, VmA, 1, VpA) + 0.25*D(18, 0, VmA, 0, VmA)
      -0.25*D(21, 0, VmA, 0, VmA) + 0.25*D(23, 0, VmA, 0, VmA);
                            
#undef D
#undef A

  }
  return controls.out;
}


template<typename Controls>
typename Controls::outputType computeKtoSigmaAmplitudeType4(const typename Controls::inputType &in){
  Controls controls(in);

  static std::map<int,int> idx_map {  {1,0}, {6,1}, {8,2}, {11,3}, {19,4}  };

  for(int t=0;t<controls.size();t++){
#define A(I) controls.A(I-1,t)
#define D(IDX,I,G1,J,G2) controls.C(idx_map[IDX],I,G1,J,G2,t)

    A(1) = -0.5*D(5, 0, VmA, 1, VpA) + 0.5*D(4, 0, VmA, 1, VpA);
    A(2) = -0.5*D(12, 0, VmA, 1, VpA) + 0.5*D(9, 0, VmA, 1, VpA);
    A(3) = -0.5*D(5, 0, VmA, 1, VpA) -0.5*D(5, 0, VmA, 0, VmA) + 0.5*D(4, 0, VmA, 1, VpA)
      + 0.5*D(4, 0, VmA, 0, VmA) -0.5*D(13, 0, VmA, 0, VmA) + 0.5*D(15, 0, VmA, 0, VmA);
    A(4) = -1*D(12, 0, VmA, 1, VpA) + 0.5*D(9, 0, VmA, 1, VpA) + 0.5*D(17, 0, VmA, 0, VmA)
      -0.5*D(20, 0, VmA, 0, VmA) + 0.5*D(22, 0, VmA, 0, VmA);
    A(5) = -1*D(5, 0, VmA, 1, VmA) + 0.5*D(4, 0, VmA, 1, VmA) + 0.5*D(4, 0, VmA, 0, VpA)
      -0.5*D(13, 0, VmA, 0, VpA) + 0.5*D(15, 0, VmA, 0, VpA);
    A(6) = -1*D(12, 0, VmA, 1, VmA) + 0.5*D(9, 0, VmA, 1, VmA) + 0.5*D(17, 0, VmA, 0, VpA)
      -0.5*D(20, 0, VmA, 0, VpA) + 0.5*D(22, 0, VmA, 0, VpA);
    A(7) = -0.25*D(5, 0, VmA, 1, VmA) + 0.5*D(4, 0, VmA, 1, VmA) -0.25*D(4, 0, VmA, 0, VpA)
      + 0.25*D(13, 0, VmA, 0, VpA) -0.25*D(15, 0, VmA, 0, VpA);
    A(8) = -0.25*D(12, 0, VmA, 1, VmA) + 0.5*D(9, 0, VmA, 1, VmA) -0.25*D(17, 0, VmA, 0, VpA)
      + 0.25*D(20, 0, VmA, 0, VpA) -0.25*D(22, 0, VmA, 0, VpA);
    A(9) = -0.25*D(5, 0, VmA, 1, VpA) + 0.5*D(4, 0, VmA, 1, VpA) -0.25*D(4, 0, VmA, 0, VmA)
      + 0.25*D(13, 0, VmA, 0, VmA) -0.25*D(15, 0, VmA, 0, VmA);
    A(10) = -0.25*D(12, 0, VmA, 1, VpA) + 0.5*D(9, 0, VmA, 1, VpA) -0.25*D(17, 0, VmA, 0, VmA)
      + 0.25*D(20, 0, VmA, 0, VmA) -0.25*D(22, 0, VmA, 0, VmA);

#undef D
#undef A

  }
  return controls.out;
}


//We collectively use index 2 for the type1/2 K->sigma contractions as it makes looping indices 2,3,4 easier
template<typename Controls>
typename Controls::outputType computeKtoSigmaAmplitudeType(const int i, const typename Controls::inputType &in){
  switch(i){
  case 2: 
    return computeKtoSigmaAmplitudeType1_2<Controls>(in);
  case 3: 
    return computeKtoSigmaAmplitudeType3<Controls>(in);
  case 4: 
    return computeKtoSigmaAmplitudeType4<Controls>(in);
  default:
    assert(0);
  }
};

CPSFIT_END_NAMESPACE

#endif
