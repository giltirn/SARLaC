#ifndef _FIT_KTOPIPI_GPARITY_COMPUTE_AMPLITUDE_H
#define _FIT_KTOPIPI_GPARITY_COMPUTE_AMPLITUDE_H

#include<config.h>
#include<utils/macros.h>

#include<distribution.h>
#include<common.h>

#include "data_containers.h"

SARLAC_START_NAMESPACE

enum LR{ VpA, VmA }; 
rawDataDistributionD computeLRcontraction(const int cidx, const int i, const LR g1, const int j, const LR g2, const contractions &from){
  //(V+aA)(V+bA) = V bA + aA V
  double a = g1 == VpA ? 1 : -1;
  double b = g2 == VpA ? 1 : -1;
  return b * from.C(cidx)(i,'V',j,'A') + a * from.C(cidx)(i,'A',j,'V');
}

//For a given type, what is the first diagram index
inline int idxOffset(const int type){
  static const int offs[4] = {1,7,13,23}; 
  return offs[type-1];
}

struct computeAmplitudeAlltKtensorControls{
  struct inputType{
    type1234Data const& data;
    std::vector<int> const &nonzero_tK;
    const int tsep_k_pi;
    inputType(type1234Data const& data, std::vector<int> const &nonzero_tK, const int tsep_k_pi): data(data), nonzero_tK(nonzero_tK), tsep_k_pi(tsep_k_pi){}
  };
  typedef NumericTensor<rawDataDistributionD,3> outputType; //(Qidx,tK_idx,t)

  int ntK;
  int nt;
  outputType out;
  const inputType &in;
  
  computeAmplitudeAlltKtensorControls(const inputType &_in): nt(_in.tsep_k_pi+1), 
							     ntK(_in.nonzero_tK.size()), 
							     in(_in), 
							     out( { 10, int(_in.nonzero_tK.size()), _in.tsep_k_pi+1 } ){}

  inline int size(){ return ntK*nt; }
  
  inline rawDataDistributionD & A(const int I, const int tt){ return out({I,tt/nt, tt%nt}); }
  
  inline rawDataDistributionD C(const int IDX, const int i, const LR g1, const int j, const LR g2, const int tt){
    int tK_idx=tt/nt, t=tt%nt;
    return computeLRcontraction(IDX, i,g1,j,g2, in.data(in.nonzero_tK[tK_idx], t));
  }

  inline void normalize(const double nrm, const int tt){
    int tK_idx=tt/nt, t=tt%nt;
    for(int i=0;i<10;i++) out({i,tK_idx,t}) = out({i,tK_idx,t}) * nrm;
  }
};  



//Note we ignore an overall factor of i

template<typename Controls>
typename Controls::outputType computeAmplitudeType1(const typename Controls::inputType &in){
  Controls controls(in);
  int off = idxOffset(1);
    
  for(int t=0;t<controls.size();t++){
#define A(I) controls.A(I-1,t)
#define C(IDX,I,G1,J,G2) controls.C(IDX-off,I,G1,J,G2,t)

    A(1) = C(1, 1,VpA,0,VmA) - C(4, 0,VmA,1,VpA) -2.*C(4, 0,VmA,0,VmA);
    A(2) = C(2, 1,VpA,0,VmA) - C(5, 0,VmA,1,VpA) -2.*C(6, 0,VmA,0,VmA);
    A(3) = -3.*C(4, 0,VmA,1,VpA) -3.*C(4, 0,VmA,0,VmA);
    A(4) = C(3, 0,VmA,0,VmA) -3.*C(6, 0,VmA,0,VmA) +C(2, 1,VpA,0,VmA) -3.*C(5, 0,VmA,1,VpA);      
    A(5) = -3.*C(4, 0,VmA,1,VmA) -3.*C(4, 0,VmA,0,VpA);
    A(6) = C(3, 0,VpA,0,VmA) -3.*C(6, 0,VmA,0,VpA) + C(2, 1,VmA,0,VmA) -3.*C(5, 0,VmA,1,VmA);
    A(7) = C(1, 1,VmA,0,VmA) -0.5*C(1, 0,VpA,0,VmA) -1.5*C(4, 0,VmA,0,VpA);
    A(8) = -0.5*C(3, 0,VpA,0,VmA) -1.5*C(6, 0,VmA,0,VpA) + C(2, 1,VmA,0,VmA);
    A(9) = C(1, 1,VpA,0,VmA) -0.5*C(1, 0,VmA,0,VmA) -1.5*C(4, 0,VmA,0,VmA);
    A(10) = -0.5*C(3, 0,VmA,0,VmA) -1.5*C(6, 0,VmA,0,VmA) + C(2, 1,VpA,0,VmA);
#undef C
#undef A

    controls.normalize(1./sqrt(6.),t);
  }
  return controls.out;
}

template<typename Controls>
typename Controls::outputType computeAmplitudeType2(const typename Controls::inputType &in){
  Controls controls(in);
  int off = idxOffset(2);
    
  for(int t=0;t<controls.size();t++){
#define A(I) controls.A(I-1,t)
#define C(IDX,I,G1,J,G2) controls.C(IDX-off,I,G1,J,G2,t)

    A(1) = 3.*C(7, 0,VmA,1,VpA) -3.*C(10, 1,VpA,0,VmA);
    A(2) = 3.*C(8, 0,VmA,1,VpA) -3.*C(11, 1,VpA,0,VmA);
    A(3) = 3.*C(7, 0,VmA,1,VpA) +3.*C(7, 0,VmA,0,VmA) -3.*C(10, 1,VpA,0,VmA) -3.*C(10, 0,VmA,0,VmA);
    A(4) = 3.*C(8, 0,VmA,1,VpA) +3.*C(9, 0,VmA,0,VmA) -3.*C(11, 1,VpA,0,VmA) -3.*C(12, 0,VmA,0,VmA);
    A(5) = 3.*C(7, 0,VmA,1,VmA) +3.*C(7, 0,VmA,0,VpA) -3.*C(10, 1,VmA,0,VmA) -3.*C(10, 0,VpA,0,VmA);
    A(6) = 3.*C(8, 0,VmA,1,VmA) +3.*C(9, 0,VmA,0,VpA) -3.*C(11, 1,VmA,0,VmA) -3.*C(12, 0,VpA,0,VmA);
    A(7) = 3.*C(7, 0,VmA,1,VmA) -1.5*C(7, 0,VmA,0,VpA) -3.*C(10, 1,VmA,0,VmA) +1.5*C(10, 0,VpA,0,VmA);
    A(8) = 3.*C(8, 0,VmA,1,VmA) -1.5*C(9, 0,VmA,0,VpA) -3.*C(11, 1,VmA,0,VmA) +1.5*C(12, 0,VpA,0,VmA);
    A(9) = 3.*C(7, 0,VmA,1,VpA) -1.5*C(7, 0,VmA,0,VmA) -3.*C(10, 1,VpA,0,VmA) +1.5*C(10, 0,VmA,0,VmA); 
    A(10) = 3.*C(8, 0,VmA,1,VpA) -1.5*C(9, 0,VmA,0,VmA) -3.*C(11, 1,VpA,0,VmA) +1.5*C(12, 0,VmA,0,VmA);
#undef C
#undef A
    controls.normalize(1./sqrt(6.),t);
  }
  return controls.out;
}

template<typename Controls>
typename Controls::outputType computeAmplitudeType3(const typename Controls::inputType &in){
  Controls controls(in);
  int off = idxOffset(3);
    
  for(int t=0;t<controls.size();t++){
#define A(I) controls.A(I-1,t)
#define C(IDX,I,G1,J,G2) controls.C(IDX-off,I,G1,J,G2,t)

    A(1) = 3.*C(13, 0,VmA,1,VpA) -3.*C(16, 1,VpA,0,VmA);      
    A(2) = 3.*C(14, 0,VmA,1,VpA) -3.*C(17, 1,VpA,0,VmA);      
    A(3) = 3.*C(13, 0,VmA,1,VpA) + 3.*C(13, 0,VmA,0,VmA) -3.*C(16, 1,VpA,0,VmA) -3.*C(16, 0,VmA,0,VmA) + 3.*C(19, 0,VmA,0,VmA) -3.*C(21, 0,VmA,0,VmA);    
    A(4) = 3.*C(14, 0,VmA,1,VpA) +3.*C(15, 0,VmA,0,VmA) -3.*C(17, 1,VpA,0,VmA) -3.*C(18, 0,VmA,0,VmA) +3.*C(20, 0,VmA,0,VmA) -3.*C(22, 0,VmA,0,VmA); 
    A(5) = 3.*C(13, 0,VmA,1,VmA) + 3.*C(13, 0,VmA,0,VpA) -3.*C(16, 1,VmA,0,VmA) -3.*C(16, 0,VpA,0,VmA) + 3.*C(19, 0,VmA,0,VpA) -3.*C(21, 0,VmA,0,VpA);      
    A(6) = 3.*C(14, 0,VmA,1,VmA) +3.*C(15, 0,VmA,0,VpA) -3.*C(17, 1,VmA,0,VmA) -3.*C(18, 0,VpA,0,VmA) +3.*C(20, 0,VmA,0,VpA) -3.*C(22, 0,VmA,0,VpA);
    A(7) = 3.*C(13, 0,VmA,1,VmA) -1.5*C(13, 0,VmA,0,VpA) -3.*C(16, 1,VmA,0,VmA) +1.5*C(16, 0,VpA,0,VmA) -1.5*C(19, 0,VmA,0,VpA) +1.5*C(21, 0,VmA,0,VpA);                                             
    A(8) = 3.*C(14, 0,VmA,1,VmA) -1.5*C(15, 0,VmA,0,VpA) -3.*C(17, 1,VmA,0,VmA) +1.5*C(18, 0,VpA,0,VmA) -1.5*C(20, 0,VmA,0,VpA) +1.5*C(22, 0,VmA,0,VpA);
    A(9) = 3.*C(13, 0,VmA,1,VpA) -1.5*C(13, 0,VmA,0,VmA) -3.*C(16, 1,VpA,0,VmA) +1.5*C(16, 0,VmA,0,VmA) -1.5*C(19, 0,VmA,0,VmA) +1.5*C(21, 0,VmA,0,VmA);
    A(10) = 3.*C(14, 0,VmA,1,VpA) -1.5*C(15, 0,VmA,0,VmA) -3.*C(17, 1,VpA,0,VmA) +1.5*C(18, 0,VmA,0,VmA) -1.5*C(20, 0,VmA,0,VmA) +1.5*C(22, 0,VmA,0,VmA);
#undef C
#undef A
    controls.normalize(1./sqrt(6.),t);
  }
  return controls.out;
}


template<typename Controls>
typename Controls::outputType computeAmplitudeType4(const typename Controls::inputType &in){
  Controls controls(in);
  int off = idxOffset(4);
    
  for(int t=0;t<controls.size();t++){
#define A(I) controls.A(I-1,t)
#define C(IDX,I,G1,J,G2) controls.C(IDX-off,I,G1,J,G2,t)

    A(1) = 3.*C(23, 0,VmA,1,VpA) -3.*C(26, 1,VpA,0,VmA);
    A(2) = 3.*C(24, 0,VmA,1,VpA) -3.*C(27, 1,VpA,0,VmA);
    A(3) = 3.*C(23, 0,VmA,1,VpA) +3.*C(23, 0,VmA,0,VmA) -3.*C(26, 1,VpA,0,VmA) - 3.*C(26, 0,VmA,0,VmA) +3.*C(29, 0,VmA,0,VmA) -3.*C(31, 0,VmA,0,VmA);
    A(4) = 3.*C(24, 0,VmA,1,VpA) +3.*C(25, 0,VmA,0,VmA) -3.*C(27, 1,VpA,0,VmA) -3.*C(28, 0,VmA,0,VmA) + 3.*C(30, 0,VmA,0,VmA) -3.*C(32, 0,VmA,0,VmA);
    A(5) = 3.*C(23, 0,VmA,1,VmA) +3.*C(23, 0,VmA,0,VpA) -3.*C(26, 1,VmA,0,VmA) -3.*C(26, 0,VpA,0,VmA) +3.*C(29, 0,VmA,0,VpA) -3.*C(31, 0,VmA,0,VpA);
    A(6) = 3.*C(24, 0,VmA,1,VmA) +3.*C(25, 0,VmA,0,VpA) -3.*C(27, 1,VmA,0,VmA) -3.*C(28, 0,VpA,0,VmA) +3.*C(30, 0,VmA,0,VpA) -3.*C(32, 0,VmA,0,VpA);
    A(7) = 3.*C(23, 0,VmA,1,VmA) -1.5*C(23, 0,VmA,0,VpA) -3.*C(26, 1,VmA,0,VmA) +1.5*C(26, 0,VpA,0,VmA) -1.5*C(29, 0,VmA,0,VpA) + 1.5*C(31, 0,VmA,0,VpA);
    A(8) = 3.*C(24, 0,VmA,1,VmA) -1.5*C(25, 0,VmA,0,VpA) -3.*C(27, 1,VmA,0,VmA) +1.5*C(28, 0,VpA,0,VmA) -1.5*C(30, 0,VmA,0,VpA) +1.5*C(32,0,VmA,0,VpA);
    A(9) = 3.*C(23, 0,VmA,1,VpA) -1.5*C(23, 0,VmA,0,VmA) -3.*C(26, 1,VpA,0,VmA) +1.5*C(26, 0,VmA,0,VmA) -1.5*C(29, 0,VmA,0,VmA) +1.5*C(31,0,VmA,0,VmA);
    A(10) = 3.*C(24, 0,VmA,1,VpA) -1.5*C(25, 0,VmA,0,VmA) -3.*C(27, 1,VpA,0,VmA) +1.5*C(28,0,VmA,0,VmA) -1.5*C(30, 0,VmA,0,VmA) +1.5*C(32, 0,VmA,0,VmA);
#undef C
#undef A
    controls.normalize(1./sqrt(6.),t);
  }
  return controls.out;
}

template<typename Controls>
typename Controls::outputType computeAmplitudeType(const int i, const typename Controls::inputType &in){
  switch(i){
  case 1: 
    return computeAmplitudeType1<Controls>(in);
  case 2: 
    return computeAmplitudeType2<Controls>(in);
  case 3: 
    return computeAmplitudeType3<Controls>(in);
  case 4: 
    return computeAmplitudeType4<Controls>(in);
  default:
    assert(0);
  }
};

SARLAC_END_NAMESPACE

#endif
