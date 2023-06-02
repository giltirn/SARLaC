#pragma once
#include "pipi_bubble_mom_data.h"

CPSFIT_START_NAMESPACE

//Averages bubble over pion momenta to produce a rotationally invariant state

template<typename RealOrComplex>
struct _getReOrZ{
  static inline double get(const std::complex<double> &val){ return std::real(val); }
};
template<typename T>
struct _getReOrZ<std::complex<T> >{
  static inline std::complex<double> get(const std::complex<double> &val){ return val; }
};
template<typename BubbleDataType>
struct projectSourcePiPiBubble_getCoeff{
  typedef typename BubbleDataType::DistributionType::DataType DataType;
  static inline auto getCoeff(const std::complex<double> &val)->decltype(_getReOrZ<DataType>::get(val)){ return _getReOrZ<DataType>::get(val); }
};

//Read and project pipi bubble
template<typename ContainerType, typename ReadPolicy>
void getProjectedSourcePiPiBubble(ContainerType &out, const int Lt, const int tsep_pipi, const PiPiProjectorBase &proj_pipi, const ReadPolicy &rp){
  out.setup(Source, Lt, tsep_pipi);
  
  ContainerType temp(Source, Lt, tsep_pipi);
  for(int p=0;p<proj_pipi.nMomenta();p++){
    readPiPiBubbleSingleMom(temp, Lt, proj_pipi.momentum(p), rp);
    for(int t=0;t<Lt;t++){
      typename ContainerType::DistributionType tmp =  projectSourcePiPiBubble_getCoeff<ContainerType>::getCoeff(proj_pipi.coefficient(p)) * temp(t);
      out(t) = p==0 ? tmp : out(t) + tmp;
    }
  }
}
//Call the above with the default filename policy
template<typename ContainerType>
void getProjectedSourcePiPiBubble(ContainerType &out, const std::string &data_dir, const int traj_start, const int traj_inc, const int traj_lessthan, 
				  const int Lt, const int tsep_pipi,  const PiPiProjectorBase &proj_pipi){
  readBubbleStationaryPolicy fp(false,Source);
  PiPiBubbleBasicReadPolicy<readBubbleStationaryPolicy> rp(fp, data_dir, traj_start, traj_inc, traj_lessthan, tsep_pipi);
  getProjectedSourcePiPiBubble(out, Lt, tsep_pipi, proj_pipi, rp);
}

//Get the projected pipi bubble using bubble data stored in an existing container
template<typename AllMomentaContainerType>
typename AllMomentaContainerType::ContainerType projectSourcePiPiBubble(const AllMomentaContainerType &pipi_self_data, const PiPiProjectorBase &proj_pipi){
  int Lt = pipi_self_data.getLt();
  int tsep_pipi = pipi_self_data.getTsepPiPi();
  typedef typename AllMomentaContainerType::ContainerType OutType;
  OutType out(Source,Lt,tsep_pipi);

  for(int t=0;t<Lt;t++)
    for(int p=0;p<proj_pipi.nMomenta();p++){
      typename AllMomentaContainerType::DistributionType tmp = projectSourcePiPiBubble_getCoeff<OutType>::getCoeff(proj_pipi.coefficient(p)) * pipi_self_data(Source,proj_pipi.momentum(p))(t);
      out(t) = p == 0 ? tmp : out(t) + tmp;
    }
  
  return out;
}

CPSFIT_END_NAMESPACE
