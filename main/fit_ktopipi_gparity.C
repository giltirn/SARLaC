#include<array>
#include<vector>
#include<complex>
#include<iostream>
#include<fstream>

#include <boost/timer/timer.hpp>

#include<utils.h>
#include<distribution.h>
#include<common_defs.h>
#include<numeric_tensors.h>
#include<correlationfunction.h>

#include <fit_pipi_gparity/data_containers.h>
#include <fit_pipi_gparity/mom_data_containers.h>
#include <fit_pipi_gparity/read_data.h>

#include <fit_ktopipi_gparity/data_containers.h>
#include <fit_ktopipi_gparity/read_data.h>
#include <fit_ktopipi_gparity/compute_amplitude.h>
#include <fit_ktopipi_gparity/main.h>

int main(void){
  int Lt = 64;
  int traj_start = 636;
  int traj_inc = 4;
  int traj_lessthan = 648;
  int tsep_k_pi = 10;
  int tsep_pipi = 4;
  std::string data_dir = ".";
  
  type1234Data type1 = readType(1, traj_start, traj_inc, traj_lessthan, tsep_k_pi, tsep_pipi, Lt, data_dir);
  type1234Data type2 = readType(2, traj_start, traj_inc, traj_lessthan, tsep_k_pi, tsep_pipi, Lt, data_dir);
  type1234Data type3 = readType(3, traj_start, traj_inc, traj_lessthan, tsep_k_pi, tsep_pipi, Lt, data_dir);
  type1234Data type4 = readType(4, traj_start, traj_inc, traj_lessthan, tsep_k_pi, tsep_pipi, Lt, data_dir);

  std::vector<int> type1_nonzerotK = type1.getNonZeroKaonTimeslices();
  std::vector<int> type2_nonzerotK = type2.getNonZeroKaonTimeslices();
  std::vector<int> type3_nonzerotK = type3.getNonZeroKaonTimeslices();
  std::vector<int> type4_nonzerotK = type4.getNonZeroKaonTimeslices();

  NumericTensor<rawDataDistributionD,1> bubble = readA2projectedBubble(traj_start,traj_inc,traj_lessthan,tsep_pipi,Lt,data_dir);
  NumericTensor<doubleJackknifeDistributionD,1> bubble_dj = bubble.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());

  //Get the type 4 contribution without the bubble and double-jackknife resample
  NumericTensor<rawDataDistributionD,3> A0_type4_alltK = computeAmplitudeType4<computeAmplitudeAlltKtensorControls>(type4); //[Qidx][tK][t]
  NumericTensor<doubleJackknifeDistributionD,3> A0_type4_alltK_dj = A0_type4_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());

  //Get and double-jackknife resample the mix4 (no bubble) diagram
  NumericTensor<rawDataDistributionD,2> mix4_alltK({Lt,Lt}); //[tK][t]
  for(int tK=0;tK<Lt;tK++) for(int t=0;t<Lt;t++) mix4_alltK({tK,t}) = type4(tK,t).mix();
  NumericTensor<doubleJackknifeDistributionD,2> mix4_alltK_dj = mix4_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());
  
  //tK-average the type4 and mix4 (no bubble) double-jackknife data then double-jackknife resample  
  NumericTensor<doubleJackknifeDistributionD,2> A0_type4_srcavg_dj = A0_type4_alltK_dj.reduce(1, averageDimensionFunctor<doubleJackknifeDistributionD,3>(1, &type4_nonzerotK)); //[Qidx][t]
  NumericTensor<doubleJackknifeDistributionD,1> mix4_srcavg_dj = mix4_alltK_dj.reduce(0, averageDimensionFunctor<doubleJackknifeDistributionD,2>(0, &type4_nonzerotK)); //[t]
   
  //Compute alpha from the above
  NumericTensor<doubleJackknifeDistributionD,2> alpha({10,Lt}, [&](int const *c){ return A0_type4_srcavg_dj({c[0],c[1]})/mix4_srcavg_dj({c[1]}); }); //[Qidx][t]

  //Compute vacuum subtractions
  NumericTensor<doubleJackknifeDistributionD,3> A0_type4_alltK_vacsub = A0_type4_alltK_dj.transform(multiplyBubbleFunctor<doubleJackknifeDistributionD>(bubble_dj,tsep_k_pi,tsep_pipi,Lt,1));//[Qidx][tK][t]
  NumericTensor<doubleJackknifeDistributionD,2> mix4_alltK_vacsub = mix4_alltK_dj.transform(multiplyBubbleFunctor<doubleJackknifeDistributionD>(bubble_dj,tsep_k_pi,tsep_pipi,Lt,0)); //[tK][t]

  //Source average the vacuum subtractions
  NumericTensor<doubleJackknifeDistributionD,2> A0_type4_srcavg_vacsub = A0_type4_alltK_vacsub.reduce(1, averageDimensionFunctor<doubleJackknifeDistributionD,3>(1, &type4_nonzerotK)); //[Qidx][t]
  NumericTensor<doubleJackknifeDistributionD,1> mix4_srcavg_vacsub = mix4_alltK_vacsub.reduce(0, averageDimensionFunctor<doubleJackknifeDistributionD,2>(0, &type4_nonzerotK)); //[t]
  
  //Apply the bubble diagram to the type4 and mix4 data
  A0_type4_alltK = A0_type4_alltK.transform(multiplyBubbleFunctor<rawDataDistributionD>(bubble,tsep_k_pi,tsep_pipi,Lt,1));
  mix4_alltK = mix4_alltK.transform(multiplyBubbleFunctor<rawDataDistributionD>(bubble,tsep_k_pi,tsep_pipi,Lt,0));

  //Recompute double-jackknife resamplings of full type4 and mix4
  A0_type4_alltK_dj = A0_type4_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());
  mix4_alltK_dj = mix4_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());

  //Recompute double-jackknife resamplings of full, source-averaged type4 and mix4
  A0_type4_srcavg_dj = A0_type4_alltK_dj.reduce(1, averageDimensionFunctor<doubleJackknifeDistributionD,3>(1, &type4_nonzerotK)); //[Qidx][t]
  mix4_srcavg_dj = mix4_alltK_dj.reduce(0, averageDimensionFunctor<doubleJackknifeDistributionD,2>(0, &type4_nonzerotK)); //[t]

  //Get, double-jackknife resample and source-average the type1, type2, type3 and mix3 contributions
  NumericTensor<rawDataDistributionD,3> A0_type1_alltK = computeAmplitudeType1<computeAmplitudeAlltKtensorControls>(type1); //[Qidx][tK][t]
  NumericTensor<doubleJackknifeDistributionD,3> A0_type1_alltK_dj = A0_type1_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());
  NumericTensor<doubleJackknifeDistributionD,2> A0_type1_srcavg_dj = A0_type1_alltK_dj.reduce(1, averageDimensionFunctor<doubleJackknifeDistributionD,3>(1, &type1_nonzerotK)); //[Qidx][t]
  
  NumericTensor<rawDataDistributionD,3> A0_type2_alltK = computeAmplitudeType2<computeAmplitudeAlltKtensorControls>(type2); //[Qidx][tK][t]
  NumericTensor<doubleJackknifeDistributionD,3> A0_type2_alltK_dj = A0_type2_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());
  NumericTensor<doubleJackknifeDistributionD,2> A0_type2_srcavg_dj = A0_type2_alltK_dj.reduce(1, averageDimensionFunctor<doubleJackknifeDistributionD,3>(1, &type2_nonzerotK)); //[Qidx][t]

  NumericTensor<rawDataDistributionD,3> A0_type3_alltK = computeAmplitudeType3<computeAmplitudeAlltKtensorControls>(type3); //[Qidx][tK][t]
  NumericTensor<doubleJackknifeDistributionD,3> A0_type3_alltK_dj = A0_type3_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());
  NumericTensor<doubleJackknifeDistributionD,2> A0_type3_srcavg_dj = A0_type3_alltK_dj.reduce(1, averageDimensionFunctor<doubleJackknifeDistributionD,3>(1, &type3_nonzerotK)); //[Qidx][t]

  NumericTensor<rawDataDistributionD,2> mix3_alltK({Lt,Lt}); //[tK][t]
  for(int tK=0;tK<Lt;tK++) for(int t=0;t<Lt;t++) mix3_alltK({tK,t}) = type3(tK,t).mix();
  NumericTensor<doubleJackknifeDistributionD,2> mix3_alltK_dj = mix3_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());
  NumericTensor<doubleJackknifeDistributionD,1> mix3_srcavg_dj = mix3_alltK_dj.reduce(0, averageDimensionFunctor<doubleJackknifeDistributionD,2>(0, &type3_nonzerotK)); //[t]
  
  //Subtract the pseudoscalar operators and mix4 vacuum term
  A0_type3_srcavg_dj = A0_type3_srcavg_dj.transform([&](int const* coord, const doubleJackknifeDistributionD &from){ return doubleJackknifeDistributionD(from - alpha(coord)*mix3_srcavg_dj(coord+1)); }); 
  A0_type4_srcavg_dj = A0_type4_srcavg_dj.transform(
						    [&](int const* coord, const doubleJackknifeDistributionD &from){
						      return doubleJackknifeDistributionD(from - alpha(coord)*( mix4_srcavg_dj(coord+1) + mix4_srcavg_vacsub(coord+1) ) );
						    }
						    ); 
  
  //Perform the type 4 vacuum subtraction
  A0_type4_srcavg_dj = A0_type4_srcavg_dj - A0_type4_srcavg_vacsub;

  //Get the full double-jackknife amplitude
  NumericTensor<doubleJackknifeDistributionD,2> A0_full_srcavg_dj = A0_type1_srcavg_dj + A0_type2_srcavg_dj + A0_type3_srcavg_dj + A0_type4_srcavg_dj;
  
  //Get the single-elimination jackknife distributions
  NumericTensor<jackknifeDistributionD,2> A0_type1_srcavg_j = A0_type1_srcavg_dj.transform([](int const* coord, const doubleJackknifeDistributionD &from){ return from.toJackknife(); });
  NumericTensor<jackknifeDistributionD,2> A0_type2_srcavg_j = A0_type2_srcavg_dj.transform([](int const* coord, const doubleJackknifeDistributionD &from){ return from.toJackknife(); });
  NumericTensor<jackknifeDistributionD,2> A0_type3_srcavg_j = A0_type3_srcavg_dj.transform([](int const* coord, const doubleJackknifeDistributionD &from){ return from.toJackknife(); });
  NumericTensor<jackknifeDistributionD,2> A0_type4_srcavg_j = A0_type4_srcavg_dj.transform([](int const* coord, const doubleJackknifeDistributionD &from){ return from.toJackknife(); });

  NumericTensor<jackknifeDistributionD,2> A0_full_srcavg_j = A0_type1_srcavg_j + A0_type2_srcavg_j + A0_type3_srcavg_j + A0_type4_srcavg_j;

  NumericTensor<jackknifeDistributionD,2> sigma = A0_full_srcavg_dj.transform([](int const *c, const doubleJackknifeDistributionD &d){ return jackknifeDistributionD(sqrt(doubleJackknifeDistributionD::covariance(d,d))); });

  

  



    
  
  return 0;
}


