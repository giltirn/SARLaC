#ifndef _FIT_KTOPIPI_KTOSIGMA_GPARITY_SIMFIT_H
#define _FIT_KTOPIPI_KTOSIGMA_GPARITY_SIMFIT_H

#include<config.h>
#include<utils/macros.h>

#include<fit_ktopipi_ktosigma_combined_gparity/simfit_plot.h>

CPSFIT_START_NAMESPACE

void simultaneousFit2state(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &ktopipi_A0_all_j,
			   const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &ktopipi_A0_all_dj,
			   const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &ktosigma_A0_all_j,
			   const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &ktosigma_A0_all_dj,
			   const jackknifeDistributionD &mK, const jackknifeDistributionD &cK, 
			   const jackknifeDistributionD &E0, const jackknifeDistributionD &E1,
			   const NumericSquareMatrix<jackknifeDistributionD> &coeffs, //row = (0=pipi, 1=sigma)  col = (0=gnd state, 1=exc state)
			   const int Lt, const int tmin_k_op, const int tmin_op_snk, bool correlated){
  int nsample = mK.size();

  std::unordered_map<std::string,size_t> param_idx_map;
  int idx = 0;
#define DEFMAP(NM) param_idx_map[#NM] = idx++
  DEFMAP(AK);
  DEFMAP(mK);
  DEFMAP(Apipi0);
  DEFMAP(Apipi1);
  DEFMAP(Asigma0);
  DEFMAP(Asigma1);
  DEFMAP(E0);
  DEFMAP(E1);
  DEFMAP(M0);
  DEFMAP(M1);
#undef DEFMAP

  std::unordered_map<std::string, std::string> inner_param_map_ktopipi;
  inner_param_map_ktopipi["Asnk0"] = "Apipi0";
  inner_param_map_ktopipi["Asnk1"] = "Apipi1";

  std::unordered_map<std::string, std::string> inner_param_map_ktosigma;
  inner_param_map_ktosigma["Asnk0"] = "Asigma0";
  inner_param_map_ktosigma["Asnk1"] = "Asigma1";

  std::vector<std::string> dcp = { "AK", "mK", "E0", "E1", "M0", "M1" };
  for(auto it = dcp.begin(); it != dcp.end(); it++)
    inner_param_map_ktopipi[*it] = inner_param_map_ktosigma[*it] = *it;
    

  std::vector<correlationFunction<SimFitCoordGen, jackknifeDistributionD> > A0_sim_j(10);
  std::vector<correlationFunction<SimFitCoordGen, doubleJackknifeDistributionD> > A0_sim_dj(10);
    
  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > const* jacks[2] = { &ktopipi_A0_all_j, &ktosigma_A0_all_j };
  std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> > const* djacks[2] = { &ktopipi_A0_all_dj, &ktosigma_A0_all_dj };
  std::unordered_map<std::string, std::string> const* pmaps[2] = { &inner_param_map_ktopipi, &inner_param_map_ktosigma };

  for(int q=0;q<10;q++){
    for(int o=0;o<2;o++){
      const auto &jack = (*jacks[o])[q];
      const auto &djack = (*djacks[o])[q];

      for(int d=0;d<jack.size();d++){
	const int t = int(jack.coord(d).t);
	const int tsep_k_snk = jack.coord(d).tsep_k_pi;
	const int tsep_op_snk = tsep_k_snk - t;
	if(t <= tsep_k_snk && t >= tmin_k_op && tsep_op_snk >= tmin_op_snk){
	  SimFitCoordGen c(t, tsep_k_snk, pmaps[o]);
	  A0_sim_j[q].push_back(c, jack.value(d));
	  A0_sim_dj[q].push_back(c, djack.value(d));
	}
      }
    }
  }
  
  typedef taggedValueContainer<double,std::string> Params;
  Params guess(&param_idx_map);
  for(int i=0;i<guess.size();i++) guess(i) = 1.;
  guess("M0") = 0.5;
  guess("M1") = 0.5;

  typedef FitSimGenTwoState FitFunc;

  FitFunc fitfunc(param_idx_map.size());

  typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, correlatedFitPolicy>::type FitPolicies;
    
  std::vector<jackknifeDistribution<Params> > params(10, jackknifeDistribution<Params>(nsample, guess));
  std::vector<jackknifeDistributionD> chisq(10, jackknifeDistributionD(nsample));
  std::vector<jackknifeDistributionD> chisq_per_dof(10, jackknifeDistributionD(nsample));

  std::vector<int> freeze_params = { 0,1,2,3,4,5,6,7 };
  jackknifeDistribution<Params> freeze_vals(nsample, Params(&param_idx_map));
  for(int s=0;s<nsample;s++){
    freeze_vals.sample(s)("AK") = cK.sample(s);
    freeze_vals.sample(s)("mK") = mK.sample(s);
    freeze_vals.sample(s)("Apipi0") = coeffs(0,0).sample(s);
    freeze_vals.sample(s)("Apipi1") = coeffs(0,1).sample(s);
    freeze_vals.sample(s)("Asigma0") = coeffs(1,0).sample(s);
    freeze_vals.sample(s)("Asigma1") = coeffs(1,1).sample(s);
    freeze_vals.sample(s)("E0") = E0.sample(s);
    freeze_vals.sample(s)("E1") = E1.sample(s);
  }

  for(int q=0;q<10;q++){
    std::cout << "Performing " << q+1 << " fit" << std::endl;
    fitter<FitPolicies> fit;
    fit.importFitFunc(fitfunc);
    fit.freeze(freeze_params, freeze_vals);
      
    importCostFunctionParameters<correlatedFitPolicy, FitPolicies> import(fit, A0_sim_dj[q]);
    if(!correlated) import.setUncorrelated();
          
    fit.fit(params[q], chisq[q], chisq_per_dof[q], A0_sim_j[q]);
  }

  for(int q=0;q<10;q++){
    std::cout << "Q" << q+1 << std::endl;
    std::cout << "Params:\n" << params[q] << std::endl;
    std::cout << "Chisq: " << chisq[q] << std::endl;
    std::cout << "Chisq/dof: " << chisq_per_dof[q] << std::endl;
  }
  writeParamsStandard(params, "params.hdf5");
  writeParamsStandard(chisq, "chisq.hdf5");
  writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");

  plotErrorWeightedData2expFlat(ktopipi_A0_all_j, ktosigma_A0_all_j, params, Lt, tmin_k_op, tmin_op_snk, fitfunc);
}




void simultaneousFit2stateDiffM1(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &ktopipi_A0_all_j,
				 const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &ktopipi_A0_all_dj,
				 const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &ktosigma_A0_all_j,
				 const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &ktosigma_A0_all_dj,
				 const jackknifeDistributionD &mK, const jackknifeDistributionD &cK, 
				 const jackknifeDistributionD &E0, const jackknifeDistributionD &E1,
				 const NumericSquareMatrix<jackknifeDistributionD> &coeffs, //row = (0=pipi, 1=sigma)  col = (0=gnd state, 1=exc state)
				 const int Lt, const int tmin_k_op, const int tmin_op_snk, bool correlated){
  int nsample = mK.size();

  std::unordered_map<std::string,size_t> param_idx_map;
  int idx = 0;
#define DEFMAP(NM) param_idx_map[#NM] = idx++
  DEFMAP(AK);
  DEFMAP(mK);
  DEFMAP(Apipi0);
  DEFMAP(Apipi1);
  DEFMAP(Asigma0);
  DEFMAP(Asigma1);
  DEFMAP(E0);
  DEFMAP(E1);
  DEFMAP(M0);
  DEFMAP(M1pipi);
  DEFMAP(M1sigma);
#undef DEFMAP

  std::unordered_map<std::string, std::string> inner_param_map_ktopipi;
  inner_param_map_ktopipi["Asnk0"] = "Apipi0";
  inner_param_map_ktopipi["Asnk1"] = "Apipi1";
  inner_param_map_ktopipi["M1"] = "M1pipi";

  std::unordered_map<std::string, std::string> inner_param_map_ktosigma;
  inner_param_map_ktosigma["Asnk0"] = "Asigma0";
  inner_param_map_ktosigma["Asnk1"] = "Asigma1";
  inner_param_map_ktosigma["M1"] = "M1sigma";

  std::vector<std::string> dcp = { "AK", "mK", "E0", "E1", "M0"};
  for(auto it = dcp.begin(); it != dcp.end(); it++)
    inner_param_map_ktopipi[*it] = inner_param_map_ktosigma[*it] = *it;
    

  std::vector<correlationFunction<SimFitCoordGen, jackknifeDistributionD> > A0_sim_j(10);
  std::vector<correlationFunction<SimFitCoordGen, doubleJackknifeDistributionD> > A0_sim_dj(10);
    
  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > const* jacks[2] = { &ktopipi_A0_all_j, &ktosigma_A0_all_j };
  std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> > const* djacks[2] = { &ktopipi_A0_all_dj, &ktosigma_A0_all_dj };
  std::unordered_map<std::string, std::string> const* pmaps[2] = { &inner_param_map_ktopipi, &inner_param_map_ktosigma };

  for(int q=0;q<10;q++){
    for(int o=0;o<2;o++){
      const auto &jack = (*jacks[o])[q];
      const auto &djack = (*djacks[o])[q];

      for(int d=0;d<jack.size();d++){
	const int t = int(jack.coord(d).t);
	const int tsep_k_snk = jack.coord(d).tsep_k_pi;
	const int tsep_op_snk = tsep_k_snk - t;
	if(t <= tsep_k_snk && t >= tmin_k_op && tsep_op_snk >= tmin_op_snk){
	  SimFitCoordGen c(t, tsep_k_snk, pmaps[o]);
	  A0_sim_j[q].push_back(c, jack.value(d));
	  A0_sim_dj[q].push_back(c, djack.value(d));
	}
      }
    }
  }
  
  typedef taggedValueContainer<double,std::string> Params;
  Params guess(&param_idx_map);
  for(int i=0;i<guess.size();i++) guess(i) = 1.;
  guess("M0") = 0.5;
  guess("M1pipi") = 0.5;
  guess("M1sigma") = 0.5;

  typedef FitSimGenTwoState FitFunc;

  FitFunc fitfunc(param_idx_map.size());

  typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, correlatedFitPolicy>::type FitPolicies;
    
  std::vector<jackknifeDistribution<Params> > params(10, jackknifeDistribution<Params>(nsample, guess));
  std::vector<jackknifeDistributionD> chisq(10, jackknifeDistributionD(nsample));
  std::vector<jackknifeDistributionD> chisq_per_dof(10, jackknifeDistributionD(nsample));

  std::vector<int> freeze_params = { 0,1,2,3,4,5,6,7 };
  jackknifeDistribution<Params> freeze_vals(nsample, Params(&param_idx_map));
  for(int s=0;s<nsample;s++){
    freeze_vals.sample(s)("AK") = cK.sample(s);
    freeze_vals.sample(s)("mK") = mK.sample(s);
    freeze_vals.sample(s)("Apipi0") = coeffs(0,0).sample(s);
    freeze_vals.sample(s)("Apipi1") = coeffs(0,1).sample(s);
    freeze_vals.sample(s)("Asigma0") = coeffs(1,0).sample(s);
    freeze_vals.sample(s)("Asigma1") = coeffs(1,1).sample(s);
    freeze_vals.sample(s)("E0") = E0.sample(s);
    freeze_vals.sample(s)("E1") = E1.sample(s);
  }

  for(int q=0;q<10;q++){
    std::cout << "Performing " << q+1 << " fit" << std::endl;
    fitter<FitPolicies> fit;
    fit.importFitFunc(fitfunc);
    fit.freeze(freeze_params, freeze_vals);
      
    importCostFunctionParameters<correlatedFitPolicy, FitPolicies> import(fit, A0_sim_dj[q]);
    if(!correlated) import.setUncorrelated();
          
    fit.fit(params[q], chisq[q], chisq_per_dof[q], A0_sim_j[q]);
  }

  for(int q=0;q<10;q++){
    std::cout << "Q" << q+1 << std::endl;
    std::cout << "Params:\n" << params[q] << std::endl;
    std::cout << "Chisq: " << chisq[q] << std::endl;
    std::cout << "Chisq/dof: " << chisq_per_dof[q] << std::endl;
  }
  writeParamsStandard(params, "params.hdf5");
  writeParamsStandard(chisq, "chisq.hdf5");
  writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");

  plotErrorWeightedData2expFlatDiffM1(ktopipi_A0_all_j, ktosigma_A0_all_j, params, Lt, tmin_k_op, tmin_op_snk, fitfunc);
}



CPSFIT_END_NAMESPACE

#endif
