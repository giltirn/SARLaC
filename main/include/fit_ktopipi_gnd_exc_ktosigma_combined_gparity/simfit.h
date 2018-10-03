#ifndef _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_H
#define _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_H

#include<config.h>
#include<utils/macros.h>

#include "args.h"

CPSFIT_START_NAMESPACE

struct simultaneousFitBase{
  typedef correlationFunction<amplitudeDataCoord, jackknifeDistributionD> CorrFuncJack;
  typedef correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> CorrFuncDJack;

  typedef std::vector<CorrFuncJack> CorrFuncJackAllQ;
  typedef std::vector<CorrFuncDJack> CorrFuncDJackAllQ;

  typedef correlationFunction<SimFitCoordGen, jackknifeDistributionD> SimFitCorrFuncJack;
  typedef correlationFunction<SimFitCoordGen, doubleJackknifeDistributionD> SimFitCorrFuncDJack;

  typedef std::unordered_map<std::string, std::string> InnerParamMap;

  #define COPY_BASE_TYPEDEFS \
    typedef simultaneousFitBase::CorrFuncJack CorrFuncJack; \
    typedef simultaneousFitBase::CorrFuncDJack CorrFuncDJack; \
    typedef simultaneousFitBase::CorrFuncJackAllQ CorrFuncJackAllQ; \
    typedef simultaneousFitBase::CorrFuncDJackAllQ CorrFuncDJackAllQ; \
    typedef simultaneousFitBase::SimFitCorrFuncJack SimFitCorrFuncJack; \
    typedef simultaneousFitBase::SimFitCorrFuncDJack SimFitCorrFuncDJack; \
    typedef simultaneousFitBase::InnerParamMap InnerParamMap

  virtual void load2ptFitParams(const InputParamArgs &iargs, const int nsample) = 0;
  virtual void fit(const CorrFuncJackAllQ &ktopipi_A0_all_j,
		   const CorrFuncDJackAllQ &ktopipi_A0_all_dj,
		   const CorrFuncJackAllQ &ktopipi_exc_A0_all_j,
		   const CorrFuncDJackAllQ &ktopipi_exc_A0_all_dj,
		   const CorrFuncJackAllQ &ktosigma_A0_all_j,
		   const CorrFuncDJackAllQ &ktosigma_A0_all_dj,
		   const int Lt, const int tmin_k_op, const int tmin_op_snk, bool correlated) = 0;


  void generateSimData(std::vector<SimFitCorrFuncJack> &A0_sim_j,
		       std::vector<SimFitCorrFuncDJack> &A0_sim_dj,
		       const std::vector<CorrFuncJackAllQ const* > &jacks,
		       const std::vector<CorrFuncDJackAllQ const* > &djacks,
		       const std::vector<InnerParamMap const *> &pmaps,
		       const int tmin_k_op, const int tmin_op_snk){
    int nops = jacks.size();
    assert(djacks.size() == nops && pmaps.size() == nops);

    for(int q=0;q<10;q++){
      for(int o=0;o<nops;o++){
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
  }

  virtual ~simultaneousFitBase(){}
};


class simultaneousFit2state: public simultaneousFitBase{
public:
  COPY_BASE_TYPEDEFS;    
    
private:
  bool loaded_frzparams;
  jackknifeDistributionD mK;
  jackknifeDistributionD cK;
  
  jackknifeDistributionD E0, E1;
  //row = (0=pipi, 1=pipi_exc 2=sigma)  col = (0=gnd state, 1=exc state)
  NumericTensor<jackknifeDistributionD,2> coeffs;  
public:

  void load2ptFitParams(const InputParamArgs &iargs, const int nsample){
    {
      std::vector<jackknifeDistributionD> p;
      readParamsStandard(p,  iargs.kaon2pt_fit_result);
      mK = p[iargs.idx_mK];
      cK = sqrt( p[iargs.idx_cK] );
    }

    //row = (0=pipi, 1=pipi_exc 2=sigma)  col = (0=gnd state, 1=exc state)
    {
      double scale = sqrt(iargs.pipi_sigma_sim_fit_Ascale);
      std::vector<jackknifeDistributionD> p;
      readParamsStandard(p, iargs.pipi_sigma_sim_fit_result);
      for(int i=0;i<p.size();i++) assert(p[i].size() == nsample);
      coeffs({0,0}) = p[iargs.idx_coeff_pipi_state0] * scale;
      coeffs({0,1}) = p[iargs.idx_coeff_pipi_state1] * scale;
      coeffs({1,0}) = p[iargs.idx_coeff_pipi_exc_state0] * scale;
      coeffs({1,1}) = p[iargs.idx_coeff_pipi_exc_state1] * scale;
      coeffs({2,0}) = p[iargs.idx_coeff_sigma_state0] * scale;
      coeffs({2,1}) = p[iargs.idx_coeff_sigma_state1] * scale;
      E0 = p[iargs.idx_E0];
      E1 = p[iargs.idx_E1];
    }

    std::cout << "cK = " << cK << std::endl;
    std::cout << "mK = " << mK << std::endl;
    std::cout << "E0 = " << E0 << std::endl;
    std::cout << "E1 = " << E1 << std::endl;
    loaded_frzparams = true;
  }
  
  simultaneousFit2state(): loaded_frzparams(false), coeffs({3,2}){}

  void fit(const CorrFuncJackAllQ &ktopipi_A0_all_j,
	   const CorrFuncDJackAllQ &ktopipi_A0_all_dj,
	   const CorrFuncJackAllQ &ktopipi_exc_A0_all_j,
	   const CorrFuncDJackAllQ &ktopipi_exc_A0_all_dj,
	   const CorrFuncJackAllQ &ktosigma_A0_all_j,
	   const CorrFuncDJackAllQ &ktosigma_A0_all_dj,
	   const int Lt, const int tmin_k_op, const int tmin_op_snk, bool correlated){
    assert(loaded_frzparams);

    int nsample = mK.size();

    std::unordered_map<std::string,size_t> param_idx_map;
    int idx = 0;
#define DEFMAP(NM) param_idx_map[#NM] = idx++
    DEFMAP(AK);
    DEFMAP(mK);
    DEFMAP(Apipi0);
    DEFMAP(Apipi1);
    DEFMAP(Apipi_exc_0);
    DEFMAP(Apipi_exc_1);
    DEFMAP(Asigma0);
    DEFMAP(Asigma1);
    DEFMAP(E0);
    DEFMAP(E1);
    DEFMAP(M0);
    DEFMAP(M1);
#undef DEFMAP

    InnerParamMap inner_param_map_ktopipi;
    inner_param_map_ktopipi["Asnk0"] = "Apipi0";
    inner_param_map_ktopipi["Asnk1"] = "Apipi1";

    InnerParamMap inner_param_map_ktopipi_exc;
    inner_param_map_ktopipi_exc["Asnk0"] = "Apipi_exc_0";
    inner_param_map_ktopipi_exc["Asnk1"] = "Apipi_exc_1";

    InnerParamMap inner_param_map_ktosigma;
    inner_param_map_ktosigma["Asnk0"] = "Asigma0";
    inner_param_map_ktosigma["Asnk1"] = "Asigma1";

    std::vector<std::string> dcp = { "AK", "mK", "E0", "E1", "M0", "M1" };
    for(auto it = dcp.begin(); it != dcp.end(); it++)
      inner_param_map_ktopipi[*it] = inner_param_map_ktopipi_exc[*it] = inner_param_map_ktosigma[*it] = *it;
    
    std::vector<SimFitCorrFuncJack> A0_sim_j(10);
    std::vector<SimFitCorrFuncDJack> A0_sim_dj(10);
    
    generateSimData(A0_sim_j,A0_sim_dj,
		    {&ktopipi_A0_all_j, &ktopipi_exc_A0_all_j, &ktosigma_A0_all_j},
		    {&ktopipi_A0_all_dj, &ktopipi_exc_A0_all_dj, &ktosigma_A0_all_dj},
		    {&inner_param_map_ktopipi, &inner_param_map_ktopipi_exc, &inner_param_map_ktosigma},
		    tmin_k_op, tmin_op_snk);
  
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

    std::vector<int> freeze_params = { 0,1,2,3,4,5,6,7,8,9 };
    jackknifeDistribution<Params> freeze_vals(nsample, Params(&param_idx_map));
    int opidx_ktopipi = 0;
    int opidx_ktopipi_exc = 1;
    int opidx_ktosigma = 2;

    for(int s=0;s<nsample;s++){
      freeze_vals.sample(s)("AK") = cK.sample(s);
      freeze_vals.sample(s)("mK") = mK.sample(s);
      freeze_vals.sample(s)("Apipi0") = coeffs({opidx_ktopipi,0}).sample(s);
      freeze_vals.sample(s)("Apipi1") = coeffs({opidx_ktopipi,1}).sample(s);
      freeze_vals.sample(s)("Apipi_exc_0") = coeffs({opidx_ktopipi_exc,0}).sample(s);
      freeze_vals.sample(s)("Apipi_exc_1") = coeffs({opidx_ktopipi_exc,1}).sample(s);
      freeze_vals.sample(s)("Asigma0") = coeffs({opidx_ktosigma,0}).sample(s);
      freeze_vals.sample(s)("Asigma1") = coeffs({opidx_ktosigma,1}).sample(s);
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

    plotErrorWeightedData2expFlat(ktopipi_A0_all_j, ktopipi_exc_A0_all_j, ktosigma_A0_all_j, params, Lt, tmin_k_op, tmin_op_snk, fitfunc);
  }
};




class simultaneousFit3state: public simultaneousFitBase{
public:
  COPY_BASE_TYPEDEFS;    
    
private:
  bool loaded_frzparams;
  jackknifeDistributionD mK;
  jackknifeDistributionD cK;
  
  jackknifeDistributionD E0, E1, E2;
  //row = (0=pipi, 1=pipi_exc 2=sigma)  col = (0=gnd state, 1=exc state 1, 2=exc state 2)
  NumericTensor<jackknifeDistributionD,2> coeffs;  
public:

  void load2ptFitParams(const InputParamArgs &iargs, const int nsample){
    {
      std::vector<jackknifeDistributionD> p;
      readParamsStandard(p,  iargs.kaon2pt_fit_result);
      mK = p[iargs.idx_mK];
      cK = sqrt( p[iargs.idx_cK] );
    }

    {
      double scale = sqrt(iargs.pipi_sigma_sim_fit_Ascale);
      std::vector<jackknifeDistributionD> p;
      readParamsStandard(p, iargs.pipi_sigma_sim_fit_result);
      for(int i=0;i<p.size();i++) assert(p[i].size() == nsample);
      coeffs({0,0}) = p[iargs.idx_coeff_pipi_state0] * scale;
      coeffs({0,1}) = p[iargs.idx_coeff_pipi_state1] * scale;
      coeffs({0,2}) = p[iargs.idx_coeff_pipi_state2] * scale;
      coeffs({1,0}) = p[iargs.idx_coeff_pipi_exc_state0] * scale;
      coeffs({1,1}) = p[iargs.idx_coeff_pipi_exc_state1] * scale;
      coeffs({1,2}) = p[iargs.idx_coeff_pipi_exc_state2] * scale;
      coeffs({2,0}) = p[iargs.idx_coeff_sigma_state0] * scale;
      coeffs({2,1}) = p[iargs.idx_coeff_sigma_state1] * scale;
      coeffs({2,2}) = p[iargs.idx_coeff_sigma_state2] * scale;
      E0 = p[iargs.idx_E0];
      E1 = p[iargs.idx_E1];
      E2 = p[iargs.idx_E2];
    }

    std::cout << "cK = " << cK << std::endl;
    std::cout << "mK = " << mK << std::endl;
    std::cout << "E0 = " << E0 << std::endl;
    std::cout << "E1 = " << E1 << std::endl;
    std::cout << "E2 = " << E2 << std::endl;

    std::cout << "Coeffs:\n" << coeffs << std::endl;

    loaded_frzparams = true;
  }
  
  simultaneousFit3state(): loaded_frzparams(false), coeffs({3,3}){}

  void fit(const CorrFuncJackAllQ &ktopipi_A0_all_j,
	   const CorrFuncDJackAllQ &ktopipi_A0_all_dj,
	   const CorrFuncJackAllQ &ktopipi_exc_A0_all_j,
	   const CorrFuncDJackAllQ &ktopipi_exc_A0_all_dj,
	   const CorrFuncJackAllQ &ktosigma_A0_all_j,
	   const CorrFuncDJackAllQ &ktosigma_A0_all_dj,
	   const int Lt, const int tmin_k_op, const int tmin_op_snk, bool correlated){
    assert(loaded_frzparams);

    int nsample = mK.size();

    std::unordered_map<std::string,size_t> param_idx_map;
    int idx = 0;
#define DEFMAP(NM) param_idx_map[#NM] = idx++
    DEFMAP(AK);
    DEFMAP(mK);
    DEFMAP(Apipi0);
    DEFMAP(Apipi1);
    DEFMAP(Apipi2);
    DEFMAP(Apipi_exc_0);
    DEFMAP(Apipi_exc_1);
    DEFMAP(Apipi_exc_2);
    DEFMAP(Asigma0);
    DEFMAP(Asigma1);
    DEFMAP(Asigma2);
    DEFMAP(E0);
    DEFMAP(E1);
    DEFMAP(E2);
    DEFMAP(M0);
    DEFMAP(M1);
    DEFMAP(M2);
#undef DEFMAP

    InnerParamMap inner_param_map_ktopipi;
    inner_param_map_ktopipi["Asnk0"] = "Apipi0";
    inner_param_map_ktopipi["Asnk1"] = "Apipi1";
    inner_param_map_ktopipi["Asnk2"] = "Apipi2";

    InnerParamMap inner_param_map_ktopipi_exc;
    inner_param_map_ktopipi_exc["Asnk0"] = "Apipi_exc_0";
    inner_param_map_ktopipi_exc["Asnk1"] = "Apipi_exc_1";
    inner_param_map_ktopipi_exc["Asnk2"] = "Apipi_exc_2";

    InnerParamMap inner_param_map_ktosigma;
    inner_param_map_ktosigma["Asnk0"] = "Asigma0";
    inner_param_map_ktosigma["Asnk1"] = "Asigma1";
    inner_param_map_ktosigma["Asnk2"] = "Asigma2";

    std::vector<std::string> dcp = { "AK", "mK", "E0", "E1", "E2", "M0", "M1", "M2" };
    for(auto it = dcp.begin(); it != dcp.end(); it++)
      inner_param_map_ktopipi[*it] = inner_param_map_ktopipi_exc[*it] = inner_param_map_ktosigma[*it] = *it;
    
    std::vector<SimFitCorrFuncJack> A0_sim_j(10);
    std::vector<SimFitCorrFuncDJack> A0_sim_dj(10);
    
    generateSimData(A0_sim_j,A0_sim_dj,
		    {&ktopipi_A0_all_j, &ktopipi_exc_A0_all_j, &ktosigma_A0_all_j},
		    {&ktopipi_A0_all_dj, &ktopipi_exc_A0_all_dj, &ktosigma_A0_all_dj},
		    {&inner_param_map_ktopipi, &inner_param_map_ktopipi_exc, &inner_param_map_ktosigma},
		    tmin_k_op, tmin_op_snk);
  
    typedef taggedValueContainer<double,std::string> Params;
    Params guess(&param_idx_map);
    for(int i=0;i<guess.size();i++) guess(i) = 1.;
    guess("M0") = 0.5;
    guess("M1") = 0.5;
    guess("M2") = 0.5;

    typedef FitSimGenThreeState FitFunc;

    FitFunc fitfunc(param_idx_map.size());

    typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, correlatedFitPolicy>::type FitPolicies;
    
    std::vector<jackknifeDistribution<Params> > params(10, jackknifeDistribution<Params>(nsample, guess));
    std::vector<jackknifeDistributionD> chisq(10, jackknifeDistributionD(nsample));
    std::vector<jackknifeDistributionD> chisq_per_dof(10, jackknifeDistributionD(nsample));

    std::vector<int> freeze_params = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13 };
    jackknifeDistribution<Params> freeze_vals(nsample, Params(&param_idx_map));
    int opidx_ktopipi = 0;
    int opidx_ktopipi_exc = 1;
    int opidx_ktosigma = 2;

    for(int s=0;s<nsample;s++){
      freeze_vals.sample(s)("AK") = cK.sample(s);
      freeze_vals.sample(s)("mK") = mK.sample(s);
      freeze_vals.sample(s)("Apipi0") = coeffs({opidx_ktopipi,0}).sample(s);
      freeze_vals.sample(s)("Apipi1") = coeffs({opidx_ktopipi,1}).sample(s);
      freeze_vals.sample(s)("Apipi2") = coeffs({opidx_ktopipi,2}).sample(s);
      freeze_vals.sample(s)("Apipi_exc_0") = coeffs({opidx_ktopipi_exc,0}).sample(s);
      freeze_vals.sample(s)("Apipi_exc_1") = coeffs({opidx_ktopipi_exc,1}).sample(s);
      freeze_vals.sample(s)("Apipi_exc_2") = coeffs({opidx_ktopipi_exc,2}).sample(s);
      freeze_vals.sample(s)("Asigma0") = coeffs({opidx_ktosigma,0}).sample(s);
      freeze_vals.sample(s)("Asigma1") = coeffs({opidx_ktosigma,1}).sample(s);
      freeze_vals.sample(s)("Asigma2") = coeffs({opidx_ktosigma,2}).sample(s);
      freeze_vals.sample(s)("E0") = E0.sample(s);
      freeze_vals.sample(s)("E1") = E1.sample(s);
      freeze_vals.sample(s)("E2") = E2.sample(s);
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

    plotErrorWeightedData3expFlat(ktopipi_A0_all_j, ktopipi_exc_A0_all_j, ktosigma_A0_all_j, params, Lt, tmin_k_op, tmin_op_snk, fitfunc);
  }
};




#undef COPY_BASE_TYPEDEFS

CPSFIT_END_NAMESPACE

#endif
