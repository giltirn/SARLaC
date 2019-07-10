#ifndef _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_H
#define _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_H

#include<config.h>
#include<utils/macros.h>

#include "args.h"
#include<pipi_common/analyze_chisq.h>

CPSFIT_START_NAMESPACE

struct printCoord{
  typedef std::unordered_map<std::string, std::string> InnerParamMap;
  const SimFitCoordGen &c;
  const std::map< InnerParamMap const*, std::string> &pmap_descr;
  printCoord(const SimFitCoordGen &c, const std::map< InnerParamMap const*, std::string> &pmap_descr): c(c), pmap_descr(pmap_descr){}
};
std::ostream & operator<<(std::ostream &os, const printCoord &p){ 
  auto it = p.pmap_descr.find(p.c.param_map); assert(it != p.pmap_descr.end());
  os << "(" << it->second << ", t=" << p.c.t << " tsep_K_snk=" << p.c.tsep_k_snk << ")";
  return os;
}


template<typename FitFunc>
void analyzeChisqFF(const correlationFunction<SimFitCoordGen,  jackknifeDistributionD> &corr_comb_j,
		    const jackknifeDistribution<typename FitFunc::ParameterType> &params, const FitFunc &fitfunc,
		    const std::map<std::unordered_map<std::string, std::string> const*, std::string> &pmap_descr){
  struct PP{
    typedef std::map<std::unordered_map<std::string, std::string> const*, std::string> const* PtrType;
    inline static PtrType & descr(){ static PtrType p; return p; }

    inline static void print(std::ostream &os, const SimFitCoordGen &c){ os << printCoord(c, *descr()); }
    inline static std::string typeInfo(const SimFitCoordGen &c){ return descr()->find(c.param_map)->second; }   
  };
  PP::descr() = &pmap_descr;
  
  AnalyzeChisq<FitFunc,PP> chisq_analyze(corr_comb_j, fitfunc, params);
  chisq_analyze.printChisqContribs(Correlation);
  chisq_analyze.examineProjectedDeviationContribsEvalNorm(Correlation);
  chisq_analyze.printChisqContribs(Covariance);
  chisq_analyze.examineProjectedDeviationContribsEvalNorm(Covariance);
}

struct simultaneousFitBase{
  typedef correlationFunction<amplitudeDataCoord, jackknifeDistributionD> CorrFuncJack;
  typedef correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> CorrFuncDJack;
  
  typedef std::vector<CorrFuncJack> CorrFuncJackAllQ;
  typedef std::vector<CorrFuncDJack> CorrFuncDJackAllQ;
  
  typedef correlationFunction<SimFitCoordGen, jackknifeDistributionD> SimFitCorrFuncJack;
  typedef correlationFunction<SimFitCoordGen, doubleJackknifeDistributionD> SimFitCorrFuncDJack;
  
  typedef taggedValueContainer<double,std::string> Params;
  typedef std::unordered_map<std::string, std::string> InnerParamMap;
  
  //All operators have the same fit form, however which shared parameters the inner parameters map to is operator dependent
  //To perform a fit we need the following:
  //   A map of a shared parameter name to an index  map(string -> size_t)  aka  paramIdxMap
  //   A map between an the operator and the subset of shared parameters it uses map(op ->  map(string -> string) ) aka  operatorSubsetMap
  //For convenience we also maintain string descriptions of which data the subset mappings belong to  map(  map(string -> string)* -> string )  aka subsetMapDescr

  typedef std::unordered_map<std::string,size_t> paramIdxMap;
  typedef std::map<PiPiOperator, InnerParamMap> operatorSubsetMap;
  typedef std::map< InnerParamMap const*, std::string> subsetMapDescr;

#define CPB(A) typedef simultaneousFitBase::A A
#define COPY_BASE_TYPEDEFS						\
  CPB(CorrFuncJack); CPB(CorrFuncDJack); CPB(CorrFuncJackAllQ); CPB(CorrFuncDJackAllQ); \
  CPB(SimFitCorrFuncJack); CPB(SimFitCorrFuncDJack); CPB(InnerParamMap); CPB(Params); CPB(paramIdxMap); CPB(operatorSubsetMap); CPB(subsetMapDescr);
 
  virtual void load2ptFitParams(const std::vector<PiPiOperator> &operators, const InputParamArgs &iargs, const int nsample) = 0;

  static inline std::string opAmplitudeParamFmt(PiPiOperator op){
    switch(op){
    case PiPiOperator::PiPiGnd:
      return "Apipi%d";
    case PiPiOperator::PiPiExc:
      return "Apipi_exc_%d";
    case PiPiOperator::Sigma:
      return "Asigma%d";
    }
    assert(0);
    return "";
  }
  static inline std::string opDescr(PiPiOperator op){
    switch(op){
    case PiPiOperator::PiPiGnd:
      return "K->pipi(111)";
    case PiPiOperator::PiPiExc:
      return "K->pipi(311)";
    case PiPiOperator::Sigma:
      return "K->sigma";
    }
    assert(0);
    return "";
  }
  virtual std::vector<jackknifeDistribution<Params> > fit(const ResampledData<jackknifeDistributionD> &data_j,
							  const ResampledData<doubleJackknifeA0StorageType> &data_dj,
							  const std::vector<PiPiOperator> &operators,
							  const int Lt, const int tmin_k_op, const int tmin_op_snk, bool correlated) = 0;

  static void printWriteFitData(const std::vector<SimFitCorrFuncJack> &A0_sim_j,
			   const subsetMapDescr &pmap_descr){
    const int nQ = A0_sim_j.size();
    std::vector<std::vector<jackknifeDistributionD> > data_in_fit(nQ);
    std::ofstream data_in_fit_key("data_in_fit.key");

    std::cout << "Data in fit:\n";

    for(int q=0;q<nQ;q++){
      std::cout << "Q" << q+1 << std::endl;
      int eidx = 0;
      for(int i=0;i<A0_sim_j[q].size();i++){
	const SimFitCoordGen &c = A0_sim_j[q].coord(i);
	const jackknifeDistributionD &val = A0_sim_j[q].value(i);
	
	std::cout << printCoord(c, pmap_descr) << " " << val << std::endl;
	data_in_fit_key << "Q" << q+1 << " elem " << eidx << " " << pmap_descr.find(c.param_map)->second << " t=" << c.t << " tsep_K_snk=" << c.tsep_k_snk << std::endl;
	data_in_fit[q].push_back(val);
	eidx++;
      }
    }
    writeParamsStandard(data_in_fit,"data_in_fit.hdf5");
  }

  //Output vectors should be resized to number of matrix elements 
  template<typename DistributionType>
  static void generateSimData(std::vector<correlationFunction<SimFitCoordGen, DistributionType> > &A0_sim_r,
			      const ResampledData<DistributionType> &data_r,
			      const std::vector<PiPiOperator> &operators,   
			      const int tmin_k_op, const int tmin_op_snk,
			      const operatorSubsetMap &op_param_maps, const subsetMapDescr &pmap_descr){
    typedef std::vector< correlationFunction<amplitudeDataCoord, DistributionType> > CorrFuncAllQtype;

    int nops = operators.size();
    std::vector<CorrFuncAllQtype const*> to_include(nops);
    std::vector<InnerParamMap const *> incl_pmaps(nops);

    for(int opidx = 0; opidx < nops; opidx++){
      auto op = operators[opidx];
      assert(data_r.contains(op));
      to_include[opidx] = &data_r(op);
      incl_pmaps[opidx] = &op_param_maps.find(op)->second;
    }

    //Determine number of matrix elements from input data
    const int nQ = to_include[0]->size();
    for(int i=0;i<to_include.size();i++)
      assert(to_include[i]->size() == nQ); 
    
    //Get the data in the fit range and put in sim fit data containers
    A0_sim_r.resize(nQ);

    for(int q=0;q<nQ;q++){
      for(int o=0;o<nops;o++){
	const auto &data_r = (*to_include[o])[q];

	for(int d=0;d<data_r.size();d++){
	  const int t = int(data_r.coord(d).t);
	  const int tsep_k_snk = data_r.coord(d).tsep_k_pi;
	  const int tsep_op_snk = tsep_k_snk - t;
	  if(t <= tsep_k_snk && t >= tmin_k_op && tsep_op_snk >= tmin_op_snk){
	    SimFitCoordGen c(t, tsep_k_snk, incl_pmaps[o]);
	    A0_sim_r[q].push_back(c, data_r.value(d));
	  }
	}
      }
    }
  }

  template<typename FitPolicies>
  static void runfit(std::vector<typename FitPolicies::FitParameterDistribution> &params,
		     const std::vector<SimFitCorrFuncJack> &A0_sim_j,
		     const std::vector<SimFitCorrFuncDJack> &A0_sim_dj,
		     const subsetMapDescr &pmap_descr,
		     const typename FitPolicies::baseFitFunc &fitfunc,
		     const std::vector<int> &freeze_params,
		     const typename FitPolicies::FitParameterDistribution &freeze_vals,
		     const typename FitPolicies::baseFitFunc::ParameterType &guess,
		     const int nsample, const bool correlated){
    const int nQ = A0_sim_j.size(); assert(A0_sim_dj.size() == nQ);

    params.resize(nQ);

    std::vector<jackknifeDistributionD> chisq(nQ, jackknifeDistributionD(nsample)), 
      chisq_per_dof(nQ, jackknifeDistributionD(nsample)), 
      pvalue(nQ, jackknifeDistributionD(nsample));

    for(int q=0;q<nQ;q++){
      params[q] = jackknifeDistribution<Params>(nsample, guess);

      std::cout << "Performing " << q+1 << " fit" << std::endl;
      fitter<FitPolicies> fit;
      fit.importFitFunc(fitfunc);
      fit.freeze(freeze_params, freeze_vals);
      
      importCostFunctionParameters<correlatedFitPolicy, FitPolicies> import(fit, A0_sim_dj[q]);
      if(!correlated) import.setUncorrelated();
          
      int ndof;
      fit.fit(params[q], chisq[q], chisq_per_dof[q], ndof, A0_sim_j[q]);

      pvalue[q] = jackknifeDistributionD(nsample, [&](const int s){ return chiSquareDistribution::pvalue(ndof, chisq[q].sample(s)); });

      std::cout << "Q" << q+1 << std::endl;
      std::cout << "Params:\n" << params[q] << std::endl;
      std::cout << "Chisq: " << chisq[q] << std::endl;
      std::cout << "Chisq/dof: " << chisq_per_dof[q] << std::endl;
      std::cout << "p-value: " << pvalue[q] << std::endl;

      std::cout << "Analysis of contributions to chi^2" << std::endl;
      
      analyzeChisqFF<typename FitPolicies::baseFitFunc>(A0_sim_j[q],params[q],fitfunc,pmap_descr);
    }  
    for(int q=0;q<nQ;q++){
      std::cout << "Q" << q+1 << std::endl;
      std::cout << "Params:\n" << params[q] << std::endl;
      std::cout << "Chisq: " << chisq[q] << std::endl;
      std::cout << "Chisq/dof: " << chisq_per_dof[q] << std::endl;
      std::cout << "p-value: " << pvalue[q] << std::endl;
    }
    writeParamsStandard(params, "params.hdf5");
    writeParamsStandard(chisq, "chisq.hdf5");
    writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");
    writeParamsStandard(pvalue, "pvalue.hdf5");
  }

  static void freeze(std::vector<int> &freeze_params,
	      jackknifeDistribution<taggedValueContainer<double,std::string> > &freeze_vals,
	      const std::string &pname, const jackknifeDistributionD &pval, const std::unordered_map<std::string,size_t> &param_idx_map){
    auto it = param_idx_map.find(pname); assert(it != param_idx_map.end()); 
    freeze_params.push_back(it->second);
    int nsample = pval.size();
    for(int s=0;s<nsample;s++) freeze_vals.sample(s)(pname) = pval.sample(s);
  }

  virtual ~simultaneousFitBase(){}
};

class simultaneousFitMultiState: public simultaneousFitBase{
public:
  COPY_BASE_TYPEDEFS;    
    
protected:
  int nstate;

  //Values for freezing
  bool loaded_frzparams;
  jackknifeDistributionD mK;
  jackknifeDistributionD cK;
  std::vector<jackknifeDistributionD> E;
  std::map<PiPiOperator, std::vector<jackknifeDistributionD> > coeffs; //for each operator, the frozen amplitudes (0=gnd state, 1=exc state, ...)
public:

  void load2ptFitParams(const std::vector<PiPiOperator> &operators, const InputParamArgs &iargs, const int nsample){
    { //Load kaon parameters
      std::vector<jackknifeDistributionD> p;
      readParamsStandard(p,  iargs.kaon2pt_fit_result);
      mK = p[iargs.idx_mK];
      cK = sqrt( p[iargs.idx_cK] );
    }

    const std::vector<int> pipi_gnd_indices = { iargs.idx_coeff_pipi_state0, iargs.idx_coeff_pipi_state1, iargs.idx_coeff_pipi_state2 };
    const std::vector<int> pipi_exc_indices = { iargs.idx_coeff_pipi_exc_state0, iargs.idx_coeff_pipi_exc_state1, iargs.idx_coeff_pipi_exc_state2 };
    const std::vector<int> sigma_indices = { iargs.idx_coeff_sigma_state0, iargs.idx_coeff_sigma_state1, iargs.idx_coeff_sigma_state2 };    
    const std::vector<int> E_indices = { iargs.idx_E0, iargs.idx_E1, iargs.idx_E2 };

    { //Load coefficients and energies
      const double scale = sqrt(iargs.pipi_sigma_sim_fit_Ascale);
      std::vector<jackknifeDistributionD> p;
      readParamsStandard(p, iargs.pipi_sigma_sim_fit_result);
      for(int i=0;i<p.size();i++) assert(p[i].size() == nsample);

      //(0=gnd state, 1=exc state, ...)
      if(doOp(PiPiOperator::PiPiGnd, operators)){
	for(int s=0;s<nstate;s++) coeffs[PiPiOperator::PiPiGnd].push_back( p[pipi_gnd_indices[s]] * scale );
      }
      if(doOp(PiPiOperator::PiPiExc, operators)){
	for(int s=0;s<nstate;s++) coeffs[PiPiOperator::PiPiExc].push_back( p[pipi_exc_indices[s]] * scale );
      }
      if(doOp(PiPiOperator::Sigma, operators)){
	for(int s=0;s<nstate;s++) coeffs[PiPiOperator::Sigma].push_back( p[sigma_indices[s]] * scale );
      }
      
      //Energies
      for(int s=0;s<nstate;s++)
	E.push_back( p[E_indices[s]] );
    }

    std::cout << "cK = " << cK << std::endl;
    std::cout << "mK = " << mK << std::endl;
    for(int s=0;s<nstate;s++)
      std::cout << "E" << s <<" = " << E[s] << std::endl;
    loaded_frzparams = true;
  }
  
  simultaneousFitMultiState(const int nstate): loaded_frzparams(false), nstate(nstate){
    assert(nstate <= 3);
  }


 void constructParameterMaps(paramIdxMap &param_idx_map,
			     operatorSubsetMap &op_param_maps,
			     subsetMapDescr &pmap_descr,
			     const std::vector<PiPiOperator> &operators,
			     const bool include_kaon_params = true) const{
    int nop = operators.size();

    //Construct parameter names for <op|state>
    std::string A_op_state[nop][nstate];
    for(int o=0;o<nop;o++)
      for(int s=0;s<nstate;s++)
	A_op_state[o][s] = stringize(opAmplitudeParamFmt(operators[o]).c_str(),s);

    //Construct parameter names for state energies and matrix elements
    std::string E[nstate], M[nstate];
    for(int s=0;s<nstate;s++){
      E[s] = stringize("E%d",s);
      M[s] = stringize("M%d",s);
    }

    //Setup the mapping of global parameter to global index
    int idx = 0;
    if(include_kaon_params){
      param_idx_map["AK"] = idx++;
      param_idx_map["mK"] = idx++;
    }
    for(int o=0;o<nop;o++)
      for(int s=0;s<nstate;s++)
	param_idx_map[A_op_state[o][s]] = idx++;

    for(int s=0;s<nstate;s++)
      param_idx_map[E[s]] = idx++;
    for(int s=0;s<nstate;s++)
      param_idx_map[M[s]] = idx++;

    //For each operator, setup the mapping of local parameter to global parameter
    for(int o=0;o<nop;o++){
      auto & pmap = op_param_maps[operators[o]];

      for(int s=0;s<nstate;s++)
	pmap[stringize("Asnk%d",s)] = A_op_state[o][s];

      //The rest of the parameters are 1->1 mapping
      if(include_kaon_params){
	pmap["AK"] = "AK";
	pmap["mK"] = "mK";
      }
      for(int s=0;s<nstate;s++){
	pmap[E[s]] = E[s];
	pmap[M[s]] = M[s];
      }
    }

    //Finally, associate with each local->global parameter mapping a description
    for(int o=0;o<nop;o++){
      auto & pmap = op_param_maps[operators[o]];
      pmap_descr[&pmap] = opDescr(operators[o]); 
    }
  }


  std::vector<jackknifeDistribution<Params> > fit(const ResampledData<jackknifeDistributionD> &data_j,
						  const ResampledData<doubleJackknifeA0StorageType> &data_dj,
						  const std::vector<PiPiOperator> &operators,
						  const int Lt, const int tmin_k_op, const int tmin_op_snk, bool correlated){
    
    assert(loaded_frzparams);
    
    const int nsample = mK.size();
    
    //Get the mappings for the fit parameters
    paramIdxMap* param_idx_map_ptr = new paramIdxMap; //ensure it sticks around until end of execution
    auto &param_idx_map = *param_idx_map_ptr;
    
    operatorSubsetMap op_param_maps;
    subsetMapDescr pmap_descr;

    constructParameterMaps(param_idx_map, op_param_maps, pmap_descr, operators);

    std::vector<SimFitCorrFuncJack> A0_sim_j;
    std::vector<SimFitCorrFuncDJack> A0_sim_dj;
    generateSimData(A0_sim_j, data_j, operators, tmin_k_op, tmin_op_snk, op_param_maps, pmap_descr);
    generateSimData(A0_sim_dj, data_dj, operators, tmin_k_op, tmin_op_snk, op_param_maps, pmap_descr);
 
    //Setup guess
    Params guess(&param_idx_map);
    for(int i=0;i<guess.size();i++) guess(i) = 1.;
    for(int s=0;s<nstate;s++)
      guess(stringize("M%d",s)) = 0.5;

    typedef FitSimGenMultiState FitFunc;

    //Setup fitfunc
    FitFunc fitfunc(param_idx_map.size(), nstate);
    typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, correlatedFitPolicy>::type FitPolicies;
       
    //Setup frozen fit params
    std::vector<int> freeze_params;
    jackknifeDistribution<Params> freeze_vals(nsample, Params(&param_idx_map));
  
    freeze(freeze_params, freeze_vals, "AK", cK, param_idx_map);
    freeze(freeze_params, freeze_vals, "mK", mK, param_idx_map);
    for(int s=0; s<nstate; s++) 
      freeze(freeze_params, freeze_vals, stringize("E%d",s), E[s], param_idx_map);

    for(int opidx = 0; opidx < operators.size(); opidx++){
      auto op = operators[opidx];

      for(int s=0; s<nstate; s++) 
	freeze(freeze_params, freeze_vals, stringize(opAmplitudeParamFmt(op).c_str(),s), coeffs[op][s], param_idx_map);
    }

    //Run the actual fit
    std::vector<jackknifeDistribution<Params> > params;    
    runfit<FitPolicies>(params, A0_sim_j, A0_sim_dj, pmap_descr, fitfunc, freeze_params, freeze_vals, guess, nsample, correlated);
 
    //plotErrorWeightedData2expFlat(ktopipi_A0_all_j, ktopipi_exc_A0_all_j, ktosigma_A0_all_j, params, Lt, tmin_k_op, tmin_op_snk, fitfunc);

    return params;
  }
};

class simultaneousFitMultiStateWavg: public simultaneousFitMultiState{
public:
  COPY_BASE_TYPEDEFS;    
  
  simultaneousFitMultiStateWavg(const int nstate): simultaneousFitMultiState(nstate){}

  template<typename DistributionType>
  inline DistributionType weightedAvg(const correlationFunction<SimFitCoordGen, DistributionType> &from, const std::vector<int> &idx){
    std::vector<DistributionType const*> towavg(idx.size());
    for(int i=0;i<idx.size();i++)
      towavg[i] = &from.value(idx[i]);
    return CPSfit::weightedAvg(towavg);
  }

  //Divide out the kaon amplitude and time dependence
  template<typename DataType>
  std::vector<DataType> removeKaonDependence(const std::vector<DataType> &from) const{
    int nsample = from[0].value(0).size();
    jackknifeDistributionD one(nsample, 1.);

    std::vector<DataType> out(from);
    for(int q=0;q<out.size();q++){
      for(int d=0;d<out[q].size();d++){
	int tK_op = (int)out[q].coord(d).t;
	jackknifeDistributionD fac = one/(this->cK * exp(-this->mK * tK_op));
	for(int s=0;s<nsample;s++)
	  out[q].value(d).sample(s) = out[q].value(d).sample(s) * fac.sample(s); //note: for double-jackknife we use the same value for all inner samples - this is a minor effect
      }
    }
    return out;
  }

  std::vector<jackknifeDistribution<Params> > fit(const ResampledData<jackknifeDistributionD> &data_j,
						  const ResampledData<doubleJackknifeA0StorageType> &data_dj,
						  const std::vector<PiPiOperator> &operators,
						  const int Lt, const int tmin_k_op, const int tmin_op_snk, bool correlated){
    
    assert(this->loaded_frzparams);
    
    const int nsample = this->mK.size();
    
    //Get the mappings for the fit parameters
    paramIdxMap* param_idx_map_ptr = new paramIdxMap; //ensure it sticks around until end of execution
    auto &param_idx_map = *param_idx_map_ptr;
    
    operatorSubsetMap op_param_maps;
    subsetMapDescr pmap_descr;

    constructParameterMaps(param_idx_map, op_param_maps, pmap_descr, operators, false); //don't include kaon parameters as we divide them out

    std::vector<SimFitCorrFuncJack> A0_sim_j;
    std::vector<SimFitCorrFuncDJack> A0_sim_dj;
    generateSimData(A0_sim_j, data_j, operators, tmin_k_op, tmin_op_snk, op_param_maps, pmap_descr);
    generateSimData(A0_sim_dj, data_dj, operators, tmin_k_op, tmin_op_snk, op_param_maps, pmap_descr);

    const int nQ = A0_sim_j.size(); //could be 7 or 10 depending on basis
  
    std::vector<SimFitCorrFuncJack> A0_sim_j_noK;
    std::vector<SimFitCorrFuncDJack> A0_sim_dj_noK;
    A0_sim_j_noK = removeKaonDependence(A0_sim_j);
    A0_sim_dj_noK = removeKaonDependence(A0_sim_dj);

    //Do the weighted average
    std::vector<SimFitCorrFuncJack> A0_sim_j_wavg(nQ);
    std::vector<SimFitCorrFuncDJack> A0_sim_dj_wavg(nQ);

    for(int q=0;q<nQ;q++){
      //Collect data indices of values we wish to wavg
      std::map<std::pair<int,InnerParamMap const*>, std::vector<int> > to_avg;  //key is (top_snk, snk_op_pmap_descr)
      for(int d=0;d<A0_sim_j_noK[q].size();d++){
	auto const &c = A0_sim_j_noK[q].coord(d);
	int tK_op = (int)c.t,   top_snk = (int)c.tsep_k_snk - tK_op;
	to_avg[{top_snk, c.param_map}].push_back(d);
      }
      
      //Perform the weighted average and insert into fit data containers
      for(auto it = to_avg.begin(); it != to_avg.end(); it++){
	SimFitCoordGen c(it->first.first, -1,  it->first.second);
	const std::vector<int>& idxv = it->second;	

	A0_sim_j_wavg[q].push_back(c, weightedAvg(A0_sim_j_noK[q], idxv));
	A0_sim_dj_wavg[q].push_back(c, weightedAvg(A0_sim_dj_noK[q], idxv));

	//Print some useful information
	{
	  const jackknifeDistributionD &wavg = A0_sim_j_wavg[q].value(A0_sim_j_wavg[q].size()-1);
	  std::cout << "Q" << q+1 << " weighted avg with descr " << pmap_descr.find(c.param_map)->second << " and top_snk = " << c.t << " : " << wavg << std::endl;
	  std::cout << "Data included" << std::endl;
	  for(int i=0;i<idxv.size();i++){ 
	    jackknifeDistributionD wavg_diff = A0_sim_j_noK[q].value(idxv[i]) - wavg;
	    std::cout << printCoord(A0_sim_j_noK[q].coord(idxv[i]), pmap_descr) << " " << A0_sim_j_noK[q].value(idxv[i]) << " (diff from wavg: " << wavg_diff << ")" << std::endl;
	  }
	}
      }
    }//q loop


    //Setup guess
    Params guess(&param_idx_map);
    for(int i=0;i<guess.size();i++) guess(i) = 1.;
    for(int s=0;s<nstate;s++)
      guess(stringize("M%d",s)) = 0.5;

    typedef FitSimGenMultiStateWavg FitFunc;

    //Setup fitfunc
    FitFunc fitfunc(param_idx_map.size(), nstate);
    typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, correlatedFitPolicy>::type FitPolicies;
       
    //Setup frozen fit params
    std::vector<int> freeze_params;
    jackknifeDistribution<Params> freeze_vals(nsample, Params(&param_idx_map));
  
    for(int s=0; s<nstate; s++) 
      freeze(freeze_params, freeze_vals, stringize("E%d",s), E[s], param_idx_map);

    for(int opidx = 0; opidx < operators.size(); opidx++){
      auto op = operators[opidx];

      for(int s=0; s<nstate; s++) 
	freeze(freeze_params, freeze_vals, stringize(opAmplitudeParamFmt(op).c_str(),s), coeffs[op][s], param_idx_map);
    }

    //Run the actual fit
    std::vector<jackknifeDistribution<Params> > params;    
    runfit<FitPolicies>(params, A0_sim_j_wavg, A0_sim_dj_wavg, pmap_descr, fitfunc, freeze_params, freeze_vals, guess, nsample, correlated);
 
    plotErrorWeightedDataNexpFlat(data_j, operators, fitfunc, params, this->mK, this->cK, Lt, tmin_k_op, tmin_op_snk);

    return params;
  }
};







simultaneousFitBase* getFitter(const SimFitFunction ff, const int nstate){
  switch(ff){
  case SimFitFunction::MultiState:
    return new simultaneousFitMultiState(nstate);
  case SimFitFunction::MultiStateWavg:
    return new simultaneousFitMultiStateWavg(nstate);
  default:
    assert(0);
  }
};


#undef COPY_BASE_TYPEDEFS

CPSFIT_END_NAMESPACE

#endif
