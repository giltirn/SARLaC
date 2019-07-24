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


//Common functionality for bootstrap/jackknife
struct simultaneousFitCommon{
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

#define CPC(A) typedef simultaneousFitCommon::A A
#define COPY_COMMON_TYPEDEFS						\
  CPC(InnerParamMap); CPC(Params); CPC(paramIdxMap); CPC(operatorSubsetMap); CPC(subsetMapDescr);

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

  template<typename DistributionType>
  static void printWriteFitData(const std::vector<correlationFunction<SimFitCoordGen, DistributionType> > &A0_sim_j,
				const subsetMapDescr &pmap_descr){
    const int nQ = A0_sim_j.size();
    std::vector<std::vector<DistributionType> > data_in_fit(nQ);
    std::ofstream data_in_fit_key("data_in_fit.key");

    std::cout << "Data in fit:\n";

    for(int q=0;q<nQ;q++){
      std::cout << "Q" << q+1 << std::endl;
      int eidx = 0;
      for(int i=0;i<A0_sim_j[q].size();i++){
	const SimFitCoordGen &c = A0_sim_j[q].coord(i);
	const DistributionType &val = A0_sim_j[q].value(i);
	
	std::cout << printCoord(c, pmap_descr) << " " << val << std::endl;
	data_in_fit_key << "Q" << q+1 << " elem " << eidx << " " << pmap_descr.find(c.param_map)->second << " t=" << c.t << " tsep_K_snk=" << c.tsep_k_snk << std::endl;
	data_in_fit[q].push_back(val);
	eidx++;
      }
    }
    writeParamsStandard(data_in_fit,"data_in_fit.hdf5");
  }

  template< template<typename, template<typename> class> class DistributionType, template<typename> class V>
  static void freeze(std::vector<int> &freeze_params,
		     DistributionType<taggedValueContainer<double,std::string> ,V> &freeze_vals,
		     const std::string &pname, const DistributionType<double, V> &pval, const std::unordered_map<std::string,size_t> &param_idx_map){
    typedef iterate< DistributionType<taggedValueContainer<double,std::string> ,V> > iter_p;
    typedef iterate< DistributionType<double, V> > iter_d;

    auto it = param_idx_map.find(pname); assert(it != param_idx_map.end()); 
    freeze_params.push_back(it->second);
    for(int s=0;s<iter_d::size(pval);s++) iter_p::at(s,freeze_vals)(pname) = iter_d::at(s,pval);
  }
};


template< template<typename, template<typename> class> class DistributionType > 
struct ResampledDataContainers{};

template<>
struct ResampledDataContainers<jackknifeDistribution>{
  const ResampledData<jackknifeDistributionD> &data_j;
  const ResampledData<doubleJackknifeA0StorageType> &data_dj;
  const ResampledData<blockDoubleJackknifeA0StorageType> &data_bdj;
  
  ResampledDataContainers(const ResampledData<jackknifeDistributionD> &data_j,
			  const ResampledData<doubleJackknifeA0StorageType> &data_dj,
			  const ResampledData<blockDoubleJackknifeA0StorageType> &data_bdj): data_j(data_j), data_dj(data_dj), data_bdj(data_bdj){}
};

template< template<typename, template<typename> class> class DistributionType > 
struct SimFitDataContainers{};

template<>
struct SimFitDataContainers<jackknifeDistribution>{
  COPY_COMMON_TYPEDEFS;

  typedef correlationFunction<SimFitCoordGen, jackknifeDistributionD> SimFitCorrFuncJack;
  typedef correlationFunction<SimFitCoordGen, doubleJackknifeDistributionD> SimFitCorrFuncDJack;
  typedef correlationFunction<SimFitCoordGen, blockDoubleJackknifeDistributionD> SimFitCorrFuncBDJack;

  std::vector<SimFitCorrFuncJack> A0_sim_j;
  std::vector<SimFitCorrFuncDJack> A0_sim_dj;
  std::vector<SimFitCorrFuncBDJack> A0_sim_bdj;
  
  int getNq() const{
    int nq = A0_sim_j.size();
    return nq;
  }

  template<typename FitPolicies>
  void generateCovarianceMatrix(importCostFunctionParameters<correlatedFitPolicy, FitPolicies> &import, fitter<FitPolicies> &fit, const CovarianceMatrix covariance_matrix, const int q) const{
    bool do_dj = covariance_matrix == CovarianceMatrix::Regular;
    bool do_bdj = covariance_matrix == CovarianceMatrix::Block;

    if(do_dj) import.import(fit, A0_sim_dj[q]);
    if(do_bdj) import.import(fit, A0_sim_bdj[q]);
  }

  int getDistributionInitializer() const{
    return A0_sim_j[0].value(0).size();
  }

  void generateSimData(const ResampledDataContainers<jackknifeDistribution> &fit_data,
		       const std::vector<PiPiOperator> &operators, const int tmin_k_op, const int tmin_op_snk,
		       const operatorSubsetMap &op_param_maps, const subsetMapDescr &pmap_descr, const CovarianceMatrix covariance_matrix){
    bool do_dj = covariance_matrix == CovarianceMatrix::Regular;
    bool do_bdj = covariance_matrix == CovarianceMatrix::Block;
    simultaneousFitCommon::generateSimData(A0_sim_j, fit_data.data_j, operators, tmin_k_op, tmin_op_snk, op_param_maps, pmap_descr);
    if(do_dj) simultaneousFitCommon::generateSimData(A0_sim_dj, fit_data.data_dj, operators, tmin_k_op, tmin_op_snk, op_param_maps, pmap_descr);
    if(do_bdj) simultaneousFitCommon::generateSimData(A0_sim_bdj, fit_data.data_bdj, operators, tmin_k_op, tmin_op_snk, op_param_maps, pmap_descr);
  }
  
  const SimFitCorrFuncJack & getFitData(const int q) const{ return A0_sim_j[q]; }

  //Func:  void [](auto &corrfunc, const int q)
  template<typename Func, typename CorrelationFunctionType>
  void applyFunctionToInternalCorrFunc(const Func &func, CorrelationFunctionType &corrfunc){
    for(int q=0;q<corrfunc.size();q++)
      func(corrfunc[q],q);
  }

  //Func:  void [](auto &corrfunc, const int q)
  template<typename Func>
  void applyFunctionToCorrFunc(const Func &func){
    applyFunctionToInternalCorrFunc(func, A0_sim_j);
    applyFunctionToInternalCorrFunc(func, A0_sim_dj);
    applyFunctionToInternalCorrFunc(func, A0_sim_bdj);
  }

};

template< template<typename, template<typename> class> class DistributionType > 
struct simultaneousFitBase: public simultaneousFitCommon{
  COPY_COMMON_TYPEDEFS;

  virtual void load2ptFitParams(const std::vector<PiPiOperator> &operators, const InputParamArgs &iargs, const int dist_size) = 0;

  virtual std::vector<DistributionType<Params, basic_vector> > fit(const ResampledDataContainers<DistributionType> &fit_data,
								   const std::vector<PiPiOperator> &operators,
								   const int Lt, const int tmin_k_op, const int tmin_op_snk, bool correlated,
								   const CovarianceMatrix covariance_matrix) = 0;

  template<typename FitFunc>
  static void runfit(std::vector<DistributionType<Params, basic_vector> > &params,
		     const SimFitDataContainers<DistributionType> &fit_data,
		     const subsetMapDescr &pmap_descr,
		     const FitFunc &fitfunc,
		     const std::vector<int> &freeze_params,
		     const DistributionType<Params, basic_vector> &freeze_vals,
		     const Params &guess,
		     const bool correlated, const CovarianceMatrix covariance_matrix){

    typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, correlatedFitPolicy>::type FitPolicies;

    typedef DistributionType<double, basic_vector> DistributionTypeD;
    const int nQ = fit_data.getNq();
    params.resize(nQ);

    auto init = fit_data.getDistributionInitializer();

    std::vector<DistributionTypeD> chisq(nQ, DistributionTypeD(init)), 
      chisq_per_dof(nQ, DistributionTypeD(init)), 
      pvalue(nQ, DistributionTypeD(init));

    for(int q=0;q<nQ;q++){
      params[q] = DistributionType<Params, basic_vector>(init, guess);

      std::cout << "Performing " << q+1 << " fit" << std::endl;
      fitter<FitPolicies> fit;
      fit.importFitFunc(fitfunc);
      fit.freeze(freeze_params, freeze_vals);
      
      importCostFunctionParameters<correlatedFitPolicy, FitPolicies> import;
      fit_data.generateCovarianceMatrix(import, fit, covariance_matrix, q);

      if(!correlated) import.setUncorrelated();
          
      int ndof;
      fit.fit(params[q], chisq[q], chisq_per_dof[q], ndof, fit_data.getFitData(q));

      pvalue[q] = DistributionTypeD(init);
      for(int i=0; i<iterate<DistributionTypeD>::size(pvalue[q]); i++) 
	iterate<DistributionTypeD>::at(i, pvalue[q]) = chiSquareDistribution::pvalue(ndof, iterate<DistributionTypeD>::at(i, chisq[q]) );

      std::cout << "Q" << q+1 << std::endl;
      std::cout << "Params:\n" << params[q] << std::endl;
      std::cout << "Chisq: " << chisq[q] << std::endl;
      std::cout << "Chisq/dof: " << chisq_per_dof[q] << std::endl;
      std::cout << "p-value(chi^2): " << pvalue[q] << std::endl;

      // std::cout << "Analysis of contributions to chi^2" << std::endl;      
      // analyzeChisqFF<typename FitPolicies::baseFitFunc>(A0_sim_j[q],params[q],fitfunc,pmap_descr);
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


  virtual ~simultaneousFitBase(){}
};

template< template<typename, template<typename> class> class DistributionType > 
class simultaneousFitMultiState: public simultaneousFitBase<DistributionType>{
public:
  COPY_COMMON_TYPEDEFS;    
    
protected:
  int nstate;

  //Values for freezing
  bool loaded_frzparams;
  
  typedef DistributionType<double, basic_vector> DistributionTypeD;
  DistributionTypeD mK;
  DistributionTypeD cK;
  std::vector<DistributionTypeD> E;
  std::map<PiPiOperator, std::vector<DistributionTypeD> > coeffs; //for each operator, the frozen amplitudes (0=gnd state, 1=exc state, ...)
public:

  typedef FitSimGenMultiState FitFunc;

  void load2ptFitParams(const std::vector<PiPiOperator> &operators, const InputParamArgs &iargs, const int dist_size){
    { //Load kaon parameters
      std::vector<DistributionTypeD> p;
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
      std::vector<DistributionTypeD> p;
      readParamsStandard(p, iargs.pipi_sigma_sim_fit_result);
      for(int i=0;i<p.size();i++) assert(p[i].size() == dist_size);

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
	A_op_state[o][s] = stringize(this->opAmplitudeParamFmt(operators[o]).c_str(),s);

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
      pmap_descr[&pmap] = this->opDescr(operators[o]); 
    }
  }



  std::vector<DistributionType<Params, basic_vector> > fit(const ResampledDataContainers<DistributionType> &fit_data,
							   const std::vector<PiPiOperator> &operators,
							   const int Lt, const int tmin_k_op, const int tmin_op_snk, bool correlated, 
							   const CovarianceMatrix covariance_matrix){
    
    assert(loaded_frzparams);
    
    //Get the mappings for the fit parameters
    paramIdxMap* param_idx_map_ptr = new paramIdxMap; //ensure it sticks around until end of execution
    auto &param_idx_map = *param_idx_map_ptr;
    
    operatorSubsetMap op_param_maps;
    subsetMapDescr pmap_descr;

    constructParameterMaps(param_idx_map, op_param_maps, pmap_descr, operators);

    SimFitDataContainers<DistributionType> simfit_data;
    simfit_data.generateSimData(fit_data, operators, tmin_k_op, tmin_op_snk, op_param_maps, pmap_descr, covariance_matrix);
 
    //Setup guess
    Params guess(&param_idx_map);
    for(int i=0;i<guess.size();i++) guess(i) = 1.;
    for(int s=0;s<nstate;s++)
      guess(stringize("M%d",s)) = 0.5;

    //Setup fitfunc
    FitFunc fitfunc(param_idx_map.size(), nstate);
    auto init = simfit_data.getDistributionInitializer();
      
    //Setup frozen fit params
    std::vector<int> freeze_params;
    DistributionType<Params, basic_vector> freeze_vals(init, Params(&param_idx_map));
  
    this->freeze(freeze_params, freeze_vals, "AK", cK, param_idx_map);
    this->freeze(freeze_params, freeze_vals, "mK", mK, param_idx_map);
    for(int s=0; s<nstate; s++) 
      this->freeze(freeze_params, freeze_vals, stringize("E%d",s), E[s], param_idx_map);

    for(int opidx = 0; opidx < operators.size(); opidx++){
      auto op = operators[opidx];

      for(int s=0; s<nstate; s++) 
	this->freeze(freeze_params, freeze_vals, stringize(this->opAmplitudeParamFmt(op).c_str(),s), coeffs[op][s], param_idx_map);
    }

    //Run the actual fit
    std::vector<DistributionType<Params, basic_vector> > params;    
  
    this->runfit(params, simfit_data, pmap_descr, fitfunc, freeze_params, freeze_vals, guess, correlated, covariance_matrix);
 
    //plotErrorWeightedData2expFlat(ktopipi_A0_all_j, ktopipi_exc_A0_all_j, ktosigma_A0_all_j, params, Lt, tmin_k_op, tmin_op_snk, fitfunc);

    return params;
  }
};


template< template<typename, template<typename> class> class DistributionType > 
class simultaneousFitMultiStateWavg: public simultaneousFitMultiState<DistributionType>{
public:
  COPY_COMMON_TYPEDEFS;    

  typedef DistributionType<double, basic_vector> DistributionTypeD;

  typedef FitSimGenMultiStateWavg FitFunc;
  
  simultaneousFitMultiStateWavg(const int nstate): simultaneousFitMultiState<DistributionType>(nstate){}

  template<typename T>
  static inline T weightedAvg(const correlationFunction<SimFitCoordGen, T> &from, const std::vector<int> &idx){
    std::vector<T const*> towavg(idx.size());
    for(int i=0;i<idx.size();i++)
      towavg[i] = &from.value(idx[i]);
    return CPSfit::weightedAvg(towavg);
  }

  void weightedAverage(SimFitDataContainers<DistributionType> &fit_data, const subsetMapDescr &pmap_descr) const{
    //Collect data indices of values we wish to wavg
    std::vector< std::map<std::pair<int,InnerParamMap const*>, std::vector<int> > > toavg_q(fit_data.getNq());

    for(int q=0;q<toavg_q.size();q++){
      std::map<std::pair<int,InnerParamMap const*>, std::vector<int> > &to_avg = toavg_q[q];  //key is (top_snk, snk_op_pmap_descr)
      for(int d=0;d< fit_data.getFitData(q).size();d++){
	auto const &c = fit_data.getFitData(q).coord(d);
	int tK_op = (int)c.t,   top_snk = (int)c.tsep_k_snk - tK_op;
	to_avg[{top_snk, c.param_map}].push_back(d);

	//Print some useful information
	for(auto it = to_avg.begin(); it != to_avg.end(); it++){
	  SimFitCoordGen c(it->first.first, -1,  it->first.second);
	  const std::vector<int>& idxv = it->second;	       

	  auto wavg = weightedAvg(fit_data.getFitData(q), idxv);
	  
	  std::cout << "Q" << q+1 << " weighted avg with descr " << pmap_descr.find(c.param_map)->second << " and top_snk = " << c.t << " : " << wavg << std::endl;
	  std::cout << "Data included" << std::endl;
	  for(int i=0;i<idxv.size();i++){ 
	    typename std::decay<decltype(fit_data.getFitData(q).value(idxv[i]))>::type 
	      wavg_diff = fit_data.getFitData(q).value(idxv[i]) - wavg;
	    std::cout << printCoord(fit_data.getFitData(q).coord(idxv[i]), pmap_descr) << " " 
		      << fit_data.getFitData(q).value(idxv[i]) 
		      << " (diff from wavg: " << wavg_diff << ")" << std::endl;
	  }
	}
      }
    }

    //Apply the weighted avg to all correlation functions in the container
    fit_data.applyFunctionToCorrFunc(
      [&](auto &corrfunc, const int q){
	typedef typename std::decay<decltype(corrfunc)>::type CorrFuncType;
	CorrFuncType out;

	const std::map<std::pair<int,InnerParamMap const*>, std::vector<int> > &to_avg = toavg_q[q];
	
	//Perform the weighted average and insert into fit data containers
	for(auto it = to_avg.begin(); it != to_avg.end(); it++){
	  SimFitCoordGen c(it->first.first, -1,  it->first.second);
	  const std::vector<int>& idxv = it->second;	       
	  out.push_back(c, weightedAvg(corrfunc, idxv));
	}

	corrfunc = out;
      });
  }

  //Divide out the kaon amplitude and time dependence
  void removeKaonDependence(SimFitDataContainers<DistributionType> &fit_data) const{
    fit_data.applyFunctionToCorrFunc(
      [&](auto &corrfunc, const int q){
	for(int d=0;d<corrfunc.size();d++){
	  const SimFitCoordGen &c = corrfunc.coord(d);
	  for(int s=0;s<corrfunc.value(d).size();s++){ 
	    //do for each *outer* sample. ValueType is either double or a distribution. We use the same value for each inner sample
	    double cK_s = this->cK.sample(s);
	    double mK_s = this->mK.sample(s);
	    int tK_op = (int)c.t;
	    corrfunc.value(d).sample(s) = corrfunc.value(d).sample(s) / (cK_s * exp(-mK_s * tK_op)); 
	  }
	}
      }
   );
  }   

  std::vector<DistributionType<Params, basic_vector> > fit(const ResampledDataContainers<DistributionType> &fit_data,
							   const std::vector<PiPiOperator> &operators,
							   const int Lt, const int tmin_k_op, const int tmin_op_snk, 
							   bool correlated, const CovarianceMatrix covariance_matrix){
    assert(this->loaded_frzparams);
    
    //Get the mappings for the fit parameters
    paramIdxMap* param_idx_map_ptr = new paramIdxMap; //ensure it sticks around until end of execution
    auto &param_idx_map = *param_idx_map_ptr;
    
    operatorSubsetMap op_param_maps;
    subsetMapDescr pmap_descr;

    this->constructParameterMaps(param_idx_map, op_param_maps, pmap_descr, operators, false); //don't include kaon parameters as we divide them out

    SimFitDataContainers<DistributionType> simfit_data;
    simfit_data.generateSimData(fit_data, operators, tmin_k_op, tmin_op_snk, op_param_maps, pmap_descr, covariance_matrix);

    const int nQ = simfit_data.getNq(); //could be 7 or 10 depending on basis
  
    removeKaonDependence(simfit_data);
    weightedAverage(simfit_data, pmap_descr);

    //Setup guess
    Params guess(&param_idx_map);
    for(int i=0;i<guess.size();i++) guess(i) = 1.;
    for(int s=0;s<this->nstate;s++)
      guess(stringize("M%d",s)) = 0.5;

    //Setup fitfunc
    FitFunc fitfunc(param_idx_map.size(), this->nstate);
       
    //Setup frozen fit params
    auto init = simfit_data.getFitData(0).value(0).getInitializer();

    std::vector<int> freeze_params;
    DistributionType<Params, basic_vector> freeze_vals(init, Params(&param_idx_map));
  
    for(int s=0; s<this->nstate; s++) 
      this->freeze(freeze_params, freeze_vals, stringize("E%d",s), this->E[s], param_idx_map);

    for(int opidx = 0; opidx < operators.size(); opidx++){
      auto op = operators[opidx];

      for(int s=0; s<this->nstate; s++) 
	this->freeze(freeze_params, freeze_vals, stringize(this->opAmplitudeParamFmt(op).c_str(),s), this->coeffs[op][s], param_idx_map);
    }


    //Run the actual fit
    std::vector<DistributionType<Params, basic_vector> > params;    
  
    this->runfit(params, simfit_data, pmap_descr, fitfunc, freeze_params, freeze_vals, guess, correlated, covariance_matrix);
 
    //plotErrorWeightedDataNexpFlat(data_j, operators, fitfunc, params, this->mK, this->cK, Lt, tmin_k_op, tmin_op_snk);

    return params;
  }
};


template< template<typename, template<typename> class> class DistributionType > 
simultaneousFitBase<DistributionType>* getFitter(const SimFitFunction ff, const int nstate){
  switch(ff){
  case SimFitFunction::MultiState:
    return new simultaneousFitMultiState<DistributionType>(nstate);
  case SimFitFunction::MultiStateWavg:
    return new simultaneousFitMultiStateWavg<DistributionType>(nstate);
  default:
    assert(0);
  }
};


template<template<typename, template<typename> class> class DistributionType, typename Params, template<typename> class V>
std::vector<std::vector<DistributionType<double,V> > > convert7basisTo10basis(const int nstate,
									      const std::vector<DistributionType<Params,V> > &params){		    

  std::cout << "Converting 7 basis results to 10 basis" << std::endl;
  typedef DistributionType<Params,V> DistributionTypeP;
  typedef iterate<DistributionTypeP> iter_p;
  typedef DistributionType<double,V> DistributionTypeD;
  typedef iterate<DistributionTypeD> iter_d;
  DistributionTypeD zero(params[0].getInitializer()); //will initialize correctly jackknife or bootstrap
  zeroit(zero);

  std::vector<std::vector<DistributionTypeD> > params_10(10, std::vector<DistributionTypeD>(nstate, zero));    
    
  struct InContainer{
    size_t s;
    size_t idx;
    const std::vector<DistributionTypeP> &d;
    double operator[](const int q) const{ return iter_p::at(s, d[q])(idx); }
    InContainer(size_t s, size_t idx, const std::vector<DistributionTypeP> &d): s(s), idx(idx), d(d){}
  };
  struct OutContainer{
    size_t s;
    size_t state;
    std::vector<std::vector<DistributionTypeD> > &d;
    double & operator[](const int q){ return iter_d::at(s,d[q][state]); }
    OutContainer(size_t s, size_t state, std::vector<std::vector<DistributionTypeD> > &d): s(s), state(state), d(d){}
  };

  auto const* tag_map = params[0].sample(0).getTagMap();
  assert(tag_map != NULL);

  for(int i=0;i<nstate;i++){
    auto it = tag_map->find(stringize("M%d",i));
    assert(it != tag_map->end());
    int Midx = it->second;
      
    for(int s=0;s<iter_d::size(zero);s++){
      InContainer in(s, Midx, params);
      OutContainer out(s, i, params_10);
      convert7to10(out,in);
    }
  }
  return params_10;
}




#undef COPY_BASE_TYPEDEFS
#undef COPY_COMMON_TYPEDEFS
#undef CPC
#undef CPB

CPSFIT_END_NAMESPACE

#endif
