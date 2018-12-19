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
  
#define COPY_BASE_TYPEDEFS				    \
  typedef simultaneousFitBase::CorrFuncJack CorrFuncJack;     \
  typedef simultaneousFitBase::CorrFuncDJack CorrFuncDJack;	    \
  typedef simultaneousFitBase::CorrFuncJackAllQ CorrFuncJackAllQ;     \
  typedef simultaneousFitBase::CorrFuncDJackAllQ CorrFuncDJackAllQ;	\
  typedef simultaneousFitBase::SimFitCorrFuncJack SimFitCorrFuncJack;	\
  typedef simultaneousFitBase::SimFitCorrFuncDJack SimFitCorrFuncDJack; \
  typedef simultaneousFitBase::InnerParamMap InnerParamMap;		\
  typedef simultaneousFitBase::Params Params
  
  virtual void load2ptFitParams(const std::vector<PiPiOperator> &operators, const InputParamArgs &iargs, const int nsample) = 0;
  virtual std::vector<jackknifeDistribution<Params> > fit(const ResampledData<jackknifeDistributionD> &data_j,
							  const ResampledData<doubleJackknifeA0StorageType> &data_dj,
							  const std::vector<PiPiOperator> &operators,
							  const int Lt, const int tmin_k_op, const int tmin_op_snk, bool correlated) = 0;


  //Output vectors should be resized to number of matrix elements 
  static void generateSimData(std::vector<SimFitCorrFuncJack> &A0_sim_j,
			      std::vector<SimFitCorrFuncDJack> &A0_sim_dj,
			      const std::vector<CorrFuncJackAllQ const* > &jacks,
			      const std::vector<CorrFuncDJackAllQ const* > &djacks,
			      const std::vector<InnerParamMap const *> &pmaps,
			      const int tmin_k_op, const int tmin_op_snk,
			      const std::map< InnerParamMap const*, std::string> &pmap_descr){
    int nops = jacks.size();
    assert(djacks.size() == nops && pmaps.size() == nops);

    const int nQ = A0_sim_j.size(); assert(A0_sim_dj.size() == nQ);

    std::vector<std::vector<jackknifeDistributionD> > data_in_fit(nQ);
    std::ofstream data_in_fit_key("data_in_fit.key");

    std::cout << "Data in fit:\n";

    for(int q=0;q<nQ;q++){
      int eidx = 0;
      std::cout << "Q" << q+1 << std::endl;

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

	    std::cout << printCoord(c, pmap_descr) << " " << jack.value(d) << std::endl;
	    data_in_fit_key << "Q" << q+1 << " elem " << eidx << " " << pmap_descr.find(c.param_map)->second << " t=" << c.t << " tsep_K_snk=" << c.tsep_k_snk << std::endl;
	    data_in_fit[q].push_back(jack.value(d));
	    eidx++;
	  }
	}
      }
    }

    writeParamsStandard(data_in_fit,"data_in_fit.hdf5");
  }

  template<typename FitPolicies>
  static void runfit(std::vector<typename FitPolicies::FitParameterDistribution> &params,
		     const std::vector<SimFitCorrFuncJack> &A0_sim_j,
		     const std::vector<SimFitCorrFuncDJack> &A0_sim_dj,
		     const std::map<InnerParamMap const*, std::string> &pmap_descr,
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

  std::vector<jackknifeDistribution<Params> > fit(const ResampledData<jackknifeDistributionD> &data_j,
						  const ResampledData<doubleJackknifeA0StorageType> &data_dj,
						  const std::vector<PiPiOperator> &operators,
						  const int Lt, const int tmin_k_op, const int tmin_op_snk, bool correlated){
    
    assert(loaded_frzparams);
    
    const int nsample = mK.size();
    
    std::vector<std::string> dcp = { "AK", "mK" };
    for(int s=0;s<nstate;s++) dcp.push_back(stringize("E%d",s));
    for(int s=0;s<nstate;s++) dcp.push_back(stringize("M%d",s));

    //Parameter maps
    std::map< InnerParamMap const*, std::string> pmap_descr;  
    std::map<PiPiOperator, InnerParamMap> op_param_maps;
    
    std::unordered_map<std::string,size_t> *param_idx_map_ptr = new std::unordered_map<std::string,size_t>; //ensure it sticks around until end of execution
    auto &param_idx_map = *param_idx_map_ptr;

    //Data to be included
    std::vector<CorrFuncJackAllQ const* > incl_jacks;
    std::vector<CorrFuncDJackAllQ const* > incl_djacks;
    std::vector<InnerParamMap const *> incl_pmaps;

    //Set up the parameter maps
    int idx = 0;
#define DEFMAP(NM) param_idx_map[#NM] = idx++
    DEFMAP(AK);
    DEFMAP(mK);
#undef DEFMAP

    static const std::vector<PiPiOperator> all_ops = {PiPiOperator::PiPiGnd, PiPiOperator::PiPiExc, PiPiOperator::Sigma};
    static const std::vector<std::string> all_ops_pfmt = {"Apipi%d", "Apipi_exc_%d", "Asigma%d"};
    static const std::vector<std::string> all_ops_descr = {"K->pipi(111)", "K->pipi(311)", "K->sigma"};

    for(int opidx = 0; opidx < all_ops.size(); opidx++){
      auto op = all_ops[opidx];

      if(doOp(op, operators)){
	//Setup outer parameters and mapping to inner parameters
	auto & pmap = op_param_maps[op];
	for(int s=0;s<nstate;s++){
	  std::string pnm = stringize(all_ops_pfmt[opidx].c_str(),s);
	  param_idx_map[pnm] = idx++; //define outer parameter to index map
	  pmap[stringize("Asnk%d",s)] = pnm; //map inner parameter to outer parameter
	}
	for(auto it = dcp.begin(); it != dcp.end(); it++) pmap[*it] = *it; //those in dcp map one-to-one and are shared over all operators

	//Description of mapping
	pmap_descr[&pmap] = all_ops_descr[opidx];

	//Associate with the data
	assert(data_j.contains(op));
	incl_jacks.push_back(&data_j(op));
	incl_djacks.push_back(&data_dj(op));
	incl_pmaps.push_back(&pmap);
      }
    }
  
    //Add energies and matrix elements to end of list of outer parameters, also shared over all operators
    for(int s=0;s<nstate;s++)
      param_idx_map[stringize("E%d",s)] = idx++;
    for(int s=0;s<nstate;s++)
      param_idx_map[stringize("M%d",s)] = idx++;
    
    //Determine number of matrix elements from input data
    const int nQ = incl_jacks[0]->size();
    for(int i=0;i<incl_jacks.size();i++){ 
      assert(incl_jacks[i]->size() == nQ); 
      assert(incl_djacks[i]->size() == nQ); 
    }

    //Get the data in the fit range and put in sim fit data containers
    std::vector<SimFitCorrFuncJack> A0_sim_j(nQ);
    std::vector<SimFitCorrFuncDJack> A0_sim_dj(nQ);
    
    generateSimData(A0_sim_j,A0_sim_dj,
		    incl_jacks, incl_djacks, incl_pmaps,
		    tmin_k_op, tmin_op_snk,
		    pmap_descr);
  
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

    for(int opidx = 0; opidx < all_ops.size(); opidx++){
      auto op = all_ops[opidx];

      if(doOp(op, operators))
	for(int s=0; s<nstate; s++) 
	  freeze(freeze_params, freeze_vals, stringize(all_ops_pfmt[opidx].c_str(),s), coeffs[op][s], param_idx_map);
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

  jackknifeDistributionD weightedAvg(const std::vector<jackknifeDistributionD const*> &v){
    int nsample = v[0]->size();
    jackknifeDistributionD out(nsample,0.);
    double wsum = 0.;
    for(int i=0;i<v.size();i++){
      double w = v[i]->standardError();
      w = 1./w/w;
      wsum += w;
      out = out + w * (*v[i]);
    }
    out = out/wsum;
    return out;
  }


  std::vector<jackknifeDistribution<Params> > fit(const ResampledData<jackknifeDistributionD> &data_j,
						  const ResampledData<doubleJackknifeA0StorageType> &data_dj,
						  const std::vector<PiPiOperator> &operators,
						  const int Lt, const int tmin_k_op, const int tmin_op_snk, bool correlated){
    
    assert(this->loaded_frzparams);
    
    const int nsample = this->mK.size();
    
    std::vector<std::string> dcp;
    for(int s=0;s<nstate;s++) dcp.push_back(stringize("E%d",s));
    for(int s=0;s<nstate;s++) dcp.push_back(stringize("M%d",s));

    //Parameter maps
    std::map< InnerParamMap const*, std::string> pmap_descr;  
    std::map<PiPiOperator, InnerParamMap> op_param_maps;
    
    std::unordered_map<std::string,size_t> *param_idx_map_ptr = new std::unordered_map<std::string,size_t>; //ensure it sticks around until end of execution
    auto &param_idx_map = *param_idx_map_ptr;

    //Data to be included
    std::vector<CorrFuncJackAllQ const* > incl_jacks;
    std::vector<CorrFuncDJackAllQ const* > incl_djacks;
    std::vector<InnerParamMap const *> incl_pmaps;

    //Set up the parameter maps
    int idx = 0;

    static const std::vector<PiPiOperator> all_ops = {PiPiOperator::PiPiGnd, PiPiOperator::PiPiExc, PiPiOperator::Sigma};
    static const std::vector<std::string> all_ops_pfmt = {"Apipi%d", "Apipi_exc_%d", "Asigma%d"};
    static const std::vector<std::string> all_ops_descr = {"K->pipi(111)", "K->pipi(311)", "K->sigma"};

    for(int opidx = 0; opidx < all_ops.size(); opidx++){
      auto op = all_ops[opidx];

      if(doOp(op, operators)){
	//Setup outer parameters and mapping to inner parameters
	auto & pmap = op_param_maps[op];
	for(int s=0;s<nstate;s++){
	  std::string pnm = stringize(all_ops_pfmt[opidx].c_str(),s);
	  param_idx_map[pnm] = idx++; //define outer parameter to index map
	  pmap[stringize("Asnk%d",s)] = pnm; //map inner parameter to outer parameter
	}
	for(auto it = dcp.begin(); it != dcp.end(); it++) pmap[*it] = *it; //those in dcp map one-to-one and are shared over all operators

	//Description of mapping
	pmap_descr[&pmap] = all_ops_descr[opidx];

	//Associate with the data
	assert(data_j.contains(op));
	incl_jacks.push_back(&data_j(op));
	incl_djacks.push_back(&data_dj(op));
	incl_pmaps.push_back(&pmap);
      }
    }
  
    //Add energies and matrix elements to end of list of outer parameters, also shared over all operators
    for(int s=0;s<nstate;s++)
      param_idx_map[stringize("E%d",s)] = idx++;
    for(int s=0;s<nstate;s++)
      param_idx_map[stringize("M%d",s)] = idx++;
    
    //Determine number of matrix elements from input data
    const int nQ = incl_jacks[0]->size();
    for(int i=0;i<incl_jacks.size();i++){ 
      assert(incl_jacks[i]->size() == nQ); 
      assert(incl_djacks[i]->size() == nQ); 
    }

    //Get the data in the fit range and put in sim fit data containers
    std::vector<SimFitCorrFuncJack> A0_sim_j(nQ);
    std::vector<SimFitCorrFuncDJack> A0_sim_dj(nQ);
    
    generateSimData(A0_sim_j,A0_sim_dj,
		    incl_jacks, incl_djacks, incl_pmaps,
		    tmin_k_op, tmin_op_snk,
		    pmap_descr);
  
    std::vector<SimFitCorrFuncJack> A0_sim_j_wavg(nQ);
    std::vector<SimFitCorrFuncDJack> A0_sim_dj_wavg(nQ);

    for(int q=0;q<nQ;q++){
  
      //Transform the data to remove the kaon time dependence and amplitude
      std::map<std::pair<int,InnerParamMap const*>, std::vector<int> > to_avg;  //key is (top_snk, snk_op_pmap_descr)
      for(int d=0;d<A0_sim_j[q].size();d++){
	auto const &c = A0_sim_j[q].coord(d);
	int tK_op = (int)c.t;
	int top_snk = (int)c.tsep_k_snk - tK_op;

	to_avg[{top_snk, c.param_map}].push_back(d);

	jackknifeDistributionD fac = this->cK * exp(-this->mK * tK_op);
	
	for(int s=0;s<nsample;s++){
	  A0_sim_j[q].value(d).sample(s)  = A0_sim_j[q].value(d).sample(s) / fac.sample(s); 
	  A0_sim_dj[q].value(d).sample(s) = A0_sim_dj[q].value(d).sample(s) / fac.sample(s);  //note we use the same value for all inner samples - this is a minor effect
	}
      }
      
      //Perform the weighted average and insert into fit data containers
      for(auto it = to_avg.begin(); it != to_avg.end(); it++){
	SimFitCoordGen c; 
	c.t = it->first.first;
	c.param_map = it->first.second;

	std::cout << "Q" << q+1 << " weighted avg with descr " << pmap_descr.find(c.param_map)->second << " and top_snk = " << c.t << std::endl;

	auto const & idxv = it->second;

	std::vector<jackknifeDistributionD const*> towavg_j(idxv.size());
	
	std::cout << "Data included" << std::endl;
	//Jackknife
	for(int i=0;i<idxv.size();i++){ 
	  std::cout << printCoord(A0_sim_j[q].coord(idxv[i]), pmap_descr) << " " << A0_sim_j[q].value(idxv[i]) << std::endl;
	  towavg_j[i] = &A0_sim_j[q].value(idxv[i]);	
	}
	A0_sim_j_wavg[q].push_back(c, weightedAvg(towavg_j));
	
	//Double jackknife
	doubleJackknifeDistributionD wavg_dj(nsample);
	for(int s=0;s<nsample;s++){
	  for(int i=0;i<idxv.size();i++) towavg_j[i] = &A0_sim_dj[q].value(idxv[i]).sample(s);
	  wavg_dj.sample(s) = weightedAvg(towavg_j);
	}
	A0_sim_dj_wavg[q].push_back(c, wavg_dj);
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

    for(int opidx = 0; opidx < all_ops.size(); opidx++){
      auto op = all_ops[opidx];

      if(doOp(op, operators))
	for(int s=0; s<nstate; s++) 
	  freeze(freeze_params, freeze_vals, stringize(all_ops_pfmt[opidx].c_str(),s), coeffs[op][s], param_idx_map);
    }

    //Run the actual fit
    std::vector<jackknifeDistribution<Params> > params;    
    runfit<FitPolicies>(params, A0_sim_j_wavg, A0_sim_dj_wavg, pmap_descr, fitfunc, freeze_params, freeze_vals, guess, nsample, correlated);
 
    //plotErrorWeightedData2expFlat(ktopipi_A0_all_j, ktopipi_exc_A0_all_j, ktosigma_A0_all_j, params, Lt, tmin_k_op, tmin_op_snk, fitfunc);

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
