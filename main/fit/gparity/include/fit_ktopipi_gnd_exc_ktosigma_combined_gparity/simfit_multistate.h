#ifndef _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_MULTISTATE_H
#define _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_MULTISTATE_H

#include<config.h>
#include<utils/macros.h>

#include "simfit_base.h"

SARLAC_START_NAMESPACE




template< template<typename, template<typename> class> class DistributionType > 
class simultaneousFitMultiState: public simultaneousFitBase<DistributionType>{
public:
  COPY_COMMON_TYPEDEFS;    
  typedef FitSimGenMultiState FitFunc;    

protected:
  int nstate;

  //Values for freezing
  bool loaded_frzparams;
  
  typedef DistributionType<double, basic_vector> DistributionTypeD;
  DistributionTypeD mK;
  DistributionTypeD cK;
  std::vector<DistributionTypeD> E;
  std::map<PiPiOperator, std::vector<DistributionTypeD> > coeffs; //for each operator, the frozen amplitudes (0=gnd state, 1=exc state, ...)

  paramIdxMap param_idx_map;  //make sure you don't delete the fitter before the end of execution or it will invalidate the parameter maps used by the output fit parameters!
  operatorSubsetMap op_param_maps;
  subsetMapDescr pmap_descr;

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

private:
  FitFunc *fitfunc;

public:

  simultaneousFitMultiState(const int nstate): loaded_frzparams(false), nstate(nstate), fitfunc(NULL){ //doesn't construct fit func
    assert(nstate <= 3);
  }

  simultaneousFitMultiState(const int nstate, const std::vector<PiPiOperator> &operators): loaded_frzparams(false), nstate(nstate){
    assert(nstate <= 3);
    constructParameterMaps(this->param_idx_map, this->op_param_maps, this->pmap_descr, operators);
    this->fitfunc = new FitFunc(this->param_idx_map.size(), nstate);
  }
  virtual ~simultaneousFitMultiState(){ 
    if(fitfunc) delete fitfunc; 
  }

  double evaluateFitFunc(const SimFitCoordGen &coord, const Params &params) const{
    return fitfunc->value(coord, params);
  }
  double evaluateFitFunc(const amplitudeDataCoord &coord, PiPiOperator op, const Params &params) const{
    InnerParamMap const* op_map = &op_param_maps.find(op)->second;
    SimFitCoordGen sim_coord(coord.t, coord.tsep_k_pi, op_map);
    return evaluateFitFunc(sim_coord, params);
  }

  const subsetMapDescr & getParameterMapDescr() const{ return pmap_descr; }

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
  
  void generateSimFitData(SimFitDataContainers<DistributionType> &simfit_data,
			  const ResampledDataContainers<DistributionType> &fit_data,
			  const std::vector<PiPiOperator> &operators,
			  const int Lt, const int tmin_k_op, const int tmin_op_snk,
			  const CovarianceMatrix covariance_matrix) const{
    simfit_data.generateSimData(fit_data, operators, tmin_k_op, tmin_op_snk, op_param_maps, pmap_descr, covariance_matrix);
  }

  void fit(std::vector<DistributionType<Params, basic_vector> > &params,
	   std::vector<DistributionType<double, basic_vector> > &chisq,
	   const ResampledDataContainers<DistributionType> &fit_data,
	   const std::vector<PiPiOperator> &operators,
	   const int Lt, const int tmin_k_op, const int tmin_op_snk, bool correlated, 
	   const CovarianceMatrix covariance_matrix,
	   bool write_output = true) const{
    SimFitDataContainers<DistributionType> simfit_data;
    generateSimFitData(simfit_data, fit_data, operators, Lt, tmin_k_op, tmin_op_snk, covariance_matrix);
    
    fit(params, chisq, simfit_data, operators, Lt, correlated, covariance_matrix);

    if(write_output) plotErrorWeightedDataNexpFlat(fit_data.getFitData(), operators, params, this->mK, this->cK, this->nstate, Lt, tmin_k_op, tmin_op_snk);
  }

  void fit(std::vector<DistributionType<Params, basic_vector> > &params,
	   std::vector<DistributionType<double, basic_vector> > &chisq,
	   SimFitDataContainers<DistributionType> &simfit_data,
	   const std::vector<PiPiOperator> &operators,
	   const int Lt, bool correlated, 
	   const CovarianceMatrix covariance_matrix,
	   bool write_output = true) const{
    
    assert(loaded_frzparams);
   
    //Setup guess
    Params guess(&param_idx_map);
    for(int i=0;i<guess.size();i++) guess(i) = 1.;
    for(int s=0;s<nstate;s++)
      guess(stringize("M%d",s)) = 0.5;

    auto init = simfit_data.getFitData(0).value(0).getInitializer();

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
	this->freeze(freeze_params, freeze_vals, stringize(opAmplitudeParamFmt(op).c_str(),s), coeffs.find(op)->second[s], param_idx_map);
    }

    //Run the actual fit
    this->runfit(params, chisq, simfit_data, pmap_descr, *this->fitfunc, freeze_params, freeze_vals, guess, correlated, covariance_matrix, write_output);
  }
};




SARLAC_END_NAMESPACE

#endif
