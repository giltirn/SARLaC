#ifndef _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_MULTISTATE_WAVG_H
#define _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_MULTISTATE_WAVG_H

#include<config.h>
#include<utils/macros.h>

#include "simfit_multistate.h"

CPSFIT_START_NAMESPACE

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
 
    plotErrorWeightedDataNexpFlat(fit_data.getFitData(), operators, params, this->mK, this->cK, this->nstate, Lt, tmin_k_op, tmin_op_snk);

    return params;
  }
};

CPSFIT_END_NAMESPACE

#endif
