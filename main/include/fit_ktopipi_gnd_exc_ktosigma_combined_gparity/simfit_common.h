#ifndef _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_COMMON_H
#define _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_COMMON_H

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

struct printCoord{
  typedef std::unordered_map<std::string, std::string> InnerParamMap;
  const SimFitCoordGen &c;
  const std::map< InnerParamMap const*, std::string> &pmap_descr;
  printCoord(const SimFitCoordGen &c, const std::map< InnerParamMap const*, std::string> &pmap_descr): c(c), pmap_descr(pmap_descr){}
};
std::ostream & operator<<(std::ostream &os, const printCoord &p){ 
  auto it = p.pmap_descr.find(p.c.param_map); assert(it != p.pmap_descr.end());
  os << "(" << it->second << ", t=" << p.c.t << " tsep_K_snk=" << p.c.tsep_k_snk << " tsep_op_snk=" << p.c.tsep_k_snk - p.c.t << ")";
  return os;
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

  static inline std::string opDescrFile(PiPiOperator op){
    switch(op){
    case PiPiOperator::PiPiGnd:
      return "kpipi_111";
    case PiPiOperator::PiPiExc:
      return "kpipi_311";
    case PiPiOperator::Sigma:
      return "ksigma";
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

CPSFIT_END_NAMESPACE

#endif
