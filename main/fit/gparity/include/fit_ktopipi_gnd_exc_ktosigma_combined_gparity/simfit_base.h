#ifndef _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_BASE_H
#define _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_BASE_H

#include<config.h>
#include<utils/macros.h>
#include<pipi_common/fit/analyze_chisq.h>
#include "simfit_common.h"

CPSFIT_START_NAMESPACE

template<typename FitFunc,
	 template<typename, template<typename> class> class DistributionType>
void analyzeChisqFF(const correlationFunction<SimFitCoordGen,  DistributionType<double, basic_vector> > &corr_comb_j,
		    const DistributionType<typename FitFunc::ParameterType, basic_vector> &params, 
		    const FitFunc &fitfunc,
		    const NumericSquareMatrix<DistributionType<double, basic_vector> > &corr, 
		    const NumericVector<DistributionType<double, basic_vector> > &sigma,
		    const std::map<std::unordered_map<std::string, std::string> const*, std::string> &pmap_descr){
  struct PP{
    typedef std::map<std::unordered_map<std::string, std::string> const*, std::string> const* PtrType;
    inline static PtrType & descr(){ static PtrType p; return p; }

    inline static void print(std::ostream &os, const SimFitCoordGen &c){ os << printCoord(c, *descr()); }
    inline static std::string typeInfo(const SimFitCoordGen &c){ return descr()->find(c.param_map)->second; }   
  };
  PP::descr() = &pmap_descr;
  
  AnalyzeChisq<FitFunc,DistributionType,PP> chisq_analyze(corr_comb_j, fitfunc, params, corr, sigma);
  chisq_analyze.printChisqContribs(Correlation);
  chisq_analyze.examineProjectedDeviationContribsEvalNorm(Correlation);
  chisq_analyze.printChisqContribs(Covariance);
  chisq_analyze.examineProjectedDeviationContribsEvalNorm(Covariance);
}

template< template<typename, template<typename> class> class DistributionType > 
struct ResampledDataContainers{};

template< template<typename, template<typename> class> class DistributionType > 
struct SimFitDataContainers{};

//Note the param_idx in the FreezeParam entries should correspond to the state index starting at 0 for the ground state
struct FreezeMatrixElementOp{
#define FREEZE_MAT_ELEM_OP_MEMBERS			\
  (int, q)						\
  (std::vector<FreezeParam>, sources)

  GENERATE_MEMBERS(FREEZE_MAT_ELEM_OP_MEMBERS);
  
  FreezeMatrixElementOp(): sources(1), q(0){}
};
GENERATE_PARSER(FreezeMatrixElementOp, FREEZE_MAT_ELEM_OP_MEMBERS);
#undef FREEZE_MAT_ELEM_OP_MEMBERS

struct FreezeMatrixElements{
#define FREEZE_MATRIX_ELEMENTS_MEMBERS			\
  (std::vector<FreezeMatrixElementOp>, operators)

  GENERATE_MEMBERS(FREEZE_MATRIX_ELEMENTS_MEMBERS);
  
  FreezeMatrixElements(): operators(1){}
};
GENERATE_PARSER(FreezeMatrixElements, FREEZE_MATRIX_ELEMENTS_MEMBERS);
#undef FREEZE_MATRIX_ELEMENTS_MEMBERS

template< template<typename, template<typename> class> class DistributionType > 
class simultaneousFitBase: public simultaneousFitCommon{
public:
  COPY_COMMON_TYPEDEFS;

private:
  FreezeMatrixElements matelem_frz;
  std::vector<std::vector<DistributionType<double, basic_vector> > > matelem_frz_vals;

public:
  simultaneousFitBase(){ matelem_frz.operators.resize(0); }

  void loadFrozenMatrixElements(const FreezeMatrixElements &matelem_frz_, const int nsample){
    matelem_frz = matelem_frz_;
    matelem_frz_vals.resize(matelem_frz.operators.size());

    for(int qq=0;qq<matelem_frz.operators.size();qq++){
      matelem_frz_vals[qq].resize(matelem_frz.operators[qq].sources.size());

      for(int s=0;s<matelem_frz.operators[qq].sources.size();s++)
	readFrozenFitParam(matelem_frz_vals[qq][s], matelem_frz.operators[qq].sources[s], nsample);
    }
  }
  void loadFrozenMatrixElements(const std::string &matelem_frz_file, const int nsample){
    FreezeMatrixElements pargs;
    parse(pargs, matelem_frz_file);
    loadFrozenMatrixElements(pargs, nsample);
  }


  virtual void load2ptFitParams(const std::vector<PiPiOperator> &operators, const InputParamArgs &iargs, const int dist_size) = 0;

  virtual double evaluateFitFunc(const SimFitCoordGen &coord, const Params &params) const = 0;

  virtual double evaluateFitFunc(const amplitudeDataCoord &coord, PiPiOperator op, const Params &params) const = 0;

  virtual const subsetMapDescr & getParameterMapDescr() const = 0;

  virtual void fit(std::vector<DistributionType<Params, basic_vector> > &params,
		   std::vector<DistributionType<double, basic_vector> > &chisq,
		   SimFitDataContainers<DistributionType> &simfit_data,
		   const std::vector<PiPiOperator> &operators,
		   const int Lt, bool correlated, 
		   const CovarianceMatrix covariance_matrix,
		   bool write_output = true) const = 0;

  virtual void generateSimFitData(SimFitDataContainers<DistributionType> &simfit_data,
				  const ResampledDataContainers<DistributionType> &fit_data,
				  const std::vector<PiPiOperator> &operators,
				  const int Lt, const int tmin_k_op, const int tmin_op_snk,
				  const CovarianceMatrix covariance_matrix) const = 0;

  //Combines the above two, with possible extra output (eg plots) that require the raw data
  virtual void fit(std::vector<DistributionType<Params, basic_vector> > &params,
		   std::vector<DistributionType<double, basic_vector> > &chisq,		   
		   const ResampledDataContainers<DistributionType> &fit_data,
		   const std::vector<PiPiOperator> &operators,
		   const int Lt, const int tmin_k_op, const int tmin_op_snk, 
		   bool correlated, const CovarianceMatrix covariance_matrix, bool write_output = true) const = 0;


  template<typename FitFunc>
  void runfit(std::vector<DistributionType<Params, basic_vector> > &params,
		     std::vector<DistributionType<double, basic_vector> > &chisq,
		     const SimFitDataContainers<DistributionType> &fit_data,
		     const subsetMapDescr &pmap_descr,
		     const FitFunc &fitfunc,
		     const std::vector<int> &freeze_params,
		     const DistributionType<Params, basic_vector> &freeze_vals,
		     const Params &guess,
		     const bool correlated, const CovarianceMatrix covariance_matrix, bool write_output) const{
    
    if(write_output) simultaneousFitCommon::printWriteFitData(fit_data.getFitDataAllQ(),pmap_descr);

    typedef DistributionType<double, basic_vector> DistributionTypeD;

    typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, correlatedFitPolicy, 
				      MarquardtLevenbergMinimizerPolicy, DistributionTypeD>::type FitPolicies;

    const int nQ = fit_data.getNq();
    params.resize(nQ);
    chisq.resize(nQ);

    auto init = fit_data.getDistributionInitializer();
    DistributionTypeD base(init);

    std::vector<DistributionTypeD>  chisq_per_dof(nQ, base), pvalue(nQ, base);

    for(int q=0;q<nQ;q++){
      chisq[q] = base;
      params[q] = DistributionType<Params, basic_vector>(init, guess);

      std::cout << "Performing " << q+1 << " fit" << std::endl;
      fitter<FitPolicies> fit;
      fit.importFitFunc(fitfunc);

      std::vector<int> freeze_params_q = freeze_params;
      DistributionType<Params, basic_vector> freeze_vals_q = freeze_vals;

      for(int qq=0;qq<matelem_frz.operators.size();qq++){ 
	if(matelem_frz.operators[qq].q == q){
	  for(int src=0;src<matelem_frz.operators[qq].sources.size();src++){
	    std::string Mnm = stringize("M%d",matelem_frz.operators[qq].sources[src].param_idx);
	    int Midx = freeze_vals_q.sample(0).index(Mnm);

	    for(int s=0;s<iterate<DistributionType<Params, basic_vector>>::size(freeze_vals_q);s++)
	      iterate<DistributionType<Params, basic_vector>>::at(s, freeze_vals_q)(Midx) = iterate<DistributionType<double, basic_vector>>::at(s, matelem_frz_vals[qq][src]);

	    freeze_params_q.push_back(Midx);
	  }
	}
      }

      fit.freeze(freeze_params_q, freeze_vals_q);
      
      importCostFunctionParameters<correlatedFitPolicy, FitPolicies> import;
      fit_data.generateCovarianceMatrix(import, fit, covariance_matrix, q);

      if(!correlated) import.setUncorrelated();
          
      std::cout << "Data in fit" << std::endl;
      const auto &data_q = fit_data.getFitData(q);

      for(int i=0;i<data_q.size();i++)
	std::cout << printCoord(data_q.coord(i), pmap_descr) << " " << data_q.value(i) << std::endl;

      std::cout << "Starting fit" << std::endl;

      int ndof;
      fit.fit(params[q], chisq[q], chisq_per_dof[q], ndof, data_q);

      pvalue[q] = DistributionTypeD(init);
      for(int i=0; i<iterate<DistributionTypeD>::size(pvalue[q]); i++) 
	iterate<DistributionTypeD>::at(i, pvalue[q]) = chiSquareDistribution::pvalue(ndof, iterate<DistributionTypeD>::at(i, chisq[q]) );

      std::cout << "Q" << q+1 << std::endl;
      std::cout << "Params:\n" << params[q] << std::endl;
      std::cout << "Chisq: " << chisq[q] << std::endl;
      std::cout << "Chisq/dof: " << chisq_per_dof[q] << std::endl;
      std::cout << "p-value(chi^2): " << pvalue[q] << std::endl;

      if(write_output) analyzeChisqFF(fit_data.getFitData(q), params[q], fitfunc, import.corr, import.sigma, pmap_descr);
    }  
    std::stringstream key;

    for(int q=0;q<nQ;q++){
      int nparam = params[q].sample(0).size();
      std::cout << "Q" << q+1 << std::endl;
      std::cout << "Params:" << std::endl;
      DistributionTypeD tmp;
      for(int i=0;i<nparam;i++){
	standardIOhelper<DistributionTypeD, DistributionType<Params, basic_vector> >::extractStructEntry(tmp, params[q], i);
	std::string label = params[q].sample(0).tag(i);
	std::cout << i << " " << label << " " << tmp << std::endl;
	key << q << " " << i << " " << label << std::endl;
      }
      std::cout << "Chisq: " << chisq[q] << std::endl;
      std::cout << "Chisq/dof: " << chisq_per_dof[q] << std::endl;
      std::cout << "p-value: " << pvalue[q] << std::endl;
    }

    if(write_output){
      std::ofstream kk("params.key");
      kk << key.str();
      writeParamsStandard(params, "params.hdf5");
      writeParamsStandard(chisq, "chisq.hdf5");
      writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");
      writeParamsStandard(pvalue, "pvalue.hdf5");
    }
  }


  virtual ~simultaneousFitBase(){}
};

CPSFIT_END_NAMESPACE

#endif
