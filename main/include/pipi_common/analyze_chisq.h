#ifndef PIPI_ANALYZE_CHISQ_H
#define PIPI_ANALYZE_CHISQ_H

#include<common.h>
#include<tensors.h>

using namespace CPSfit;

template<typename CoordType>
struct CoordPrintPolicyBasic{
  inline static void print(std::ostream &os, const CoordType &c){ os << c; }
  inline static std::string typeInfo(const CoordType &c){ return "all"; } //allows data contributions to chi^2 to be sorted by a type
};

enum WhichMatrix{ Covariance, Correlation, CorrelationLedoitWolf };

template<typename FitFunc, typename CoordPrintPolicy = CoordPrintPolicyBasic<typename FitFunc::GeneralizedCoordinate>  >
class AnalyzeChisq{
  typedef typename FitFunc::GeneralizedCoordinate GeneralizedCoordinate;
  typedef typename FitFunc::ParameterType ParameterType;
  
  static NumericSquareMatrix<double> computeFrozenCovMat(const correlationFunction<GeneralizedCoordinate, jackknifeDistributionD> &data_inrange){
    NumericSquareMatrix<double> cov(data_inrange.size());
    for(int i=0;i<cov.size();i++)
      for(int j=0;j<cov.size();j++)
	cov(i,j) = jackknifeDistributionD::covariance(data_inrange.value(i),data_inrange.value(j));
    return cov;
  }
  
  static void computeFrozenCorrSigma(NumericSquareMatrix<double> &corr, NumericVector<double> &sigma, 
				     const NumericSquareMatrix<double> &cov){
    int N = cov.size();
    sigma.resize(N);
    corr.resize(N);
    for(int i=0;i<N;i++){
      sigma(i) = sqrt(cov(i,i));
      for(int j=0;j<N;j++)
	corr(i,j) = cov(i,j)/sqrt(cov(i,i))/sqrt(cov(j,j));
    }
  }

  static void computeLedoitWolfCorr(NumericSquareMatrix<double> &corrlw, const NumericSquareMatrix<double> &corr,
				    double ledoit_wolf_lambda){    
    int N = corr.size();
    corrlw = corr;
    if(ledoit_wolf_lambda != 0.){
      for(int i=0;i<N;i++){
	for(int j=0;j<N;j++)
	  if(i==j) corrlw(i,j) = ledoit_wolf_lambda + (1.-ledoit_wolf_lambda)*corrlw(i,j);
	  else corrlw(i,j) = (1.-ledoit_wolf_lambda)*corrlw(i,j);
      }
    }
  }


  static void computeEvalsEvecs(std::vector<double> &evals, std::vector<NumericVector<double> > &evecs,
			 const NumericSquareMatrix<double> &M){    
    symmetricMatrixEigensolve(evecs,evals,M);
  }
  
  static NumericVector<double> projectOntoEvecs(const NumericVector<double> &v, const std::vector<NumericVector<double> > &evecs){
    int Nevec = evecs.size();
    int Ndata = v.size();
    NumericVector<double> out(Nevec);
    for(int i=0;i<Nevec;i++){
      out(i) = dot(evecs[i], v);
    }
    return out;
  }

  //Evaluate the fit func across the data range using the input parameters
  //Need data only for the generalized coordinates
  static NumericVector<double> evalFitFunc(const FitFunc &fitfunc, const ParameterType &params, 
				    const correlationFunction<GeneralizedCoordinate, jackknifeDistributionD> &data_inrange){
    NumericVector<double> out(data_inrange.size());
    for(int i=0;i<data_inrange.size();i++)
      out(i) = fitfunc.value(data_inrange.coord(i), params);
    return out;
  }
  //Same as above but use central value of fit param jackknife
  static inline NumericVector<double> evalFitFunc(const FitFunc &fitfunc, const jackknifeDistribution<ParameterType> &params, 
					   const correlationFunction<GeneralizedCoordinate, jackknifeDistributionD> &data_inrange){
    return evalFitFunc(fitfunc, params.best(), data_inrange);
  }
   
  const correlationFunction<GeneralizedCoordinate, jackknifeDistributionD> &data_inrange;
  const FitFunc &fitfunc;
  const jackknifeDistribution<ParameterType> &params;

  NumericSquareMatrix<double> cov;
  std::vector<double> cov_evals;
  std::vector<NumericVector<double> > cov_evecs;

  NumericSquareMatrix<double> corr;
  NumericVector<double> sigma;
  std::vector<double> corr_evals;
  std::vector<NumericVector<double> > corr_evecs;

  //Ledoit-Wolf correlation matrix (same as regular correlation matrix for lambda=0)
  NumericSquareMatrix<double> corrlw;
  std::vector<double> corrlw_evals;
  std::vector<NumericVector<double> > corrlw_evecs;

  const std::vector<double> &getEvals(const WhichMatrix m) const{
    switch(m){
    case Covariance:
      return cov_evals;
    case Correlation:
      return corr_evals;
    case CorrelationLedoitWolf:
      return corrlw_evals;
    }
  }
  const std::vector<NumericVector<double> > &getEvecs(const WhichMatrix m) const{
    switch(m){
    case Covariance:
      return cov_evecs;
    case Correlation:
      return corr_evecs;
    case CorrelationLedoitWolf:
      return corrlw_evecs;
    }
  }     

  NumericVector<double> fitval_cen;
  NumericVector<double> fitval_cen_nrm;
  NumericVector<double> dataval_cen;
  NumericVector<double> dataval_cen_nrm;
  NumericVector<double> delta;
  NumericVector<double> delta_nrm;

  const NumericVector<double> &getFitVals(const WhichMatrix m) const{
    switch(m){
    case Covariance:
      return fitval_cen;
    case Correlation:
    case CorrelationLedoitWolf:
      return fitval_cen_nrm;
    }
  }     
  const NumericVector<double> &getDataVals(const WhichMatrix m) const{
    switch(m){
    case Covariance:
      return dataval_cen;
    case Correlation:
    case CorrelationLedoitWolf:
      return dataval_cen_nrm;
    }
  }  
  const NumericVector<double> &getDeltaVals(const WhichMatrix m) const{
    switch(m){
    case Covariance:
       return delta;
    case Correlation:
    case CorrelationLedoitWolf:
       return delta_nrm;
    }
  }  

  inline std::string nm(const WhichMatrix m) const{ 
    switch(m){
    case Covariance:
      return "covariance matrix";
    case Correlation:
      return "correlation matrix";
    case CorrelationLedoitWolf:
      return "Ledoit-Wolf correlation matrix";
    }
  }  

public:

  AnalyzeChisq(const correlationFunction<GeneralizedCoordinate, jackknifeDistributionD> &data_inrange,
	       const FitFunc &fitfunc,
	       const jackknifeDistribution<ParameterType> &params,
	       double ledoit_wolf_lambda = 0.): data_inrange(data_inrange), fitfunc(fitfunc), params(params){
    //Covariance matrix
    cov = computeFrozenCovMat(data_inrange);
    computeEvalsEvecs(cov_evals, cov_evecs, cov);
    
    //Correlation matrix
    computeFrozenCorrSigma(corr, sigma, cov);
    computeEvalsEvecs(corr_evals, corr_evecs, corr);

    //Ledoit-Wolf correlation matrix
    computeLedoitWolfCorr(corrlw, corr, ledoit_wolf_lambda);
    computeEvalsEvecs(corrlw_evals, corrlw_evecs, corrlw);

    //Fit and data
    fitval_cen = evalFitFunc(fitfunc, params, data_inrange);
    fitval_cen_nrm = NumericVector<double>(fitval_cen.size(), [&](const int i){ return fitval_cen(i)/sigma(i); });

    dataval_cen = NumericVector<double>(data_inrange.size(), [&](const int i){ return data_inrange.value(i).best(); });
    dataval_cen_nrm = NumericVector<double>(dataval_cen.size(), [&](const int i){ return dataval_cen(i)/sigma(i); });

    delta = NumericVector<double>(data_inrange.size(), [&](const int i){ return (dataval_cen(i) - fitval_cen(i)); });  
    delta_nrm = NumericVector<double>(data_inrange.size(), [&](const int i){ return (dataval_cen(i) - fitval_cen(i))/sigma(i); });
  }
  
  void printChisqContribs(const WhichMatrix mat, std::ostream &os = std::cout) const{
    os << "Examining contributions to chi^2 from eigenmodes of the " << nm(mat) << ":\n";

    const std::vector<double> &evals = getEvals(mat);
    const std::vector<NumericVector<double> > &evecs = getEvecs(mat);
    const NumericVector<double> &delta = getDeltaVals(mat);
    NumericVector<double> dproj = projectOntoEvecs(delta,evecs);
    
    double chisq_tot = 0.;
    for(int i=0;i<evals.size();i++){
      double chisq_contrib = dproj(i)*dproj(i)/evals[i];
      os << i << " (eval " << evals[i] << ") " << chisq_contrib << std::endl;
      chisq_tot += chisq_contrib;
    }
    os << "Total chi^2: " << chisq_tot << std::endl;
  }

  void examineEigenvectors(const WhichMatrix mat, std::ostream &os = std::cout) const{
    os << "Examining eigenvectors of the " << nm(mat) << ":\n";
    
    const std::vector<NumericVector<double> > &evecs = getEvecs(mat);
    for(int i=0;i<evecs.size();i++){
      os << "Mode " << i << ":\n";
      for(int tt=0;tt<fitval_cen.size();tt++){
	CoordPrintPolicy::print(os,data_inrange.coord(tt));
	os << " " << evecs[i](tt) << "\n";
      }
    }
  }

  void examineProjectedFitFuncContribs(const WhichMatrix mat, std::ostream &os = std::cout) const{
    os << "Examining data coordinate contributions to " << nm(mat) << " mode-projected fit function:\n";
    const std::vector<NumericVector<double> > &evecs = getEvecs(mat);
    const NumericVector<double> &fitval = getFitVals(mat);

    for(int i=0;i<evecs.size();i++){
      double tot = 0;
      os << "Mode " << i << ":\n";
      for(int tt=0;tt<fitval.size();tt++){
	double contrib = evecs[i](tt) * fitval(tt);
	CoordPrintPolicy::print(os,data_inrange.coord(tt));
	os << " " << contrib << "\n";
	tot += contrib;
      }
      os << "Projected mode total: " << tot << std::endl;
    }
  }

  void examineProjectedDeviationContribs(const WhichMatrix mat, std::ostream &os = std::cout) const{
    os << "Examining data coordinate contributions of " << nm(mat) << " to mode-projected deviations:\n";
    const std::vector<NumericVector<double> > &evecs = getEvecs(mat);
    const std::vector<double> &evals = getEvals(mat);
    const NumericVector<double> &fitval = getFitVals(mat);
    const NumericVector<double> &dataval = getDataVals(mat);

    for(int i=0;i<evecs.size();i++){
      double tot = 0;
      os << "Mode " << i << ":\n";
      for(int tt=0;tt<fitval.size();tt++){
	double data_contrib = evecs[i](tt) * dataval(tt);
	double fit_contrib = evecs[i](tt) * fitval(tt);
	double contrib = data_contrib - fit_contrib;

	CoordPrintPolicy::print(os, data_inrange.coord(tt));
	os << " " << data_contrib << " - " << fit_contrib << " = " << contrib << "\n";
	tot += contrib;
      }
      os << "Projected mode total: " << tot << std::endl;
      os << "Contrib to chi^2: (" << tot << ")^2/" << evals[i] << " = " << tot*tot/evals[i] << std::endl; 
    }
  }

  void examineProjectedDeviationContribsEvalNorm(const WhichMatrix mat, std::ostream &os = std::cout) const{
    os << "Examining data coordinate contributions of " << nm(mat) << " to mode-projected deviations normalized by sqrt(eval):\n";
    const std::vector<NumericVector<double> > &evecs = getEvecs(mat);
    const std::vector<double> &evals = getEvals(mat);
    const NumericVector<double> &fitval = getFitVals(mat);
    const NumericVector<double> &dataval = getDataVals(mat);

    std::vector< std::vector<double> > fraction_worse( fitval.size(),  std::vector<double>(evecs.size()) );

    for(int i=0;i<evecs.size();i++){
      std::vector<std::pair<int,double> > contribs(fitval.size());
      std::map<std::string, double> type_contribs;

      double tot = 0;
      os << "Mode " << i << ":\n";
      for(int tt=0;tt<fitval.size();tt++){
	double evec_nrm_i_t = evecs[i](tt)/sqrt(evals[i]);

	double data_contrib = evec_nrm_i_t *  dataval(tt);
	double fit_contrib = evec_nrm_i_t * fitval(tt);
	double contrib = data_contrib - fit_contrib;

	CoordPrintPolicy::print(os, data_inrange.coord(tt));
	os << " (" << evec_nrm_i_t << ") : " << data_contrib << " - " << fit_contrib << " = " << contrib << "\n";
	tot += contrib;

	contribs[tt] = {tt, contrib};
	
	//Sum up contributions for each data type
	std::string type = CoordPrintPolicy::typeInfo(data_inrange.coord(tt));
	auto it = type_contribs.find(type);
	if(it != type_contribs.end()) it->second += contrib;
	else type_contribs[type] = contrib;
      }
      os << "Projected eval-normalized mode total: " << tot << std::endl;
      os << "Contrib to chi^2: (" << tot << ")^2" << " = " << tot*tot << std::endl; 

      std::sort( contribs.begin(), contribs.end(), [&](const std::pair<int,double> &a, const std::pair<int,double> &b){ return fabs(a.second) < fabs(b.second); });
      os << "Contributions sorted by absolute value:\n";
      for(int tt=0;tt<fitval.size();tt++){
	int didx =  contribs[tt].first;
	CoordPrintPolicy::print(os, data_inrange.coord(didx));
	os << " " << contribs[tt].second << std::endl;

	fraction_worse[didx][i] = double(fitval.size()-tt-1)/fitval.size();	
      }

      os << "Contributions per data type (total " << tot << "):\n";
      for(auto it=type_contribs.begin(); it!=type_contribs.end(); it++){
	os << it->first << " " << it->second << std::endl;
      }
    }
    
    os << "For each data point and mode, the fraction of data points with a larger contribution than that of this point (sorted):\n";    
    
    for(int tt=0;tt<fitval.size();tt++){
      std::sort(fraction_worse[tt].begin(), fraction_worse[tt].end());
      CoordPrintPolicy::print(os, data_inrange.coord(tt));
      os << " :";
      
      double avg = 0.;
      for(int i=0;i<evecs.size();i++){
	avg += fraction_worse[tt][i];
	os << " " << fraction_worse[tt][i];
      }
      os << "  Average " << avg/evecs.size() << std::endl;
    }
  }


};



#endif
