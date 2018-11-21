#ifndef PIPI_ANALYZE_CHISQ_H
#define PIPI_ANALYZE_CHISQ_H

#include<common.h>
#include<tensors.h>

using namespace CPSfit;

template<typename CoordType>
struct CoordPrintPolicyBasic{
  inline static void print(std::ostream &os, const CoordType &c){ os << c; }
};

enum WhichMatrix{ Covariance, Correlation };

template<typename FitFunc, typename CoordPrintPolicy = CoordPrintPolicyBasic<typename FitFunc::GeneralizedCoordinate>  >
class AnalyzeChisq{
  typedef typename FitFunc::GeneralizedCoordinate GeneralizedCoordinate;
  typedef typename FitFunc::ParameterType ParameterType;
  
  NumericSquareMatrix<double> computeFrozenCovMat(const correlationFunction<GeneralizedCoordinate, jackknifeDistributionD> &data_inrange) const{
    NumericSquareMatrix<double> cov(data_inrange.size());
    for(int i=0;i<cov.size();i++)
      for(int j=0;j<cov.size();j++)
	cov(i,j) = jackknifeDistributionD::covariance(data_inrange.value(i),data_inrange.value(j));
    return cov;
  }
  
  void computeFrozenCorrSigma(NumericSquareMatrix<double> &corr, NumericVector<double> &sigma, 
			      const NumericSquareMatrix<double> &cov) const{
    int N = cov.size();
    sigma.resize(N);
    corr.resize(N);
    for(int i=0;i<N;i++){
      sigma(i) = sqrt(cov(i,i));
      for(int j=0;j<N;j++)
	corr(i,j) = cov(i,j)/sqrt(cov(i,i))/sqrt(cov(j,j));
    }
  }

  void computeEvalsEvecs(std::vector<double> &evals, std::vector<NumericVector<double> > &evecs,
			 const NumericSquareMatrix<double> &M) const{    
    symmetricMatrixEigensolve(evecs,evals,M);
  }
  
  NumericVector<double> projectOntoEvecs(const NumericVector<double> &v, const std::vector<NumericVector<double> > &evecs) const{
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
  NumericVector<double> evalFitFunc(const FitFunc &fitfunc, const ParameterType &params, 
				    const correlationFunction<GeneralizedCoordinate, jackknifeDistributionD> &data_inrange) const {
    NumericVector<double> out(data_inrange.size());
    for(int i=0;i<data_inrange.size();i++)
      out(i) = fitfunc.value(data_inrange.coord(i), params);
    return out;
  }
  //Same as above but use central value of fit param jackknife
  inline NumericVector<double> evalFitFunc(const FitFunc &fitfunc, const jackknifeDistribution<ParameterType> &params, 
					   const correlationFunction<GeneralizedCoordinate, jackknifeDistributionD> &data_inrange) const{
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

  NumericVector<double> fitval_cen;
  NumericVector<double> dataval_cen;
  NumericVector<double> delta;
  NumericVector<double> delta_covproj; //projected onto basis of cov mat evecs
  NumericVector<double> delta_nrm;
  NumericVector<double> delta_nrm_corrproj; //on basis of corr mat evecs

  inline std::string nm(const WhichMatrix m)const { return m == Correlation ? "correlation matrix" : "covariance matrix"; }
public:

  AnalyzeChisq(const correlationFunction<GeneralizedCoordinate, jackknifeDistributionD> &data_inrange,
	       const FitFunc &fitfunc,
	       const jackknifeDistribution<ParameterType> &params): data_inrange(data_inrange), fitfunc(fitfunc), params(params){
    cov = computeFrozenCovMat(data_inrange);
    computeEvalsEvecs(cov_evals, cov_evecs, cov);
    computeFrozenCorrSigma(corr, sigma, cov);
    computeEvalsEvecs(corr_evals, corr_evecs, corr);

    fitval_cen = evalFitFunc(fitfunc, params, data_inrange);
    dataval_cen = NumericVector<double>(data_inrange.size(), [&](const int i){ return data_inrange.value(i).best(); });

    delta = NumericVector<double>(data_inrange.size(), [&](const int i){ return (dataval_cen(i) - fitval_cen(i)); });
    delta_covproj = projectOntoEvecs(delta, cov_evecs);
    
    delta_nrm = NumericVector<double>(data_inrange.size(), [&](const int i){ return (dataval_cen(i) - fitval_cen(i))/sigma(i); });
    delta_nrm_corrproj = projectOntoEvecs(delta_nrm, corr_evecs);
  }
  
  void printChisqContribs(const WhichMatrix mat, std::ostream &os = std::cout) const{
    os << "Examining contributions to chi^2 from eigenmodes of the " << nm(mat) << ":\n";

    const std::vector<double> &evals = mat == Correlation ? corr_evals : cov_evals;
    const NumericVector<double> &dproj = mat == Correlation ? delta_nrm_corrproj : delta_covproj;

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
    
    const std::vector<NumericVector<double> > &evecs = mat == Correlation ? corr_evecs : cov_evecs;
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
    const std::vector<NumericVector<double> > &evecs = mat == Correlation ? corr_evecs : cov_evecs;

    for(int i=0;i<evecs.size();i++){
      double tot = 0;
      os << "Mode " << i << ":\n";
      for(int tt=0;tt<fitval_cen.size();tt++){
	double contrib = evecs[i](tt) * fitval_cen(tt);
	if(mat == Correlation) contrib /= sigma(tt);
	CoordPrintPolicy::print(os,data_inrange.coord(tt));
	os << " " << contrib << "\n";
	tot += contrib;
      }
      os << "Projected mode total: " << tot << std::endl;
    }
  }

  void examineProjectedDeviationContribs(const WhichMatrix mat, std::ostream &os = std::cout) const{
    os << "Examining data coordinate contributions of " << nm(mat) << " to mode-projected deviations:\n";
    const std::vector<NumericVector<double> > &evecs = mat == Correlation ? corr_evecs : cov_evecs;
    const std::vector<double> &evals = mat == Correlation ? corr_evals : cov_evals;

    for(int i=0;i<evecs.size();i++){
      double tot = 0;
      os << "Mode " << i << ":\n";
      for(int tt=0;tt<fitval_cen.size();tt++){
	double data_contrib = evecs[i](tt) * dataval_cen(tt);
	double fit_contrib = evecs[i](tt) * fitval_cen(tt);
	if(mat == Correlation){
	  data_contrib /= sigma(tt);
	  fit_contrib /= sigma(tt);
	}
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
    const std::vector<NumericVector<double> > &evecs = mat == Correlation ? corr_evecs : cov_evecs;
    const std::vector<double> &evals = mat == Correlation ? corr_evals : cov_evals;

    for(int i=0;i<evecs.size();i++){
      double tot = 0;
      os << "Mode " << i << ":\n";
      for(int tt=0;tt<fitval_cen.size();tt++){
	double evec_nrm_i_t = evecs[i](tt)/sqrt(evals[i]);

	double data_contrib = evec_nrm_i_t *  dataval_cen(tt);
	double fit_contrib = evec_nrm_i_t * fitval_cen(tt);
	
	if(mat == Correlation){
	  data_contrib /= sigma(tt);
	  fit_contrib /= sigma(tt);
	}
	double contrib = data_contrib - fit_contrib;

	CoordPrintPolicy::print(os, data_inrange.coord(tt));
	os << " (" << evec_nrm_i_t << ") : " << data_contrib << " - " << fit_contrib << " = " << contrib << "\n";
	tot += contrib;
      }
      os << "Projected eval-normalized mode total: " << tot << std::endl;
      os << "Contrib to chi^2: (" << tot << ")^2" << " = " << tot*tot << std::endl; 
    }
  }

  // //Compute mean and std.dev of data-point projected deviations onto normalized evecs to see if certain timeslices are systematically giving large chi^2 contributions
  // void examineProjectedDeviationDistributionEvalNorm(std::ostream &os = std::cout) const{
  //   os << "Examining distribution of absolute data-point deviations over basis of normalized eigenvectors:\n";
    
  //   std::vector<rawDataDistributionD> dev_dist(fitval_cen.size(), rawDataDistributionD(corr_evals.size()));
    
  //   for(int i=0;i<corr_evals.size();i++){
  //     for(int tt=0;tt<fitval_cen.size();tt++){
  // 	double evec_nrm_i_t = corr_evecs[i](tt)/sqrt(corr_evals[i]);
  // 	double contrib = evec_nrm_i_t * delta_nrm(tt);
	
  // 	dev_dist[tt].sample(i) = fabs(contrib);
  //     }
  //   }
  //   for(int tt=0;tt<fitval_cen.size();tt++){
  //     CoordPrintPolicy::print(os, data_inrange.coord(tt));
  //     os << " " << dev_dist[tt].mean() << " +- " << dev_dist[tt].standardDeviation() << std::endl;
  //   }
  // }

};



#endif
