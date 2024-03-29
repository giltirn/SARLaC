#pragma once

class covMatStrategyBase{
public:
  virtual NumericSquareMatrix<double> compute(const correlationFunction<double, rawDataDistributionD> &data) const = 0; //when not doing a fit
  virtual void compute(simpleSingleFitWrapper &fitter, const correlationFunction<double, rawDataDistributionD> &data) const{ fitter.importCovarianceMatrix(this->compute(data)); }
  virtual ~covMatStrategyBase(){};
};
class covMatStrategyCorrelated: public covMatStrategyBase{
public:  
  NumericSquareMatrix<double> compute(const correlationFunction<double, rawDataDistributionD> &data) const override{
    int Lt = data.size();
    NumericSquareMatrix<double> cov(Lt);
    for(int t1=0;t1<Lt;t1++)
      for(int t2=0;t2<Lt;t2++)
	cov(t1,t2) = rawDataDistributionD::covariance( data.value(t1), data.value(t2) );
    return cov;
  }
};
class covMatStrategyUncorrelated: public covMatStrategyBase{
public:
  NumericSquareMatrix<double> compute(const correlationFunction<double, rawDataDistributionD> &data) const override{
    int Lt = data.size();
    NumericSquareMatrix<double> cov(Lt, 0.);
    for(int t=0;t<Lt;t++)
      cov(t,t) = rawDataDistributionD::covariance( data.value(t), data.value(t) );
    return cov;
  }  
};

#define CMAT_PARAMS_CUTOFF (double, cutoff)
struct covMatStrategyCutoffArgs{
  GENERATE_MEMBERS(CMAT_PARAMS_CUTOFF); 
  covMatStrategyCutoffArgs(): cutoff(0.005){  }
};
GENERATE_PARSER( covMatStrategyCutoffArgs, CMAT_PARAMS_CUTOFF);

//cf https://arxiv.org/pdf/1101.2248.pdf
//Remove contribution of eigenvectors
class covMatStrategyCutoff: public covMatStrategyBase{
public:
  double cutoff; 

  covMatStrategyCutoff(const std::string &args_file){
    covMatStrategyCutoffArgs args;
    parseOrTemplate(args, args_file, "cov_cutoff_template.args");
    cutoff = args.cutoff;
  }
  
  NumericSquareMatrix<double> compute(const correlationFunction<double, rawDataDistributionD> &data) const override{
    std::cout << "WARNING: Only the inverse covariance matrix is defined for covMatStrategyCutoff (matrix is singular)" << std::endl;
    return NumericSquareMatrix<double>(data.size(),0.);
  }

  void compute(simpleSingleFitWrapper &fitter, const correlationFunction<double, rawDataDistributionD> &data) const override{    
    int Lt = data.size();
    NumericSquareMatrix<double> cov(Lt);
    for(int t1=0;t1<Lt;t1++)
      for(int t2=0;t2<Lt;t2++)
	cov(t1,t2) = rawDataDistributionD::covariance( data.value(t1), data.value(t2) );

    //Remove the contribution of eigenvalues L < cutoff from covariance matrix
    std::vector< NumericVector<double> > evecs(Lt, NumericVector<double>(Lt) );
    std::vector<double> evals(Lt);
    
    GSLsymmEigenSolver< NumericVector<double>, NumericSquareMatrix<double> >::symmetricMatrixSolve(evecs, evals, cov);
    
    std::cout << "Cutoff at " << cutoff << std::endl;

    NumericSquareMatrix<double> inv_cov(Lt, 0.);
    std::cout << "Evals : ";
    for(int l=0;l<Lt;l++){
      std::cout << evals[l];
      if(evals[l]>cutoff){
	std::cout << "*";

	for(int t=0;t<Lt;t++)
	  for(int u=0;u<Lt;u++)
	    inv_cov(t,u) += evecs[l](t)*evecs[l](u)*(1./evals[l]);
      }
      std::cout << " ";
    }
    std::cout << std::endl;

    std::vector<double> sigma_dummy(Lt,1.0);
    fitter.importInverseCorrelationMatrix(inv_cov, sigma_dummy);
  }  
};



//cf https://arxiv.org/pdf/1101.2248.pdf
//Modified covariance matrix
class covMatStrategyMCM: public covMatStrategyBase{
public:
  double cutoff; 

  covMatStrategyMCM(const std::string &args_file){
    covMatStrategyCutoffArgs args;
    parseOrTemplate(args, args_file, "cov_cutoff_template.args");
    cutoff = args.cutoff;
  }
  
  NumericSquareMatrix<double> compute(const correlationFunction<double, rawDataDistributionD> &data) const override{
    int Lt = data.size();
    NumericSquareMatrix<double> cov(Lt);
    for(int t1=0;t1<Lt;t1++)
      for(int t2=0;t2<Lt;t2++)
	cov(t1,t2) = rawDataDistributionD::covariance( data.value(t1), data.value(t2) );
    
    std::vector<double> sigma(Lt);
    NumericSquareMatrix<double> corr(Lt);

    for(int t=0;t<Lt;t++) sigma[t] = sqrt(cov(t,t));
    for(int t=0;t<Lt;t++)
      for(int u=0;u<Lt;u++)
	corr(t,u) = cov(t,u)/sigma[t]/sigma[u];

    std::vector< NumericVector<double> > evecs(Lt, NumericVector<double>(Lt) );
    std::vector<double> evals(Lt);
    
    GSLsymmEigenSolver< NumericVector<double>, NumericSquareMatrix<double> >::symmetricMatrixSolve(evecs, evals, cov);

    std::cout << "Cutoff at " << cutoff << std::endl;

    corr.zero();
    std::cout << "Evals : ";
    for(int l=0;l<Lt;l++){
      std::cout << evals[l];
      if(evals[l]>cutoff){
	std::cout << "*";

	for(int t=0;t<Lt;t++)
	  for(int u=0;u<Lt;u++)
	    corr(t,u) += evecs[l](t)*evecs[l](u)*evals[l];
      }
      std::cout << " ";
    }
    std::cout << std::endl;

    //Restore diagonal entries    
    for(int t=0;t<Lt;t++) corr(t,t) = 1.0;
    
    for(int t=0;t<Lt;t++)
      for(int u=0;u<Lt;u++)
	cov(t,u) = corr(t,u) * sigma[t] * sigma[u];

    return cov;
  }

};





//args_file need only be provided for those strategies with tunable arguments
std::unique_ptr<covMatStrategyBase> covMatStrategyFactory(CovMatStrategy strat, const std::string &args_file){
  if(strat == CovMatStrategy::Correlated){
    return std::unique_ptr<covMatStrategyBase>(new covMatStrategyCorrelated);
  }else if(strat == CovMatStrategy::Uncorrelated){
    return std::unique_ptr<covMatStrategyBase>(new covMatStrategyUncorrelated);
  }else if(strat == CovMatStrategy::Cutoff){
    return std::unique_ptr<covMatStrategyBase>(new covMatStrategyCutoff(args_file));
  }else if(strat == CovMatStrategy::MCM){
    return std::unique_ptr<covMatStrategyBase>(new covMatStrategyMCM(args_file));
  }else{
    error_exit(std::cout << "Invalid covariance matrix strategy" << std::endl);
  }
}
