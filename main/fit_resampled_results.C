//Perform a fit to a set of jackknife resampled data

#include<parser.h>
#include<common.h>
#include<fit.h>
#include<plot.h>

using namespace CPSfit;

#define DATA_INFO_MEMBERS				\
  ( std::string, filename )				\
  ( std::vector<int>, input_idx)			\
  ( double, coord )

struct DataInfo{
  GENERATE_MEMBERS(DATA_INFO_MEMBERS);

  DataInfo(): input_idx(1,0), coord(0), filename("file.hdf5"){}
};
GENERATE_PARSER(DataInfo, DATA_INFO_MEMBERS);

GENERATE_ENUM_AND_PARSER(FitFuncType, (FParabola) );

#define ARGS_MEMBERS				\
  ( std::vector<DataInfo>, data )		\
  ( FitFuncType, fitfunc)

struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS);

  Args(): data(1), fitfunc(FParabola){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS);

class FitParabola{
public:
  typedef double ValueType;
  typedef NumericVector<double> ParameterType;
  typedef NumericVector<double> ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate;

  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
    double dt = t-p(2);
    return p(0) + p(1)*dt*dt;
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    double dt = t-p(2);
    ValueDerivativeType yderivs(3);
    yderivs(0) = 1.;
    yderivs(1) = dt*dt;
    yderivs(2) = -2*p(1)*dt;
    return yderivs;
  }

  inline int Nparams() const{ return 3; }

  inline ParameterType basicGuess() const{ ParameterType out(3); out(0) = out(1) = out(2) = 1; return out; }
};

namespace CPSfit{
inline NumericVector<double> operator*(const NumericVector<double> &a, const NumericVector<double> &b){
  return NumericVector<double>(a.size(), [&](const int i){ return a(i)*b(i); });
}
};




template<typename FitFunc>
void fitCentralSpecFF(const jackknifeCorrelationFunctionD &data, const FitFunc &fitfunc, const Args &args){
  typedef correlationFunction<double,double> CentralValueCorrelationFunction;

  const int ndata = data.size();
  CentralValueCorrelationFunction data_cen(ndata, [&](const int i){ return CentralValueCorrelationFunction::ElementType(data.coord(i), data.value(i).best()); });
  
  std::vector<double> sigma(ndata,1.);
  //for(int i=0;i<ndata;i++) sigma[i] = data.value(i).standardError();
  
  typedef UncorrelatedChisqCostFunction<FitFunc, CentralValueCorrelationFunction> CostFunctionType;
  CostFunctionType cost_func(fitfunc, data_cen, sigma);
  
  MarquardtLevenbergParameters<double> mlparams;
  mlparams.delta_cost_min = 1e-14;
  mlparams.verbose = true;

  MarquardtLevenbergMinimizer<CostFunctionType> minimizer(cost_func, mlparams);
  
  typename FitFunc::ParameterType params = fitfunc.basicGuess();

  params(0) = 23;
  params(1) = 0.02;
  params(2) = 508;

  double chisq = minimizer.fit(params);
  
  std::cout << "Chisq " << chisq << std::endl;
  std::cout << "Params " << params << std::endl;

  int nfitpts = 60;
  double delta = (data_cen.coord(ndata-1) - data_cen.coord(0))/(nfitpts-1);

  CentralValueCorrelationFunction fit_pred(nfitpts, 
					   [&](const int i){ 
					     double coor = data_cen.coord(0) + delta*i;
					     return CentralValueCorrelationFunction::ElementType(coor, 
												 fitfunc.value(coor,params));
					   }
					   );



  MatPlotLibScriptGenerate plot;
  typedef DataSeriesAccessor<CentralValueCorrelationFunction,ScalarCoordinateAccessor<double> , ScalarValueAccessor<double> > accessor;

  MatPlotLibScriptGenerate::kwargsType kwargs;

  accessor dacc(data_cen);
  kwargs["marker"] = "o";
  plot.plotData(dacc, kwargs);
  

  //accessor facc(fit_pred);  

  struct FitLine{
    const CentralValueCorrelationFunction &data;
    
    FitLine(const CentralValueCorrelationFunction &_data): data(_data){}
    inline double x(const int i) const{ return data.coord(i); }
    inline double upper(const int i) const{ return data.value(i); }
    inline double lower(const int i) const{ return data.value(i); }    
    inline int size() const{ return data.size(); }
  };
  kwargs.clear();
  kwargs["alpha"] = 0.5;
  FitLine facc(fit_pred);
  plot.errorBand(facc,kwargs);

  //kwargs["marker"] = "s";
  //plot.plotData(facc, kwargs);
  
  plot.write("fit_central.py","fit_central.pdf");
}







template<typename FitFunc>
void fitSpecFF(const jackknifeCorrelationFunctionD &data, const FitFunc &fitfunc, const Args &args){
  typedef typename composeFitPolicy<FitFunc,standardFitFuncPolicy,correlatedFitPolicy>::type FitPolicies;
  
  const int nsample = data.value(0).size();
  const int ndata = data.size();
  
  NumericSquareMatrix<double> cov(ndata, [&](const int i, const int j){ return jackknifeDistributionD::covariance(data.value(i),data.value(j)); });
  NumericVector<jackknifeDistributionD> sigma_j(ndata, [&](const int i){ return jackknifeDistributionD(nsample,sqrt(cov(i,i))); });
  
  NumericSquareMatrix<double> corr(ndata, [&](const int i, const int j){ return cov(i,j)/sqrt(cov(i,i)*cov(j,j)); });
  NumericSquareMatrix<double> inv_corr(ndata);
  svd_inverse(inv_corr,corr);

  //NumericSquareMatrix<jackknifeDistributionD> inv_corr_j(ndata, [&](const int i, const int j){ return jackknifeDistributionD(nsample, inv_corr(i,j)); });

  NumericSquareMatrix<jackknifeDistributionD> inv_corr_j(ndata, [&](const int i, const int j){ return jackknifeDistributionD(nsample, i==j ? 1. : 0.); });


  std::cout << "Correlation matrix:\n" << corr << std::endl;
  std::cout << "Sigma: " << sigma_j << std::endl;

  MarquardtLevenbergParameters<double> mlparams;
  mlparams.delta_cost_min = 1e-10;

  fitter<FitPolicies> fitter(mlparams);
  fitter.importCostFunctionParameters(inv_corr_j,sigma_j);
  fitter.importFitFunc(fitfunc);
  
  typename FitFunc::ParameterType guess = fitfunc.basicGuess();

  guess(0) = 20;
  guess(1) = -1e-5;
  guess(2) = 500;

  jackknifeDistribution<typename FitFunc::ParameterType> params(nsample,guess);
  
  jackknifeDistributionD chisq, chisq_per_dof;
  fitter.fit(params,chisq,chisq_per_dof,data);
  
  std::cout << "Fit results\n";
  std::cout << "Params: " << params << std::endl;
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
}

void fit(const jackknifeCorrelationFunctionD &data, const Args &args){
  if(args.fitfunc == FParabola){
    FitParabola fitfunc;
    fitCentralSpecFF<FitParabola>(data,fitfunc,args);       
    fitSpecFF<FitParabola>(data,fitfunc,args);
  }else{
    assert(0);
  }
}
	 
int main(const int argc, const char *argv[]){
  Args args;
  if(argc < 2){
    std::ofstream of("template.args");
    (std::cout << "No parameter file provided: writing template to 'template.args' and exiting\n").flush();
    of << args;
    return 1;
  }    
  
  parse(args, argv[1]);
  
  const int ndata = args.data.size();
  jackknifeCorrelationFunctionD data(ndata);
  int nsample;
  for(int i=0;i<ndata;i++){
    data.coord(i) = args.data[i].coord;
    readHDF5file(data.value(i), args.data[i].filename, args.data[i].input_idx);
    if(i==0) nsample = data.value(i).size();
    else if(data.value(i).size() != nsample) error_exit(std::cout << "jackknife distribution from file " << args.data[i].filename << " has " << data.value(i).size() << " samples, expected "<< nsample << std::endl);

    std::cout << data.coord(i) << " " << data.value(i) << std::endl;
  }
  
  fit(data,args);
  
  return 0;
}


