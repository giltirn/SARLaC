#include <minimizer.h>
#include <cost_function.h>
#include <fitfunc.h>
#include <fitfunc_mapping.h>
#include <random.h>
#include <plot.h>
#include <distribution.h>
#include <data_series.h>


int main(void){
  RNG.initialize(1234);

  int npoints = 10;
  int nsample = 10;
  typedef dataSeries<double, distribution<double> > RawDataSeriesType;
  RawDataSeriesType data(npoints, nsample);

  std::vector<double> sigma(npoints);
  
  for(int i=0;i<npoints;i++){
    data.coord(i) = i;
    for(int j=0;j<nsample;j++){
      data.value(i).sample(j) = gaussianRandom<double>(i, 0.3);
    }
    std::cout << "Mean " << i << " " << data.value(i).mean() << " Std. Dev. " << data.value(i).standardDeviation() << " " << " Std. Err. " << data.value(i).standardError() << std::endl;

    sigma[i] = data.value(i).standardError(); //use raw standard errors as weights
  }

  typedef dataSeries<double, jackknifeDistribution<double> > JackknifeSeriesType;
  JackknifeSeriesType jackknife(npoints);
  for(int i=0;i<npoints;i++){
    jackknife.coord(i) = data.coord(i);
    jackknife.value(i).resample(data.value(i));
  }
    
  std::cout << "Resampled data:\n";
  
  for(int i=0;i<npoints;i++){
    std::cout << "Standard error " << i << " " << jackknife.value(i).standardError() << std::endl;    
  }

  typedef NumericLinearFit<double, double, 1> FitFunctionType; //a + bx
  FitFunctionType func;
  
  //Implement frozen fit through generic interface, mapping between a subset and superset of the parameters
  typedef NumericVector<double> ParamSupersetType;
  typedef NumericVector<double> ParamSubsetType;
  typedef NumericVector<double> DerivSupersetType;
  typedef NumericVector<double> DerivSubsetType;
  
  typedef FitFuncGenericRemap<FitFunctionType> FrozenFitFunction;
  
  //Setup cost function
  typedef typename JackknifeSeriesType::sampleAccessorType SampleSeriesType; //unconst access
  typedef typename JackknifeSeriesType::constSampleAccessorType SampleSeriesConstType; //const access
  
  typedef UncorrelatedChisqCostFunction<FrozenFitFunction, SampleSeriesConstType> CostFunctionType;
  typedef CostFunctionType::CostType CostType;
  
  typedef MarquardtLevenbergMinimizer<CostFunctionType> MinimizerType;

  MarquardtLevenbergParameters<CostType> mlparams;
  //mlparams.output = &null_stream;
  
  std::vector<jackknifeDistribution<double> > param(2);
  param[0] = jackknifeDistribution<double>(nsample, 0);
  param[1] = jackknifeDistribution<double>(nsample, 0.5);
  
  for(int j=0;j<nsample;j++){
    FrozenFitFunction func_frozen(func, ParamSupersetType(2), DerivSubsetType(1)); //provide defaults (can freeze here too)
    
    func_frozen.addParamMapSubsetToSuperset(ArrayElementMap<ParamSupersetType,ParamSubsetType>(1,0)); //sub 0-> super 1
    func_frozen.addParamMapSubsetToSuperset(ArrayElementFreeze<double,ParamSupersetType,ParamSubsetType>(0, param[0].sample(j))); //super 0=const from input

    func_frozen.addDerivMapSupersetToSubset(ArrayElementMap<ParamSupersetType,ParamSubsetType>(0,1)); //super 1-> sub 0
    
    SampleSeriesConstType dsample(jackknife, j);

    CostFunctionType costfunc(func_frozen, dsample, sigma);

    MinimizerType fitter(costfunc, mlparams);

    // NumericVector<double> dp(2);
    // dp(0) = param[0].sample(j);
    // dp(1) = param[1].sample(j);

    NumericVector<double> dp(1);
    dp[1] = param[1].sample(j);
    
    CostType cost = fitter.fit(dp);
    assert(fitter.hasConverged());
    
    std::cout << "Converged in " << fitter.iterations() << " iterations with final chisq " << cost << " chisq/dof " <<  cost/costfunc.Ndof() << std::endl;
    
    // param[0].sample(j) = dp[0];
    // param[1].sample(j) = dp[1];

    param[1].sample(j) = dp[0];
  }

  //Plot result
  MatPlotLibInterface mpl;
  typedef DataSeriesAccessor<JackknifeSeriesType, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistribution<double> > > PlotDataAccessor;
  PlotDataAccessor datainterface(jackknife);
  mpl.plotData(datainterface);

  //Generate fit curve
  int npt = 100;
  JackknifeSeriesType fit(npt,nsample);
  for(int i=0;i<npt;i++){
    fit.coord(i) = double(i)/10;
    for(int j=0;j<nsample;j++){
      NumericVector<double> fp(2);
      fp(0) = param[0].sample(j);
      fp(1) = param[1].sample(j);
      
      fit.value(i).sample(j) = func.value(fit.coord(i),fp);
    }
  }
  boost::python::dict args;
  args["alpha"] = 0.3;
  args["color"] = "r";
  
  PlotDataAccessor curveinterface(fit);
  mpl.errorBand(curveinterface,args);
  
  //mpl.plotData(curveinterface);
  
  
  mpl.write("test.pdf");

  std::cout << "Normal exit\n"; std::cout.flush();
  return 0;
}
