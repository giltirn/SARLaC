#include <random.h>
#include <plot.h>
#include <distribution.h>
#include <data_series.h>
#include <fit.h>
#include <tensors.h>
#include <containers/parameter_vector.h>

using namespace SARLaC;

int main(void){
  RNG.initialize(1234);

  int npoints = 10;
  int nsample = 10;
  typedef dataSeries<double, rawDataDistribution<double> > RawDataSeriesType;
  RawDataSeriesType data(npoints, nsample);

  std::vector<double> sigma(npoints);

  std::cout << "Data generated by Gaussian distribution along line  y = x\n";
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

  typedef FitFuncLinearMultiDim<double, double, 1> FitFunctionType; //a + bx
  FitFunctionType func;
  
  //Implement frozen fit through generic interface, mapping between a subset and superset of the parameters
  typedef parameterVector<double> ParamSupersetType;
  typedef parameterVector<double> ParamSubsetType;
  typedef parameterVector<double> DerivSupersetType;
  typedef parameterVector<double> DerivSubsetType;
  
  typedef FitFuncGenericRemap<FitFunctionType> FrozenFitFunction;
  
  //Setup cost function
  typedef sampleSeries<JackknifeSeriesType> SampleSeriesType;
  typedef sampleSeries<const JackknifeSeriesType> SampleSeriesConstType;

  typedef UncorrelatedChisqCostFunction<FrozenFitFunction, SampleSeriesConstType> CostFunctionType;
  typedef CostFunctionType::CostType CostType;
  
  typedef MarquardtLevenbergMinimizer<CostFunctionType> MinimizerType;

  MarquardtLevenbergParameters<CostType> mlparams;
  
  std::vector<jackknifeDistribution<double> > param(2);
  param[0] = jackknifeDistribution<double>(nsample, 0.);
  param[1] = jackknifeDistribution<double>(nsample, 0.5);

  std::cout << "2 parameters\n";
  std::cout << "param[0] is frozen, with value " << param[0] << std::endl;
  std::cout << "param[1] is varying, with guess " << param[1] << std::endl;
  
  for(int j=0;j<nsample;j++){
    FrozenFitFunction func_frozen(func, ParamSupersetType(2), DerivSubsetType(1)); //provide defaults (can freeze here too)
    
    func_frozen.addParamMapSubsetToSuperset(ArrayElementMap<ParamSupersetType,ParamSubsetType>(1,0)); //sub 0-> super 1
    func_frozen.addParamMapSubsetToSuperset(ArrayElementFreeze<double,ParamSupersetType,ParamSubsetType>(0, param[0].sample(j))); //super 0=const from input

    func_frozen.addDerivMapSupersetToSubset(ArrayElementMap<ParamSupersetType,ParamSubsetType>(0,1)); //super 1-> sub 0
    
    SampleSeriesConstType dsample(jackknife, j);

    CostFunctionType costfunc(func_frozen, dsample, sigma);

    MinimizerType fitter(costfunc, mlparams);

    parameterVector<double> dp(1);
    dp[1] = param[1].sample(j);
    
    CostType cost = fitter.fit(dp);
    assert(fitter.hasConverged());
    
    std::cout << "Converged in " << fitter.iterations() << " iterations with final chisq " << cost << " chisq/dof " <<  cost/costfunc.Ndof() << std::endl;

    param[1].sample(j) = dp[0];
  }

  std::cout << "Result:\n";
  std::cout << "param[0] " << param[0] << std::endl;
  std::cout << "param[1] " << param[1] << std::endl;
  
  //Plot result
  MatPlotLibScriptGenerate mpl;  
  typedef DataSeriesAccessor<JackknifeSeriesType, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistribution<double> > > PlotDataAccessor;
  PlotDataAccessor datainterface(jackknife);
  mpl.plotData(datainterface);

  //Generate fit curve
  int npt = 100;
  JackknifeSeriesType fit(npt,nsample);
  for(int i=0;i<npt;i++){
    fit.coord(i) = double(i)/10;
    for(int j=0;j<nsample;j++){
      parameterVector<double> fp(2);
      fp(0) = param[0].sample(j);
      fp(1) = param[1].sample(j);
      
      fit.value(i).sample(j) = func.value(fit.coord(i),fp);
    }
  }
  typename MatPlotLibScriptGenerate::kwargsType args;
  args["alpha"] = 0.3;
  args["color"] = "r";
  
  PlotDataAccessor curveinterface(fit);
  mpl.errorBand(curveinterface,args);
  
  //mpl.plotData(curveinterface);
  
  
  mpl.write("test.py");

  std::cout << "Normal exit\n"; std::cout.flush();
  return 0;
}
