#include <fstream>
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

  //Generate some random data
  int npoints = 10;
  int nsample = 100;

  publicationPrint printer;

  typedef dataSeries<double, distribution<double> > RawDataSeriesType;
  RawDataSeriesType data(npoints, nsample);
  
  for(int i=0;i<npoints;i++){
    data.coord(i) = i;
    for(int j=0;j<nsample;j++){
      data.value(i).sample(j) = gaussianRandom<double>(i, 0.3);
    }
    std::cout << "Mean " << i << " " << data.value(i).mean() << " Std. Dev. " << data.value(i).standardDeviation() << " " << " Std. Err. " << data.value(i).standardError() << std::endl;
  }
  
  //Setup fit range
  double x_min = 2.0, x_max = 8.0;
  
  typedef filteredDataSeries<RawDataSeriesType> InRangeDataSeriesType;
  filterXrange<double> xfilter(x_min,x_max);

  InRangeDataSeriesType inrange_data(data,xfilter);
  xfilter.invert();
  InRangeDataSeriesType outofrange_data(data,xfilter);
  
  typedef dataSeries<double, jackknifeDistribution<double> > JackknifeSeriesType;
  
  JackknifeSeriesType inrange_jackknife;   resample(inrange_jackknife, inrange_data);
  JackknifeSeriesType outofrange_jackknife;   resample(outofrange_jackknife, outofrange_data);

  const int ndata_in_range = inrange_jackknife.size();
  
  std::cout << "Resampled data in range:\n";
  
  for(int i=0;i<ndata_in_range;i++){
    std::cout << "Standard error " << i << " " << inrange_jackknife.value(i).standardError() << std::endl;    
  }
  
  //Setup fit function
  typedef NumericLinearFit<double, double, 1> FitFunctionType; //a + bx
  FitFunctionType func;

  typedef FrozenFitFunc<FitFunctionType> FrozenFitFunction;
  
  //Setup cost function, acts on individual jackknife samples so need accessor
  std::vector<double> sigma(ndata_in_range);
  for(int i=0;i<sigma.size();i++)
    sigma[i] = inrange_jackknife.value(i).standardError();

  //Compute double-jackknife covariance matrix
  typedef dataSeries<double, doubleJackknifeDistribution<double> > doubleJackknifeSeriesType;
  doubleJackknifeSeriesType inrange_doublejacknife; resample(inrange_doublejacknife, inrange_data);

  printer.setMinWidth(20);

  typedef NumericMatrix<jackknifeDistribution<double> >  jackknifeMatrixD;
  jackknifeMatrixD cov(ndata_in_range);
  for(int i=0;i<ndata_in_range;i++)
    for(int j=0;j<ndata_in_range;j++)
      cov(i,j) = doubleJackknifeDistribution<double>::covariance(inrange_doublejacknife.value(i), inrange_doublejacknife.value(j));
  printer << "double Jackknife covariance matrix:\n" << cov;

  //Normalize covariance matrix to get correlation matrix
  jackknifeMatrixD corr(ndata_in_range);
  for(int i=0;i<ndata_in_range;i++)
    for(int j=0;j<ndata_in_range;j++)
      corr(i,j) = cov(i,j)/sqrt(cov(i,i))/sqrt(cov(j,j));
  printer << "double Jackknife correlation matrix:\n" << corr;
  
  //Generate inverse correlation matrix
  jackknifeMatrixD inv_corr(ndata_in_range, jackknifeDistribution<double>(nsample));
  svd_inverse(inv_corr, corr);

  printer << "double Jackknife inverse correlation matrix:\n" << inv_corr;

  //Test the inverse
  jackknifeMatrixD inv_test(ndata_in_range, jackknifeDistribution<double>(nsample));
  for(int i=0;i<ndata_in_range;i++){
    for(int j=0;j<ndata_in_range;j++){
      inv_test(i,j) = jackknifeDistribution<double>(nsample,0.);
      for(int k=0;k<ndata_in_range;k++)
	inv_test(i,j) = inv_test(i,j) + corr(i,k)*inv_corr(k,j);
    }
  }
  printer << "inverse test:\n" << inv_test;

  
  typedef sampleSeries<const JackknifeSeriesType> SampleSeriesConstType; //const access
  
  //typedef UncorrelatedChisqCostFunction<FrozenFitFunction, SampleSeriesConstType> CostFunctionType;
  typedef NumericMatrixSampleView<const jackknifeMatrixD> InvCorrMatrixViewType;
  typedef CorrelatedChisqCostFunction<FrozenFitFunction, SampleSeriesConstType, InvCorrMatrixViewType> CostFunctionType;
  
  typedef CostFunctionType::CostType CostType;

  //Setup minimizer
  typedef MarquardtLevenbergMinimizer<CostFunctionType> MinimizerType;

  MarquardtLevenbergParameters<CostType> mlparams;
  mlparams.output = &null_stream;
  
  std::vector<jackknifeDistribution<double> > param(2);
  param[0] = jackknifeDistribution<double>(nsample, 5);
  param[1] = jackknifeDistribution<double>(nsample, 0.5);

#pragma omp parallel for
  for(int j=0;j<nsample;j++){
    FrozenFitFunction func_frozen(func); //provide defaults (can freeze here too)
    //func_frozen.freeze(0, param[0].sample(j));
    
    SampleSeriesConstType dsample(inrange_jackknife, j);

    //CostFunctionType costfunc(func_frozen, dsample, sigma);

    InvCorrMatrixViewType inv_corr_j(inv_corr,j);
    CostFunctionType costfunc(func_frozen, dsample, sigma, inv_corr_j);
    
    MinimizerType fitter(costfunc, mlparams);

    NumericVector<double> dp(2);
    dp(0) = param[0].sample(j);
    dp(1) = param[1].sample(j);

    NumericVector<double> dp_sub = func_frozen.mapParamsSupersetToSubset(dp);

    CostType cost = fitter.fit(dp_sub);
    assert(fitter.hasConverged());
    
    //std::cout << "Converged in " << fitter.iterations() << " iterations with final chisq " << cost << " chisq/dof " <<  cost/costfunc.Ndof() << std::endl;

    dp = func_frozen.mapParamsSubsetToSuperset(dp_sub);
    
    param[0].sample(j) = dp[0];
    param[1].sample(j) = dp[1];
  }

  printf("Parameters %g +- %g,   %g +- %g\n",param[0].mean(), param[0].standardError(), param[1].mean(), param[1].standardError());

  //Write out parameters
  std::ofstream ofs("parameters.dat");
  {
    boost::archive::text_oarchive oa(ofs);
    oa << param;
  }
    
  //Plot result
  MatPlotLibInterface mpl;
  typedef DataSeriesAccessor<JackknifeSeriesType, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistribution<double> > > PlotDataAccessor;
  {
    PlotDataAccessor datainterface(inrange_jackknife);
    mpl.plotData(datainterface);
  }
  {
    PlotDataAccessor datainterface(outofrange_jackknife);
    boost::python::dict args; args["hollowsymbol"] = 1;
    mpl.plotData(datainterface,args);
  }
  

  //Generate fit curve
  int npt = 100;
  JackknifeSeriesType fit(npt,nsample);
  double delta = (x_max - x_min)/(npt-1);
#pragma omp parallel for
  for(int i=0;i<npt;i++){
    fit.coord(i) = x_min + i*delta;
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

  typedef DataSeriesAccessor<JackknifeSeriesType, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistribution<double> > > PlotCurveAccessor;
  PlotCurveAccessor curveinterface(fit);
  mpl.errorBand(curveinterface,args);
  
  //mpl.plotData(curveinterface);
  
  
  mpl.write("test.pdf");
  
  std::cout << "Normal exit\n"; std::cout.flush();
  return 0;
}

