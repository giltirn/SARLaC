#include <fstream>
#include <sstream>

#include<fit.h>
#include<random.h>
#include<plot.h>
#include<distribution.h>
#include<data_series.h>
#include<parser.h>
#include<containers.h>
#include<common.h>

using namespace SARLaC;

#define TEST_FIT_ARGS_MEMBERS \
  ( int, nsample ) \
  ( int, npoints ) \
  ( double, x_min ) \
  ( double, x_max )

struct TestFitArgs{
  GENERATE_MEMBERS(TEST_FIT_ARGS_MEMBERS)  //you don't have to use this, instead you can define the members separately. However seeing as you have gone to the trouble of writing the ...MEMBERS definition you may as well...

  //other methods here if you want
};
GENERATE_PARSER(TestFitArgs, TEST_FIT_ARGS_MEMBERS)



int main(void){  
  RNG.initialize(1234);

  std::string fit_args_in = 
"{"
"   nsample = 100"
"   npoints = 10"
"   x_min = 2.0"
"   x_max = 8.0"
"}";

  std::istringstream is(fit_args_in); //normally you would do this from an external file!
  is >> std::noskipws;
  boost::spirit::istream_iterator fiter(is);
  
  TestFitArgs fit_args;
  fiter >> fit_args;

  std::cout << "Read arguments: \n" << fit_args << std::endl;

  int npoints = fit_args.npoints;
  int nsample = fit_args.nsample;

  ///////////////////////////////////////////////////////////
  //Generate some random data
  typedef correlationFunction<double, rawDataDistributionD > RawDataSeriesType;
  RawDataSeriesType data(npoints, nsample);
  
  for(int i=0;i<npoints;i++){
    data.coord(i) = i;
    for(int j=0;j<nsample;j++){ //random values for each data sample
      data.value(i).sample(j) = gaussianRandom<double>(1.0, 0.3);
    }
    std::cout << "Mean " << i << " " << data.value(i).mean() << " Std. Dev. " << data.value(i).standardDeviation() << " " << " Std. Err. " << data.value(i).standardError() << std::endl;
  }

  ///////////////////////////////////////////////////////////
  //Boostrap resample the data. Here I will use the version with an explicit resample table for the sake of clarity (you can also let it generate this internally using a fixed, common seed)
  int nboot = 800; //number of bootstrap resamplings
  std::vector<std::vector<int> > rtable = resampleTable(RNG, nsample, nboot);

  //In order to do a formally "correct" fit, we should generate a covariance matrix separately for each bootstrap ensemble. This can be achieved conveniently using the bootJackknife class
  //which combines an outer bootstrap resampling to generate the resamples ensembled, and an inner jackknife resampling to compute the covariance matrix
  typedef correlationFunction<double, bootstrapDistributionD > ResampledDataSeriesType;
  typedef correlationFunction<double, bootJackknifeDistributionD > ResampledCovSeriesType;
  ResampledDataSeriesType data_r(npoints);
  ResampledCovSeriesType data_cov(npoints);
  for(int p=0;p<npoints;p++){
    data_r.coord(p) = data_cov.coord(p) = data.coord(p);
    data_r.value(p).resample(data.value(p), rtable);
    data_cov.value(p).resample(data.value(p), rtable);

    std::cout << "Bootstrap std.err. " << p << " " << data_r.value(p).standardError() << std::endl; //should be in good agreement with number computed from sample standard dev. above 
  }
 
  ////////////////////////////////////////////////////////////
  //Setup fit function. It needs to be wrapped for the simple fitter
  FitConstant fitfunc;
  genericFitFuncWrapper<FitConstant> fwrap(fitfunc, fitfunc.guess()); //the 2nd parameter tells the wrapper how to initialize an object containing the parameters of the wrapped fit function (contents uninmportant, only size!)

  //Setup fitter
  MinimizerType min_type = MinimizerType::MarquardtLevenberg;
  simpleFitWrapper<bootstrapDistributionD> fitter(fwrap, min_type);

  //We can use the fitter's in-built functionality to compute the covariance matrices from our bootJackknife data
  fitter.generateCovarianceMatrix(data_cov);

  /////////////////////////////////////////////////////////////
  //Now we are ready to fit. We need to prepare our outputs. The parameters should be initialized to a good "guess" from which the minimizer starts
  auto binit = data_r.value(0).getInitializer(); //convenient way of getting the initializer for bootstraps that tells them how to resize
  bootstrapDistributionD guess(binit, 0.5); //guess is 0.5

  std::vector<bootstrapDistributionD> params(1, guess); //only 1 parameter!
  bootstrapDistributionD chisq(binit), chisq_per_dof(binit);
  int dof;
  fitter.fit(params, chisq, chisq_per_dof, dof, data_r);

  std::cout << "Fit result " << params[0] << std::endl;
  
  ///////////////////////////////////////////////////////////////
  //Store outputs
  writeParamsStandard(params, "params.hdf5"); //write output in a format that hdf5_print, hdf5_calc understands
  writeParamsStandard(chisq, "chisq.hdf5");
  writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");

  std::cout << "Done" << std::endl;
  return 0;
}
