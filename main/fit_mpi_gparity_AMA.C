//#define BOOST_SPIRIT_X3_DEBUG
#include <fstream>
#include <algorithm>
#include <minimizer.h>
#include <cost_function.h>
#include <fitfunc.h>
#include <fitfunc_mapping.h>
#include <random.h>
#include <plot.h>
#include <distribution.h>
#include <data_series.h>
#include <parser.h>
#include <common_defs.h>
#include <sstream>

#include <fit_mpi_gparity_AMA/fitfunc.h>
#include <fit_mpi_gparity_AMA/args.h>
#include <fit_mpi_gparity_AMA/read_data.h>
#include <fit_mpi_gparity_AMA/data_manipulations.h>
#include <fit_mpi_gparity_AMA/cmdline.h>
#include <fit_mpi_gparity_AMA/plot.h>

#ifdef USE_CALLGRIND
#include <valgrind/callgrind.h>
#else
#define CALLGRIND_START_INSTRUMENTATION
#define CALLGRIND_STOP_INSTRUMENTATION
#define CALLGRIND_DUMP_STATS
#endif

//Fit the pion mass from a simultaneous fit to multiple pseudoscalar two-point functions using AMA

int main(const int argc, const char** argv){
  Args args;
  if(argc < 2){
    std::ofstream of("template.args");
    of << args;
    std::cout << "Wrote template argument file to template.args\n";
    return 0;
  }
  
  const std::string arg_file = argv[1];
  CMDline cmdline(argc,argv,2);
    
  std::ifstream arg_f(arg_file.c_str());
  
  arg_f >> args;

  std::cout << "Read arguments: \n" << args << std::endl;

  const int ntraj = (args.traj_lessthan - args.traj_start)/args.traj_inc;
  assert(ntraj > 0);
  
  distributionMatrix PP_LW_FF_exact_data(args.Lt, distributionD(ntraj));
  distributionMatrix PP_LW_FF_sloppy_data(args.Lt, distributionD(ntraj));
  distributionMatrix PP_LW_BB_exact_data(args.Lt, distributionD(ntraj));
  distributionMatrix PP_LW_BB_sloppy_data(args.Lt, distributionD(ntraj));
  
#pragma omp parallel for
  for(int i=0;i<ntraj;i++){
    const int c = args.traj_start + i*args.traj_inc;
    read(PP_LW_FF_exact_data, PP_LW_FF_sloppy_data, args.PP_LW.FF_data, c, i);
    read(PP_LW_BB_exact_data, PP_LW_BB_sloppy_data, args.PP_LW.BB_data, c, i);  
  }
  PP_LW_BB_exact_data = timeReflect(PP_LW_BB_exact_data);
  PP_LW_BB_sloppy_data = timeReflect(PP_LW_BB_sloppy_data);

  distributionVector PP_LW_FF_sloppy_avg = sourceTimeSliceAverage(PP_LW_FF_sloppy_data);
  distributionVector PP_LW_BB_sloppy_avg = sourceTimeSliceAverage(PP_LW_BB_sloppy_data);

  std::cout << "Computing PP_LW_FF_correction\n";
  distributionVector PP_LW_FF_correction = computeAMAcorrection(PP_LW_FF_sloppy_data, PP_LW_FF_exact_data);
  std::cout << "Computing PP_LW_BB_correction\n";
  distributionVector PP_LW_BB_correction = computeAMAcorrection(PP_LW_BB_sloppy_data, PP_LW_BB_exact_data);
  publicationPrint printer;
  for(int t=0;t<args.Lt;t++){
    printer << "PP_LW_FF_sloppy_avg[" << t << "] = " << PP_LW_FF_sloppy_avg[t] << std::endl;
    printer << "PP_LW_BB_sloppy_avg[" << t << "] = " << PP_LW_BB_sloppy_avg[t] << std::endl;
    printer << "PP_LW_FF_correction[" << t << "] = " << PP_LW_FF_correction[t] << std::endl;
    printer << "PP_LW_BB_correction[" << t << "] = " << PP_LW_BB_correction[t] << std::endl;
  }

  distributionVector PP_LW = (PP_LW_FF_correction + PP_LW_FF_sloppy_avg + PP_LW_BB_correction + PP_LW_BB_sloppy_avg)/double(2.);
  for(int t=0;t<args.Lt;t++){
    printer << "PP_LW[" << t << "] = " << PP_LW[t] << std::endl;
  }

  const int ntypes = 1;
  const DataType type_map[1] = {PP_LW_data};
  const int nx = ntypes * args.Lt; //t + Lt*type
  typedef dataSeries<Coord, doubleJackknifeDistributionD> doubleJackknifeTimeSeries;
  typedef dataSeries<Coord, jackknifeDistributionD> jackknifeTimeSeries;
  
  doubleJackknifeTimeSeries data_dj(nx, ntraj);
  jackknifeTimeSeries data_j(nx, ntraj);
  for(int tt=0;tt<nx;tt++){
    const int t = tt % args.Lt;
    const DataType d = type_map[tt / args.Lt];    
    data_dj.coord(tt) = data_j.coord(tt) = Coord(t,d);
    data_dj.value(tt).resample(PP_LW[t]);
    data_j.value(tt).resample(PP_LW[t]);
    printer << "Resampled data " << data_j.coord(tt) << "  " << data_j.value(tt) << std::endl;
  }

  //Setup fit range
  typedef filteredDataSeries<doubleJackknifeTimeSeries> filteredDoubleJackknifeTimeSeries;
  typedef filteredDataSeries<jackknifeTimeSeries> filteredJackknifeTimeSeries;

  filterCoordTrange tfilter(args.t_min,args.t_max);
  
  filteredDoubleJackknifeTimeSeries inrange_data_dj(data_dj,tfilter);
  filteredJackknifeTimeSeries inrange_data_j(data_j,tfilter);
  
  tfilter.invert();

  filteredJackknifeTimeSeries outofrange_data_j(data_j,tfilter);
  
  const int ndata_fit = inrange_data_j.size();
    
  //Uncorrelated fit
  typedef NumericMatrix<jackknifeDistributionD> jackknifeMatrix;
  jackknifeMatrix cov(ndata_fit,jackknifeDistributionD(ntraj,0.));
  std::vector<jackknifeDistributionD> sigma(ndata_fit);
  NumericMatrix<double> corr(ndata_fit,0.);
  for(int x=0;x<ndata_fit;x++){
    cov(x,x) = doubleJackknifeDistributionD::covariance(data_dj.value(x), data_dj.value(x));
    sigma[x] = sqrt(cov(x,x));
    corr(x,x) = 1.;
  }

  //Load guesses if applicable
  AllFitParams guess;
  if(cmdline.load_guess){
    std::ifstream f(cmdline.guess_file.c_str());
    assert(f.good());
    f >> guess;
    f.close();
  }else{
    for(int i=1;i<guess.size();i++) guess(i) = 1e3;
    guess.m = 0.5;
  }
  
  std::cout << "Using guess: " << guess << std::endl;
  
  
  //Do the fit
  //FitMpi fitfunc(2*args.Lt); //FF and BB are cosh-like in 2*Lt
  FitMpiFrozen fitfunc(2*args.Lt);
  for(int i=2;i<AllFitParams::size();i++) fitfunc.freeze(i,-1);
  
  typedef sampleSeries<const filteredJackknifeTimeSeries> sampleSeriesConstType; //const access
  typedef UncorrelatedChisqCostFunction<decltype(fitfunc), sampleSeriesConstType, double, NumericVector<double> > CostFunctionType;
  typedef CostFunctionType::CostType CostType;

  //Setup minimizer
  typedef MarquardtLevenbergMinimizer<CostFunctionType> MinimizerType;

  MarquardtLevenbergParameters<CostType> mlparams;
  mlparams.output = &null_stream;

  jackknifeDistribution<AllFitParams> params(ntraj, guess);

#pragma omp parallel for
  for(int j=0;j<ntraj;j++){    
    sampleSeriesConstType dsample(inrange_data_j, j);
    std::vector<double> sigma_j(sigma.size()); for(int p=0;p<sigma.size();p++) sigma_j[p] = sigma[p].sample(j);
    CostFunctionType costfunc(fitfunc, dsample, sigma_j);

    MinimizerType fitter(costfunc, mlparams);
    
    AllFitParams &pj = params.sample(j);
    typename FitMpiFrozen::ParameterType pjr;
    fitfunc.reduce(pjr, pj); //get the subset of parameters that aren't frozen (those that are are stored internally when freeze is called above)
    CostType cost = fitter.fit(pjr);
    assert(fitter.hasConverged());
    fitfunc.expand(pj,pjr);
  }

  std::cout << "Params: " << params.mean() << " " << params.standardError() << std::endl;

  //Plot some nice things
  CALLGRIND_START_INSTRUMENTATION;
  jackknifeTimeSeriesD effmass = effectiveMass(data_j, PP_LW_data, 2*args.Lt);
  CALLGRIND_STOP_INSTRUMENTATION;
  CALLGRIND_DUMP_STATS;
  
  MatPlotLibScriptGenerate plotter;
  {
    typedef DataSeriesAccessor<jackknifeTimeSeriesD, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistributionD> > Accessor;
    typedef MatPlotLibScriptGenerate::handleType Handle;
    Accessor a(effmass);
    Handle ah = plotter.plotData(a);

    plotter.setLegend(ah, "PP_LW");
    plotter.createLegend();
    plotter.write("mpi_plot.py");
  }
    
  return 0;
}
