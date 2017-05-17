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

  distributionVector PP_LW_raw = readCombine(args, PP_LW_data);
  distributionVector AP_LW_raw = readCombine(args, AP_LW_data);

  jackknifeTimeSeriesD PP_LW_jack = resampleVector(PP_LW_raw, args.Lt);
  jackknifeTimeSeriesD AP_LW_jack = resampleVector(AP_LW_raw, args.Lt);
  
  publicationPrint printer;

  const int ntypes = 2;
  const DataType type_map[ntypes] = {PP_LW_data, AP_LW_data};
  const distributionVector* data_map_raw[ntypes] = {&PP_LW_raw, &AP_LW_raw};
  const jackknifeTimeSeriesD* data_map_jack[ntypes] = {&PP_LW_jack, &AP_LW_jack};
  const int nx = ntypes * args.Lt; //t + Lt*type

  typedef dataSeries<Coord, doubleJackknifeDistributionD> doubleJackknifeAllData;
  typedef dataSeries<Coord, jackknifeDistributionD> jackknifeAllData;
  
  doubleJackknifeAllData data_dj(nx, ntraj);
  jackknifeAllData data_j(nx, ntraj);
  for(int tt=0;tt<nx;tt++){
    const int t = tt % args.Lt;
    const DataType d = type_map[tt / args.Lt];
    
    const distributionVector &vraw = *data_map_raw[tt / args.Lt];
    if(vraw.size() == 0) continue;
    assert(vraw.size() == args.Lt);    
    data_dj.coord(tt) = data_j.coord(tt) = Coord(t,d);
    data_dj.value(tt).resample(vraw[t]);

    const jackknifeTimeSeriesD &vjack = *data_map_jack[tt / args.Lt];
    assert(vjack.size() == args.Lt);  
    data_j.value(tt) = vjack.value(t);

    printer << "Resampled data " << data_j.coord(tt) << "  " << data_j.value(tt) << std::endl;
  }

  //Setup fit range
  typedef filteredDataSeries<doubleJackknifeAllData> filteredDoubleJackknifeAllData;
  typedef filteredDataSeries<jackknifeAllData> filteredJackknifeAllData;

  filterCoordTrange tfilter(args.t_min,args.t_max);
  
  filteredDoubleJackknifeAllData inrange_data_dj(data_dj,tfilter);
  filteredJackknifeAllData inrange_data_j(data_j,tfilter);
  
  tfilter.invert();

  filteredJackknifeAllData outofrange_data_j(data_j,tfilter);
  
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
  
  typedef sampleSeries<const filteredJackknifeAllData> sampleSeriesConstType; //const access
  typedef UncorrelatedChisqCostFunction<decltype(fitfunc), sampleSeriesConstType, double, NumericVector<double> > CostFunctionType;
  typedef CostFunctionType::CostType CostType;

  //Setup minimizer
  typedef MarquardtLevenbergMinimizer<CostFunctionType> MinimizerType;

  MarquardtLevenbergParameters<CostType> mlparams;
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
  MatPlotLibScriptGenerate plotter;
  for(int type_idx=0;type_idx<ntypes;type_idx++){
    if(data_map_jack[type_idx]->size() == 0) continue;
    typedef DataSeriesAccessor<jackknifeTimeSeriesD, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistributionD> > Accessor;
    typedef MatPlotLibScriptGenerate::handleType Handle;
    jackknifeTimeSeriesD effmass = effectiveMass(*data_map_jack[type_idx], type_map[type_idx], args.Lt);
    Accessor a(effmass);
    Handle ah = plotter.plotData(a);
    std::string nm = toStr(type_map[type_idx]);
    nm = nm.substr(0,nm.size()-5);    
    plotter.setLegend(ah, nm);
  }
  plotter.createLegend();
  plotter.write("mpi_plot.py");

  
  return 0;
}
