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

#define FIT_MPI_JACKKNIFE_PROPAGATE_CENTRAL

#ifdef FIT_MPI_JACKKNIFE_PROPAGATE_CENTRAL
constexpr int first_sample = -1;
template<typename T>
using jackknifeDistributionT = jackknifeCdistribution<T>;
typedef jackknifeCdistributionD jackknifeDistributionType;
typedef jackknifeCtimeSeriesD jackknifeTimeSeriesType;
#else
template<typename T>
using jackknifeDistributionT = jackknifeDistribution<T>;
constexpr int first_sample = 0;
typedef jackknifeDistribution jackknifeDistributionType;
typedef jackknifeTimeSeriesD jackknifeTimeSeriesType;
#endif

#include <fit_mpi_gparity_AMA/args.h>
#include <fit_mpi_gparity_AMA/fitfunc.h>
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
template<typename DistributionType>
class ParameterPrint: public OstreamHook{
  int idx;
  const DistributionType &d;
public:
  ParameterPrint(const DistributionType &_d, const int _idx): d(_d), idx(_idx){}  
  void write(std::ostream &os) const{
    os << "(" << printStats<DistributionType>::centralValue(d)(idx) << " +- " << printStats<DistributionType>::error(d)(idx) << ")";
  }
};



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


  AllParamMap param_map(args);
  
  const int ntraj = (args.traj_lessthan - args.traj_start)/args.traj_inc;
  assert(ntraj > 0);

  //Read and resample the data
  std::vector<DataType> types;
  std::vector<rawDataDistributionVector> raw;
  std::vector<jackknifeTimeSeriesType> jack;
  for(int i=0;i<args.data.size();i++){
    if(args.data[i].FF_data.include_data || args.data[i].BB_data.include_data){
      raw.push_back( readCombine(args, i) );
      jack.push_back( resampleVector(raw.back(), args.Lt) );
      types.push_back( args.data[i].type );
    }
  }
  const int ndata_types = raw.size();
  
  //Combine data and double-jackknife resample
  const int nx = ndata_types * args.Lt; //t + Lt*type

  typedef dataSeries<Coord, doubleJackknifeDistributionD> doubleJackknifeAllData;
  typedef dataSeries<Coord, jackknifeDistributionType> jackknifeAllData;
  
  doubleJackknifeAllData data_dj(nx, ntraj);
  jackknifeAllData data_j(nx, ntraj);
  for(int tt=0;tt<nx;tt++){
    const int t = tt % args.Lt;
    const int type_idx = tt / args.Lt;    
    const DataType d = types[type_idx];
    
    const rawDataDistributionVector &vraw = raw[type_idx];
    assert(vraw.size() == args.Lt);    
    data_dj.coord(tt) = data_j.coord(tt) = Coord(t,d);
    data_dj.value(tt).resample(vraw[t]);

    const jackknifeTimeSeriesType &vjack = jack[type_idx];
    assert(vjack.size() == args.Lt);  
    data_j.value(tt) = vjack.value(t);
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
  typedef NumericSquareMatrix<jackknifeDistributionType> jackknifeMatrix;
  jackknifeMatrix cov(ndata_fit,jackknifeDistributionType(ntraj,0.));
  std::vector<jackknifeDistributionType> sigma(ndata_fit);
  NumericSquareMatrix<double> corr(ndata_fit,0.);
  for(int x=0;x<ndata_fit;x++){
#ifdef FIT_MPI_JACKKNIFE_PROPAGATE_CENTRAL
    jackknifeDistributionD tmp = doubleJackknifeDistributionD::covariance(inrange_data_dj.value(x), inrange_data_dj.value(x));
    cov(x,x).import(tmp);
#else
    cov(x,x) = doubleJackknifeDistributionD::covariance(inrange_data_dj.value(x), inrange_data_dj.value(x));
#endif
    cov(x,x) = cov(x,x) * double(ntraj-1); //properly normalize

    sigma[x] = sqrt(cov(x,x));
    corr(x,x) = 1.;
    std::cout << inrange_data_dj.coord(x) << " cov " << cov(x,x) << " sigma " << sigma[x] << std::endl; 
  }

  std::cout << "Data included in fit:\n";
  for(int i=0;i<ndata_fit;i++){
    std::cout << inrange_data_j.coord(i) << "  " << inrange_data_j.value(i) << " with sigma " << sigma[i] << std::endl;
  }


  //Load guesses if applicable
  typedef typename FitMpi::ParameterType ParamContainer;
  ParamContainer guess(param_map);

  if(cmdline.load_guess){
    std::ifstream f(cmdline.guess_file.c_str());
    assert(f.good());
    while(!f.eof()){
      Params p; double v;
      f >> p >> v;
      assert(!f.bad() && !f.fail());
      guess(p) = v;
    }
    f.close();
  }else{
    for(int i=1;i<guess.size();i++) guess(i) = 1e3;
    guess(Mass) = 0.5;
  }
  
  std::cout << "Using guess: " << guess << std::endl;
  
  
  //Do the fit
  FitMpi fitfunc(2*args.Lt, param_map);
  
  typedef sampleSeries<const filteredJackknifeAllData> sampleSeriesConstType; //const access
  typedef UncorrelatedChisqCostFunction<decltype(fitfunc), sampleSeriesConstType, double, NumericVector<double> > CostFunctionType;
  typedef CostFunctionType::CostType CostType;

  //Setup minimizer
  typedef MarquardtLevenbergMinimizer<CostFunctionType> MinimizerType;

  MarquardtLevenbergParameters<CostType> mlparams;
  jackknifeDistributionT<ParamContainer> params(ntraj, guess);

  jackknifeDistributionT<CostType> chisq(ntraj);
  jackknifeDistributionT<CostType> chisqperdof(ntraj);
  int dof = ndata_fit - fitfunc.Nparams();

#pragma omp parallel for
  for(int j=first_sample;j<ntraj;j++){    
    sampleSeriesConstType dsample(inrange_data_j, j);
    std::vector<double> sigma_j(sigma.size()); for(int p=0;p<sigma.size();p++) sigma_j[p] = sigma[p].sample(j);
    CostFunctionType costfunc(fitfunc, dsample, sigma_j);

    MinimizerType fitter(costfunc, mlparams);
    
    ParamContainer &pj = params.sample(j);
    CostType cost = fitter.fit(pj);
    assert(fitter.hasConverged());
    chisq.sample(j) = cost;
    chisqperdof.sample(j) = cost/dof;
  }
  std::cout << "chi^2 = " << chisq << "\nchi^2/dof = " << chisqperdof << " (" << dof << " degrees of freedom)\n";

  std::cout << "Fit results:\n";
  for(int i=0;i<param_map.nParams();i++){
    std::cout << param_map.paramName(i) << " = " << ParameterPrint<decltype(params)>(params,i) << std::endl;
  }

  //Plot some nice things  
  MatPlotLibScriptGenerate plotter;
  std::string pallete[] = {
    "#ff0000",
    "#00ff00",
    "#0000ff",
    "#ff9c00",
    "#a3007f",
    "#00a000",
    "#ff00ff",
    "#00ffff",
    "#d38d4e"};

  typename MatPlotLibScriptGenerate::kwargsType plot_args;
  typedef MatPlotLibScriptGenerate::handleType Handle;
  
  for(int type_idx=0;type_idx<ndata_types;type_idx++){
    typedef DataSeriesAccessor<jackknifeTimeSeriesType, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistributionType> > Accessor;

    jackknifeTimeSeriesType effmass = effectiveMass(jack[type_idx], types[type_idx], args.Lt);
    Accessor a(effmass);
    plot_args["color"] = pallete[type_idx];
    Handle ah = plotter.plotData(a,plot_args);
    std::string nm = toString(types[type_idx]);
    nm = nm.substr(0,nm.size()-5);    
    plotter.setLegend(ah, nm);
  }
  //   Plot the fitted mass as constant
  {
    ParamContainer mn = params.best();
    ParamContainer err = params.standardError();
    const double m = mn(Mass);
    const double dm = err(Mass);
    
    std::vector<double> x = {double(args.t_min), double(args.t_max)};
    std::vector<double> upper = {m + dm, m + dm};
    std::vector<double> lower = {m - dm, m - dm};    
    BandVectorAccessor band(x,upper,lower);
    plot_args["color"] = "r";
    plot_args["alpha"] = 0.2;
    Handle ah = plotter.errorBand(band, plot_args);
    plotter.setLegend(ah,"Fit");
  }
  
  plotter.createLegend();
  plotter.setXlabel("$t$");
  plotter.setYlabel("$m_{\\rm eff}(t)$");
  plotter.setXaxisBounds(-0.2,args.Lt+0.2);

  const double ymid = params.best()(Mass);
  const double yw = params.standardError()(Mass) * 20;
  
  plotter.setYaxisBounds(ymid-yw, ymid+yw);
  
  plotter.write("mpi_plot.py");

#ifdef HAVE_HDF5
  writeParamsStandard(params, "params.hdf5");  
#endif
  
  return 0;
}
