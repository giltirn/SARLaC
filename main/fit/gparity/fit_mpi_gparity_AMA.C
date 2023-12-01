//#define BOOST_SPIRIT_X3_DEBUG
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <fstream>
#include <algorithm>
#include <sstream>

#include <random.h>
#include <plot.h>
#include <distribution.h>
#include <data_series.h>
#include <parser.h>
#include <fit.h>
#include <common.h>

using namespace CPSfit;

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

#include <fit_mpi_gparity_AMA/parse_data.h>
#include <fit_mpi_gparity_AMA/args.h>
#include <fit_mpi_gparity_AMA/fitfunc.h>
#include <fit_mpi_gparity_AMA/data_manipulations.h>
#include <fit_mpi_gparity_AMA/read_data.h>
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

typedef dataSeries<Coord, doubleJackknifeDistributionD> doubleJackknifeAllData;
typedef dataSeries<Coord, jackknifeDistributionType> jackknifeAllData;

typedef filteredDataSeries<doubleJackknifeAllData> filteredDoubleJackknifeAllData;
typedef filteredDataSeries<jackknifeAllData> filteredJackknifeAllData;

struct fitTypeDefs{
  typedef jackknifeDistributionType DistributionType;
  typedef filteredJackknifeAllData CorrelationFunctionDistribution;
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

  parse(args,arg_file);

  std::cout << "Read arguments: \n" << args << std::endl;

  const int ntraj = (args.traj_lessthan - args.traj_start)/args.traj_inc;
  assert(ntraj > 0);

  //Read and resample the data
  TypeInfo type_map;
  std::vector<std::string> types;
  std::vector<rawDataDistributionVector> raw;
  std::vector<jackknifeTimeSeriesType> jack;
  for(int i=0;i<args.data.size();i++){
    if(args.data[i].FF_data.include_data || args.data[i].BB_data.include_data){
      raw.push_back( readCombine(args, i) );
      jack.push_back( resampleVector<jackknifeTimeSeriesType>(raw.back(), args.Lt) );
      types.push_back( args.data[i].type );
      type_map.registerType( args.data[i].type, args.data[i].fitfunc );
    }
  }
  const int ndata_types = raw.size();
  
  //Combine data and double-jackknife resample
  const int nx = ndata_types * args.Lt; //t + Lt*type

  doubleJackknifeAllData data_dj(nx, ntraj);
  jackknifeAllData data_j(nx, ntraj);
  for(int tt=0;tt<nx;tt++){
    const int t = tt % args.Lt;
    const int type_idx = tt / args.Lt;    
    const std::string &d = types[type_idx];
    
    const rawDataDistributionVector &vraw = raw[type_idx];
    assert(vraw.size() == args.Lt);    
    data_dj.coord(tt) = data_j.coord(tt) = Coord(t,d,type_map);
    data_dj.value(tt).resample(vraw[t]);

    const jackknifeTimeSeriesType &vjack = jack[type_idx];
    assert(vjack.size() == args.Lt);  
    data_j.value(tt) = vjack.value(t);
  }

  //Setup fit range
  filterCoordTrange tfilter(args.t_min,args.t_max);
  
  filteredDoubleJackknifeAllData inrange_data_dj(data_dj,tfilter);
  filteredJackknifeAllData inrange_data_j(data_j,tfilter);
  
  tfilter.invert();

  filteredJackknifeAllData outofrange_data_j(data_j,tfilter);
  
  const int ndata_fit = inrange_data_j.size();
    
  //Uncorrelated fit
  typedef NumericSquareMatrix<jackknifeDistributionType> jackknifeMatrix;
  jackknifeMatrix cov(ndata_fit,jackknifeDistributionType(ntraj,0.));
  NumericVector<jackknifeDistributionType> sigma(ndata_fit);
  for(int x=0;x<ndata_fit;x++){
#ifdef FIT_MPI_JACKKNIFE_PROPAGATE_CENTRAL
    jackknifeDistributionD tmp = doubleJackknifeDistributionD::covariance(inrange_data_dj.value(x), inrange_data_dj.value(x));
    cov(x,x).import(tmp);
#else
    cov(x,x) = doubleJackknifeDistributionD::covariance(inrange_data_dj.value(x), inrange_data_dj.value(x));
#endif
    sigma(x) = sqrt(cov(x,x));
    std::cout << inrange_data_dj.coord(x) << " cov " << cov(x,x) << " sigma " << sigma(x) << std::endl; 
  }

  std::cout << "Data included in fit:\n";
  for(int i=0;i<ndata_fit;i++){
    std::cout << inrange_data_j.coord(i) << "  " << inrange_data_j.value(i) << " with sigma " << sigma[i] << std::endl;
  }

  FitMpi fitfunc(2*args.Lt, type_map);

  //Load guesses if applicable
  typedef typename FitMpi::ParameterType ParamContainer;
  ParamContainer guess(fitfunc.Nparams());

  if(cmdline.load_guess){
    std::ifstream f(cmdline.guess_file.c_str());
    std::string line; 

    std::string param;
    double v;
    while(std::getline(f,line)){
      std::cout << "Parsing guess file line '" << line << "'\n"; std::cout.flush();
      std::vector<std::string> tok;
      boost::split(tok, line, [&](const char c){ return c == ' '; });
      assert(tok.size() == 2);

      param = tok[0];
      v = strToAny<double>(tok[1]);

      if(param == "Mass") guess[0] = v;
      else{
	if(!type_map.containsType(param)){
	  std::cout << "Warning: Guess provided for fit parameter " << param << " but this parameter is not in the map\n";
	  std::cout << "Known types:\n";
	  for(int i=0;i<type_map.nTypes();i++) std::cout << "\t" << type_map.typeName(i) << std::endl;
	}else{
	  int idx = type_map.typeIdx(param);
	  guess[idx+1] = v;
	}
      }
    }
  }else{
    guess(0) = 0.5;
    for(int i=1;i<guess.size();i++) guess(i) = 1e3;
  }
  
  std::cout << "Using guess: " << guess << std::endl;
  
  
  //Do the fit
  jackknifeDistributionT<ParamContainer> params(ntraj, guess);
  
  typedef composeFitPolicy<FitMpi, standardFitFuncPolicy, uncorrelatedFitPolicy, MarquardtLevenbergMinimizerPolicy, jackknifeDistributionType, fitTypeDefs>::type FitPolicies;

  jackknifeDistributionType chisq;
  jackknifeDistributionType chisq_per_dof;
  
  fitter<FitPolicies> fit;
  fit.importCostFunctionParameters(sigma);
  fit.importFitFunc(fitfunc);
  fit.fit(params, chisq, chisq_per_dof, inrange_data_j);

  std::cout << "Fit results:\n";
  std::cout << "Mass = " << ParameterPrint<decltype(params)>(params,0) << std::endl;
  for(int i=0;i<type_map.nTypes();i++){
    std::cout << type_map.typeName(i) << " = " << ParameterPrint<decltype(params)>(params,i+1) << std::endl;
  }

  std::cout << "chi^2 = " << chisq << std::endl;
  std::cout << "chi^2/dof = " << chisq_per_dof << std::endl;

  //Plot some nice things  
  MatPlotLibScriptGenerate plotter;
  std::vector<std::string> pallete = {
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
    std::cout << "Computing effective mass for type " << types[type_idx] << std::endl;
    typedef DataSeriesAccessor<jackknifeTimeSeriesType, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistributionType> > Accessor;

    jackknifeTimeSeriesType effmass = effectiveMass(jack[type_idx], types[type_idx], args.Lt, type_map);
    Accessor a(effmass);
    if(type_idx < pallete.size()) plot_args["color"] = pallete[type_idx];
    Handle ah = plotter.plotData(a,plot_args);
    plotter.setLegend(ah, types[type_idx]);
  }
  //   Plot the fitted mass as constant
  {
    ParamContainer mn = params.best();
    ParamContainer err = params.standardError();
    const double m = mn(0);
    const double dm = err(0);
    
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

  const double ymid = params.best()(0);
  const double yw = params.standardError()(0) * 20;
  
  plotter.setYaxisBounds(ymid-yw, ymid+yw);
  
  plotter.write("mpi_plot.py");

  std::cout << "Writing results" << std::endl;
#ifdef HAVE_HDF5
  writeParamsStandard(chisq, "chisq.hdf5");
  writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");
  writeParamsStandard(params, "params.hdf5");  
#endif
  
  std::cout << "Done" << std::endl;
  return 0;
}
