#include<common.h>
#include<fit.h>
#include<parser.h>
#include<plot.h>

using namespace CPSfit;

#include <fit_mpi_gparity_AMA/parse_data.h>
#include <fit_mpi_gparity_AMA/data_manipulations.h>


#define BK_GPARITY_AMA_ARGS_MEMBERS		\
  ( std::string, threepoint_fmt_sloppy )		\
  ( std::string, threepoint_fmt_exact )		\
  ( std::string, AT_P_F_fmt_sloppy )		\
  ( std::string, AT_P_F_fmt_exact )		\
  ( std::string, AT_P_B_fmt_sloppy )		\
  ( std::string, AT_P_B_fmt_exact )		\
  ( int, L)							\
  ( int, Lt)							\
  ( int, t_min)							\
  ( int, t_max)							\
  ( int, traj_start )						\
  ( int, traj_inc )						\
  ( int, traj_lessthan ) 

struct BKgparityAMAargs{
  GENERATE_MEMBERS(BK_GPARITY_AMA_ARGS_MEMBERS);

  BKgparityAMAargs(): L(16), Lt(32), t_min(6), t_max(26), traj_start(0), traj_inc(40), traj_lessthan(1000){}
} ;
GENERATE_PARSER(BKgparityAMAargs, BK_GPARITY_AMA_ARGS_MEMBERS);

void process(jackknifeTimeSeriesD &out_j, doubleJackknifeTimeSeriesD &out_dj, const std::string &sloppy_fmt, const std::string &exact_fmt, const int traj_start, const int traj_inc, const int traj_lessthan, const int Lt, const std::string &descr){
  int nsample = (traj_lessthan - traj_start) / traj_inc;
  rawDataDistributionMatrix sloppy_raw(32, rawDataDistribution<double>(nsample));
  rawDataDistributionMatrix exact_raw(32, rawDataDistribution<double>(nsample));

  for(int s=0;s<nsample;s++){
    int traj = traj_start + traj_inc * s;

    read(sloppy_raw, s, sloppy_fmt, traj, Real);
    read(exact_raw, s, exact_fmt, traj, Real);
  }
  rawDataDistributionVector sloppy_avg = sourceTimeSliceAverage(sloppy_raw);
  rawDataDistributionVector correction = computeAMAcorrection(sloppy_raw,exact_raw);
  
  jackknifeTimeSeriesD sloppy_avg_j = resampleVector<jackknifeTimeSeriesD>(sloppy_avg,Lt);
  jackknifeTimeSeriesD correction_j = resampleVector<jackknifeTimeSeriesD>(correction,Lt);

  out_j = jackknifeTimeSeriesD(Lt, [&](const int t){ return jackknifeTimeSeriesD::ElementType(t, sloppy_avg_j.value(t) + correction_j.value(t)); });
  
  std::cout << descr << ": t sloppy correction corrected\n";
  for(int t=0;t<Lt;t++){
    std::cout << t << " " << sloppy_avg_j.value(t) << " " << correction_j.value(t) << " " << out_j.value(t) << std::endl;
  }

  doubleJackknifeTimeSeriesD sloppy_avg_dj = resampleVector<doubleJackknifeTimeSeriesD>(sloppy_avg,Lt);
  doubleJackknifeTimeSeriesD correction_dj = resampleVector<doubleJackknifeTimeSeriesD>(correction,Lt);

  out_dj = doubleJackknifeTimeSeriesD(Lt, [&](const int t){ return doubleJackknifeTimeSeriesD::ElementType(t, sloppy_avg_dj.value(t) + correction_dj.value(t)); });
}


int main(const int argc, const char* argv[]){
  BKgparityAMAargs args;
  if(argc < 2){
    std::ofstream of("template.args");
    (std::cout << "No parameter file provided: writing template to 'template.args' and exiting\n").flush();
    of << args;
    return 1;
  }    
 
  parse(args, argv[1]);

  jackknifeTimeSeriesD num_j, den_F_j, den_B_j;
  doubleJackknifeTimeSeriesD num_dj, den_F_dj, den_B_dj;

  process(num_j, num_dj, args.threepoint_fmt_sloppy, args.threepoint_fmt_exact, args.traj_start, args.traj_inc, args.traj_lessthan, args.Lt, "Numerator");
  process(den_F_j, den_F_dj, args.AT_P_F_fmt_sloppy, args.AT_P_F_fmt_exact, args.traj_start, args.traj_inc, args.traj_lessthan, args.Lt, "Den F");
  process(den_B_j, den_B_dj, args.AT_P_B_fmt_sloppy, args.AT_P_B_fmt_exact, args.traj_start, args.traj_inc, args.traj_lessthan, args.Lt, "Den B");

  int V = args.L * args.L * args.L;

  jackknifeCorrelationFunctionD BK_j(args.Lt, [&](const int t){ return jackknifeCorrelationFunctionD::ElementType(t, V * 3./8 * num_j.value(t)/den_F_j.value(t)/den_B_j.value(t)); });
  doubleJackknifeCorrelationFunctionD BK_dj(args.Lt, [&](const int t){ return doubleJackknifeCorrelationFunctionD::ElementType(t, V * 3./8 * num_dj.value(t)/den_F_dj.value(t)/den_B_dj.value(t)); });
  
  std::cout << "BK: t value\n";
  for(int t=0;t<args.Lt;t++){
    std::cout << t << " " << BK_j.value(t) << std::endl;
  }

  int nsample = BK_j.value(0).size();
  int nfit = args.t_max-args.t_min+1;
  
  jackknifeCorrelationFunctionD BK_j_inrange(nfit,[&](const int i){ return BK_j[args.t_min+i]; });
  doubleJackknifeCorrelationFunctionD BK_dj_inrange(nfit,[&](const int i){ return BK_dj[args.t_min+i]; });

  typedef FitConstant FitFunc;
  typedef composeFitPolicy<FitFunc,standardFitFuncPolicy,uncorrelatedFitPolicy>::type FitPolicies;
  typedef importCostFunctionParameters<uncorrelatedFitPolicy,FitPolicies> Importer;
  
  FitFunc fitfunc;

  fitter<FitPolicies> fit;
  fit.importFitFunc(fitfunc);
  
  Importer import(fit, BK_dj_inrange);

  jackknifeDistributionD chisq, chisq_per_dof;
  FitPolicies::FitParameterDistribution params(nsample,fitfunc.guess());
  fit.fit(params, chisq, chisq_per_dof, BK_j_inrange);

  typedef standardIOhelper<jackknifeDistributionD, FitPolicies::FitParameterDistribution> Unwrapper;
  std::vector<jackknifeDistributionD> params_unwrap(fitfunc.Nparams());

  std::cout << "Params:\n";
  for(int i=0;i<fitfunc.Nparams();i++){
    Unwrapper::extractStructEntry(params_unwrap[i], params, i);
    std::cout << i << " " << params_unwrap[i] << std::endl;
  }
  
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
  
  //Plot result
  MatPlotLibScriptGenerate plotter;
  typedef MatPlotLibScriptGenerate::handleType Handle;
  typename MatPlotLibScriptGenerate::kwargsType plot_args;
    
  typedef DataSeriesAccessor<jackknifeCorrelationFunctionD, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistributionD> > Accessor;
  Accessor a(BK_j);
  Handle ah = plotter.plotData(a);

  typename FitFunc::ParameterType mn = params.best();
  typename FitFunc::ParameterType err = params.standardError();
  const double m = mn(0);
  const double dm = err(0);
  
  std::vector<double> x = {double(args.t_min), double(args.t_max)};
  std::vector<double> upper = {m + dm, m + dm};
  std::vector<double> lower = {m - dm, m - dm};    
  BandVectorAccessor band(x,upper,lower);
  plot_args["alpha"] = 0.2;
  ah = plotter.errorBand(band, plot_args);
   
  plotter.createLegend();
  plotter.setXlabel("$t$");
  plotter.setYlabel("$B_{K}^{ {\\rm latt} }(t)$");
  plotter.setXaxisBounds(-0.2,args.Lt+0.2);

  const double ymid = m;
  const double yw = 20 * dm;
  
  plotter.setYaxisBounds(ymid-yw, ymid+yw);

  std::cout << "Writing plot to 'bk.py'\n";  
  plotter.write("bk.py", "bk.pdf");

  std::cout << "Finished\n";
  return 0;
}
