//If we perform an N-state/N-operator fit we can compute the coefficients for an operators that optimally projects onto the ground-sate
//This can then be re-fit with a one-parameter fit to the energy for a single correlator

#include <utils.h>
#include <random.h>
#include <fit/bootstrap_pvalue.h>
#include <pipi_common/pipi_common.h>
using namespace SARLaC;

#include<fit_pipi_gnd_exc_sim_gparity/enums.h>
#include<fit_pipi_gnd_exc_sim_gparity/fit.h>

#include<fit_pipi_gnd_exc_sigma_sim_gparity/fitfunc.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/filters.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/raw_data.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/resampled_data.h>

#include<fit_pipi_optimal_bootstrap/args.h>
#include<fit_pipi_optimal_bootstrap/cmdline.h>

template<typename CorrFuncType>
CorrFuncType computeProjectedCorrelator(const ResampledData<CorrFuncType> &data, const std::vector<Operator> &ops, const int Lt, const std::vector<bootstrapDistributionD> &r){
  int nop = ops.size();
  auto zero = data.correlator(ops[0],ops[0]).value(0);
  zeroit(zero);

  CorrFuncType Copt(Lt);

  for(int t=0;t<Lt;t++){
    Copt.coord(t) = t;
    auto &v = Copt.value(t);
    v = zero;

    for(int a=0;a<nop;a++)
      for(int b=0;b<nop;b++){
	int aa=a, bb=b;
	if(!data.haveData(ops[a],ops[b]) && data.haveData(ops[b],ops[a])){ //only triangular part of matrix is saved
	  aa=b; bb=a;
	}
	assert(data.haveData(ops[aa],ops[bb]));
	       
	auto c = data.correlator(ops[aa],ops[bb]).value(t);

	c.best() = c.best() * r[a].best() * r[b].best();
	for(int s=0;s<c.size();s++) c.sample(s) = c.sample(s) * r[a].sample(s) * r[b].sample(s);
	
	v = v + c;

	//v = v + data_b.correlator(ops[aa],ops[bb]).value(t) * r[a] * r[b];
      }
  }
  return Copt;
}

template<typename CorrFuncType>
CorrFuncType extractFitData(const CorrFuncType &data, const int t_min, const int t_max){
  CorrFuncType Copt_inrange(t_max - t_min + 1);
  int i=0;
  for(int t=t_min; t<=t_max; t++){
    Copt_inrange.coord(i) = t;
    Copt_inrange.value(i) = data.value(t);
    i++;
  }
  return Copt_inrange;
}

class FitMultiExp{
  int nstate;
public:
  typedef parameterVector<double> Params;
  typedef double ValueType;
  typedef Params ParameterType;
  typedef Params ValueDerivativeType;
  typedef double GeneralizedCoordinate;

  template<typename T, typename Boost>
  static inline T boostit(const int idx, const ParameterType &p, const Boost &b){
    return b(p(idx),idx);
  }

  template<typename T, typename Boost>
  T eval(const GeneralizedCoordinate &x, const ParameterType &p, const Boost &b) const{  
    T out(0.);
    
    for(int s=0;s<nstate;s++){
      auto A = boostit<T,Boost>(2*s, p,b);
      auto E = boostit<T,Boost>(2*s+1, p,b);
      out = out + A*exp(-E*x);
    }
    return out;
  }

  inline ValueType value(const GeneralizedCoordinate &x, const ParameterType &p) const{
    return eval<double>(x,p,[&](const double a, const int i){ return a; });
  }
  inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &x, const ParameterType &p) const{
    ValueDerivativeType d = p;
    for(int i=0;i<2*nstate;i++) d(i) = eval<dual>(x,p,[&](const double v, const int j){ return dual(v, j==i ? 1.:0.); }).xp;
    return d;
  }

  FitMultiExp(const int nstate): nstate(nstate){}

  inline int Nparams() const{ return 2*nstate; }

  inline int Nstate() const{ return nstate; }
};





int main(const int argc, const char* argv[]){
  RNG.initialize(1234);

  Args args;
  if(argc < 2){
    std::ofstream of("template.args");
    of << args;
    std::cout << "Wrote template argument file to template.args\n";
    return 0;
  }
  {
    const std::string arg_file = argv[1];
    parse(args,arg_file);
    std::cout << "Read arguments: \n" << args << std::endl;
  }
  CMDline cmdline(argc,argv,2);

  if(cmdline.save_guess_template){
    parameterVector<double> v(2*args.nstate_fit);
    std::ofstream of("guess_template.dat");
    of << v;
    of.close();
    exit(0);
  }

  const std::vector<Operator> &ops = args.operators;  
  int nop = ops.size();

  assert(args.op_amplitudes.size() == nop);
  for(int i=0;i<nop;i++) assert(args.op_amplitudes[i].idx.size() == nop); //coupling of operators to states

  std::vector<bootstrapDistributionD> r = computeR(args.op_amplitudes, args.Ascale);

  //Read data
  RawData raw_data;
  if(!cmdline.load_combined_data){
    if(cmdline.load_raw_data){
      std::cout << "Reading raw data from " << cmdline.load_raw_data_file << std::endl;
      HDF5reader rd(cmdline.load_raw_data_file);  raw_data.read(rd, "raw_data");
      for(int i=0;i<ops.size();i++)
	for(int j=i;j<ops.size();j++)
	  assert(raw_data.haveData(ops[i],ops[j]));
    }else{
      raw_data.read(args.Lt, args.data_dir, args.traj_start, args.traj_inc, args.traj_lessthan,
		    args.pipi_figure_file_format, args.pipi_bubble_file_format, args.tsep_pipi, args.tstep_pipi,
		    args.pipi_to_sigma_file_format, args.tstep_pipi_to_sigma,
		    args.sigma2pt_file_format, args.sigma_bubble_file_format,
		    ops, 0);
    }
    if(cmdline.save_raw_data){
      std::cout << "Saving raw data to " << cmdline.save_raw_data_file << std::endl;
      HDF5writer wr(cmdline.save_raw_data_file);  raw_data.write(wr, "raw_data");
    }
  }

  ResampledData<bootstrapCorrelationFunctionD> data_b;
  ResampledData<bootJackknifeCorrelationFunctionD> data_bj;
  if(cmdline.load_combined_data){
    HDF5reader rd(cmdline.load_combined_data_file);
    read(rd, data_b, "data_b");
    read(rd, data_bj, "data_bj");
    for(int i=0;i<ops.size();i++)
      for(int j=i;j<ops.size();j++)
	if(!data_b.haveData(ops[i],ops[j]))
	  error_exit(std::cout << "Loaded checkpoint does not contain data for (" << ops[i] << ", " << ops[j] << ")\n");
  }else{
    int nsample = raw_data.correlator(ops[0],ops[0]).value(0).size();

    //Generate the resample table
    resampleTableOptions ropt;
    ropt.read_from_file = cmdline.load_boot_resample_table;
    ropt.read_file = cmdline.load_boot_resample_table_file;
    ropt.write_to_file = cmdline.save_boot_resample_table;
    ropt.write_file = cmdline.save_boot_resample_table_file;
    
    std::vector<std::vector<int> > rtable = generateResampleTable(nsample, args.nboot, args.resample_table_type, args.block_size, RNG, ropt);
    if(rtable.size() != args.nboot) error_exit(std::cout << "Expected resample table of size " << args.nboot << ", got " << rtable.size() << std::endl);

    bootstrapBlockResampler resampler(rtable);
    data_b.generatedResampledData(raw_data, resampler, args.Lt, args.tsep_pipi, args.do_vacuum_subtraction, args.timeslice_avg_vac_sub);
    data_bj.generatedResampledData(raw_data, resampler, args.Lt, args.tsep_pipi, args.do_vacuum_subtraction, args.timeslice_avg_vac_sub);
  }  

  if(cmdline.save_combined_data){
    HDF5writer wr(cmdline.save_combined_data_file);
    write(wr, data_b, "data_b");
    write(wr, data_bj, "data_bj");
  }

  bootstrapDistributionD zero(data_b.correlator(ops[0],ops[0]).value(0)); zeroit(zero);
  int nboot = zero.size();

  //Compute projected correlator
  bootstrapCorrelationFunctionD Copt_b = computeProjectedCorrelator(data_b, ops, args.Lt, r);
  bootJackknifeCorrelationFunctionD Copt_bj = computeProjectedCorrelator(data_bj, ops, args.Lt, r);

  //Filter by fit range
  bootstrapCorrelationFunctionD Copt_inrange_b = extractFitData(Copt_b, args.t_min, args.t_max);
  bootJackknifeCorrelationFunctionD Copt_inrange_bj = extractFitData(Copt_bj, args.t_min, args.t_max);

  //Due to the pi-pi separation in the sources, the correlators don't have exactly the same time dependence. However if we restrict ourselves to t_max << Lt we can ignore the backwards propagating state

  typedef FitMultiExp FitFunc;
  typedef parameterVector<double> Params;

  FitFunc fitfunc(args.nstate_fit);
  Params guess;

  if(cmdline.load_guess)  {
    parse(guess, cmdline.guess_file);
    std::cout << "Loaded guess: " << guess << std::endl;
  }else{
    guess.resize(2*args.nstate_fit+1);
    for(int i=0;i<args.nstate_fit;i++){
      guess(2*i) = 1.;
      guess(2*i+1) = (1+i)*0.35;
    }
  }

  genericFitFuncWrapper<FitFunc> fwrap(fitfunc,guess);
 
  simpleFitWrapper<bootstrapDistributionD> fitter(fwrap, MinimizerType::MarquardtLevenberg);

  if(cmdline.load_minimizer_params) fitter.setMinimizer(MinimizerType::MarquardtLevenberg, 
							getMinimizerParams(MinimizerType::MarquardtLevenberg, cmdline.load_minimizer_params, cmdline.minimizer_params_file));

  if(cmdline.load_frozen_fit_params) readFrozenParams<bootstrapDistributionD>(fitter, cmdline.load_frozen_fit_params_file, nboot);

  fitter.generateCovarianceMatrix(Copt_inrange_bj, args.correlated ? CostType::Correlated : CostType::Uncorrelated);

  bootstrapDistribution<Params> params(zero.getInitializer(),guess);
  bootstrapDistributionD chisq(zero), chisq_per_dof(zero);
  int dof;

  fitter.fit(params, chisq, chisq_per_dof, dof, Copt_inrange_b);
  
  std::cout << "Params:\n";
  {
    for(int p=0;p<params.sample(0).size();p++){
      int s = p/2;

      bootstrapDistributionD tmp;
      standardIOhelper<bootstrapDistributionD, bootstrapDistribution<Params> >::extractStructEntry(tmp, params, p);
      std::cout << (p%2 == 0 ? 'A' : 'E') << s << " = " << tmp << std::endl;
    }
  }

  std::cout << std::endl;
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;

#ifdef HAVE_HDF5
  writeParamsStandard(params, "params.hdf5");
  writeParamsStandard(chisq, "chisq.hdf5");
  writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");
#endif

  //Compute bootstrap p-value
  {
    //Recenter
    bootstrapCorrelationFunctionD Copt_inrange_b_rc(Copt_inrange_b);
    for(int tt=0;tt<Copt_inrange_b_rc.size();tt++){
      double t = Copt_inrange_b_rc.coord(tt);
      double fval = fitfunc.value(t, params.best());
      Copt_inrange_b_rc.value(tt) = Copt_inrange_b_rc.value(tt) + (fval - Copt_inrange_b.value(tt).best());
    }
    bootstrapDistribution<Params> params_rc(zero.getInitializer(),guess);
    bootstrapDistributionD chisq_rc(zero), chisq_per_dof_rc(zero);
  
    fitter.fit(params_rc, chisq_rc, chisq_per_dof_rc, dof, Copt_inrange_b_rc);

    //Compute p-value

    //Note we shouldn't use the "best" value here because it is not obtained using a bootstrap resampling; instead use the mean
    std::cout << "Bootstrap distribution of q^2 has mean " << chisq_rc.mean() << " and std.dev " << chisq_rc.standardDeviation() << std::endl;

    std::vector<double> q2_dist(nboot); for(int i=0;i<nboot;i++) q2_dist[i] = chisq_rc.sample(i);
    std::sort(q2_dist.begin(), q2_dist.end(), [&](const double a, const double b){ return a<b; });

    double p_boot = computePvalue(chisq.best(), q2_dist);
  
    {
      std::ofstream of("p_boot.dat");
      of << p_boot << std::endl;
    }

    //Save sorted bootstrap q2 distribution as rawDataDistribution
    rawDataDistributionD q2_dist_r(nboot, [&](const int i){ return q2_dist[i]; });
    writeParamsStandard(q2_dist_r, "boot_q2.hdf5");
  
    std::cout << "Bootstrap p-value " << p_boot << std::endl;
  }
  
  //Plot single-exponential effective mass
  {
    const int Eidx = 1;
    StandardFitParams base(1., 0.5);
    FitExp fitfunc_eff;

    correlationFunction<double, bootstrapDistributionD> Ceff = effectiveMass2pt(Copt_b, fitfunc_eff, base, 1, args.Lt, 0.4);

    MatPlotLibScriptGenerate plot;
  
    typedef DataSeriesAccessor< correlationFunction<double, bootstrapDistributionD>, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<bootstrapDistributionD> > acc;
  
    plot.plotData(acc(Ceff), "E_eff");
  
    bootstrapDistributionD tmp;
    standardIOhelper<bootstrapDistributionD, bootstrapDistribution<Params> >::extractStructEntry(tmp, params, Eidx);

    typedef BandRangeConstantDistributionValue<bootstrapDistributionD, DistributionPlotAccessor<bootstrapDistributionD> > bacc;
    MatPlotLibScriptGenerateBase::kwargsType kwargs;
    kwargs["alpha"] = 0.6;

    plot.errorBand(bacc(args.t_min, args.t_max, tmp), kwargs, "Fit");

    plot.write("plot.py", "plot.pdf");
  }
  
  std::cout << "Done\n";
  return 0;
}
