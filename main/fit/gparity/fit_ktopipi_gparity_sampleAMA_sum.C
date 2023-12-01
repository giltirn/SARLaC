#if 1

//Code needs to be updated and fixed

int main(void){
  return 0;
}


#else

#include <ktopipi_common/ktopipi_common.h>
#include <ktopipi_sampleAMA_common/ktopipi_sampleAMA_common.h>

using namespace CPSfit;

//#define PRINT_CORRECTION

#include <fit_ktopipi_gparity_sampleAMA/cmdline.h>
#include <fit_ktopipi_gparity_sampleAMA/args.h>
#include <fit_ktopipi_gparity_sampleAMA/fit_sama_expand.h>

struct allInputs{
  SampleAMAargs args;
  SampleAMAcmdLine cmdline;

  sampleAMA_resamplers resamplers;

  allInputs(const SampleAMAargs &args, const SampleAMAcmdLine &cmdline): args(args),cmdline(cmdline),
									 resamplers(args.traj_start_S, args.traj_lessthan_S, args.traj_start_C, args.traj_lessthan_C, args.traj_inc, args.bin_size){}
};

//Params   M0 M1 cK mK cpipi Epipi cexc Eexc
class FitSum{

public:
  typedef parameterVector<double> Params;
  typedef double ValueType;
  typedef Params ParameterType;
  typedef Params ValueDerivativeType;
  typedef double GeneralizedCoordinate;

  double a;
  double b;
  
  FitSum(double a, double b): a(a),b(b){}

  template<typename T, typename Boost>
  T eval(const GeneralizedCoordinate t, const ParameterType &p, const Boost &boost) const{
    parameterVector<T> bp(p.size());
    for(int i=0;i<p.size();i++) bp(i) = boost(p(i),i);
    
    T M0 = bp(0), M1 = bp(1), cK = bp(2), mK = bp(3), cpipi = bp(4), Epipi = bp(5), cexc = bp(6), Eexc = bp(7);

    //w(t') = \sum_{ij} M_{ij}\frac{  \left( e^{-E_it'}e^{-(m_j-E_i) a } - e^{-m_j(t'-b+1)}e^{-E_i(b-1)} \right)  }{  \left(1 - e^{-(m_j-E_i)}\right)  }

    //Include only 0,0 and 1,0 terms

    T fac00_num_1 = exp(-Epipi*t)*exp(-(mK-Epipi)*a);
    T fac00_num_2 = exp(-mK*(t-b+1)) * exp(-Epipi*(b-1));
    T fac00_den = 1. - exp(-(mK-Epipi));
    
    T fac10_num_1 = exp(-Eexc*t)*exp(-(mK-Eexc)*a);
    T fac10_num_2 = exp(-mK*(t-b+1)) * exp(-Eexc*(b-1));
    T fac10_den = 1. - exp(-(mK-Eexc));

    return 
      cK * cpipi * M0 * (fac00_num_1 - fac00_num_2)/fac00_den
      +
      cK * cexc * M1 * (fac10_num_1 - fac10_num_2)/fac10_den;
  }
  inline ValueType value(const GeneralizedCoordinate &x, const ParameterType &p) const{
    return eval<double>(x,p,[&](const double v, const int i){ return v; });
  }
  inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &x, const ParameterType &p) const{
    ValueDerivativeType d(8);
    for(int i=0;i<8;i++) d(i) = eval<dual>(x,p,[&](const double v, const int j){ return dual(v, j==i ? 1.:0.); }).xp;
    return d;
  }

  FitSum(){}

  static inline int Nparams(){ return 8; }
};

template<typename DistributionType>
correlationFunction<double, DistributionType> getw(const std::vector<int> &tseps, const std::map<int, std::vector<int> > &data_map,
						   const std::vector<correlationFunction<amplitudeDataCoord, DistributionType> > &A0_all,
						   const int q, const int a, const int b){
  DistributionType zero = A0_all[q].value(0); zeroit(zero);
  int nsample = zero.size();

  enum { isJack = std::is_same<DistributionType,jackknifeDistributionD>::value };

  correlationFunction<double, DistributionType> summed_corr;

  for(auto tsep = tseps.begin(); tsep != tseps.end(); tsep++){
    if(isJack) std::cout << "tsep = " << *tsep << std::endl;
    DistributionType into = zero;

    const std::vector<int> &elems_tsep = data_map.find(*tsep)->second;
    for(int i=0;i<elems_tsep.size();i++){
      double top = A0_all[q].coord(elems_tsep[i]).t;
      if( (int)top < a || (int)top > (*tsep - b) ) continue;
      
      const DistributionType &val = A0_all[q].value(elems_tsep[i]);
      if(isJack) std::cout << "Adding top=" << top << " val=" << val << std::endl;
      into = into + val;
    }        
    summed_corr.push_back(*tsep, into);
  }//tsep loop

  return summed_corr;
}


jackknifeDistributionD groundTimeDependence(const jackknifeDistributionD &cK,
					    const jackknifeDistributionD &mK,
					    const jackknifeDistributionD &cpipi,
					    const jackknifeDistributionD &Epipi,
					    double a, double b, double t){
  jackknifeDistributionD num_1 = exp(-Epipi * t) * exp(- (mK-Epipi)*a);
  jackknifeDistributionD num_2 = exp(-mK *( t - b + 1))*exp(-Epipi*(b-1));
  jackknifeDistributionD den = 1 - exp(-(mK-Epipi) );
  return jackknifeDistributionD(cK*cpipi*(num_1 - num_2)/den);
}

int main(const int argc, const char* argv[]){
  RNG.initialize(1234);

  std::vector<char const*> argv_filtered;

  int a = 0;
  int b = 1;
  bool fake_data = false;
  bool fake_data_errinflate = false;
  bool fake_1exp =false;
  bool fit_1exp = false; //freeze excited state matrix element to 0
  bool fake_exc_neg = false; //excited state has negative amplitude

  int i=0;
  while(i < argc){
    if(std::string(argv[i]) == "-ab"){
      a = strToAny<int>(argv[i+1]);
      b = strToAny<int>(argv[i+2]);
      i+=3;
    }else if(std::string(argv[i]) == "-fake"){
      fake_data = true;
      i++;
    }else if(std::string(argv[i]) == "-fake_errinflate"){
      fake_data_errinflate = true;
      i++;
    }else if(std::string(argv[i]) == "-fake_1exp"){
      fake_1exp = true;
      i++;
    }else if(std::string(argv[i]) == "-fit_1exp"){
      fit_1exp = true;
      i++;
    }else if(std::string(argv[i]) == "-fake_exc_neg"){
      fake_exc_neg = true;
      i++;
    }else{
      argv_filtered.push_back(argv[i++]);
    }
  }
  
  jackknifeDistributionD Epipi, Eexc, cpipi, cexc;
  {
    std::vector<jackknifeDistributionD> tmp;
    readParamsStandard(tmp, "pipi_fitparams.hdf5");
    cpipi = tmp[0];
    cexc = tmp[1];
    Epipi = tmp[4];
    Eexc = tmp[5];
  }  
  jackknifeDistributionD cK, mK;
  {
    std::vector<jackknifeDistributionD> tmp;
    readParamsStandard(tmp, "mK_fitparams.hdf5");
    cK = sqrt(tmp[0]);
    mK = tmp[1];
  }

  std::cout << "cpipi = " << cpipi << std::endl;
  std::cout << "Epipi = " << Epipi << std::endl;
  std::cout << "cexc = " << cexc << std::endl;
  std::cout << "Eexc = " << Eexc << std::endl;
  std::cout << "cK = " << cK << std::endl;
  std::cout << "mK = " << mK << std::endl;

  int nsample = mK.size();

  printMem("Beginning of execution");
  
  SampleAMAargs args;
  if(argc < 2){
    std::ofstream of("template.args");
    of << args;
    std::cout << "Wrote template argument file to template.args\n";
    return 0;
  }
  const std::string arg_file = argv[1];
  parse(args, arg_file);
  
  SampleAMAcmdLine cmdline(argv_filtered.size(),argv_filtered.data(),2);
  freezeCheck<>(cmdline);

  allInputs inputs(args,cmdline);

  readKtoPiPiAllDataSampleAMAoptions data_opt;
  data_opt.importGlobalOptions(cmdline);
  data_opt.read_opts = cmdline.getSampleAMAreadOptions();

  if(cmdline.plot_only){
    std::cout << "Plotting results and exiting\n";
    assert(cmdline.load_amplitude_data);
    PlotOnlyInputs pi(args.Lt, args.tmin_k_op, args.tmin_op_pi, cmdline.load_amplitude_data_file);
    fitfuncCall<PlotOnlyInputs,PlotOnlyCall>(args.fitfunc,pi);
    std::cout << "Done" << std::endl;
    return 0;
  }

  std::vector<std::string> data_file_fmt_sloppy =
    { "traj_<TRAJ>_type1_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>_mom<MOM>",
      "traj_<TRAJ>_type2_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>",
      "traj_<TRAJ>_type3_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>",
      "traj_<TRAJ>_type4" };

  std::vector<std::string> data_file_fmt_exact(data_file_fmt_sloppy);
  for(int i=0;i<data_file_fmt_exact.size();i++) data_file_fmt_exact[i] = data_file_fmt_exact[i] + cmdline.symmetric_quark_momenta_figure_file_extension;

  std::vector<std::pair<threeMomentum, double> > type1_pimom_proj = {  { {1,1,1}, 1.0/8.0 }, { {-1,-1,-1}, 1.0/8.0 },  { {-1,1,1}, 3.0/8.0 }, { {1,-1,-1}, 3.0/8.0 }  };

  std::string bubble_file_fmt_sloppy = "traj_<TRAJ>_FigureVdis_sep<TSEP_PIPI>_mom<PB>";
  std::string bubble_file_fmt_exact =  bubble_file_fmt_sloppy + "_symm";

  std::vector<std::pair<threeMomentum, double> > bubble_pimom_proj =  {  { {1,1,1}, 1.0/8.0 }, { {-1,-1,-1}, 1.0/8.0 },  
								         { {-1,1,1}, 1.0/8.0 }, { {1,-1,-1}, 1.0/8.0 }, 
								         { {1,-1,1}, 1.0/8.0 }, { {-1,1,-1}, 1.0/8.0 }, 
								         { {1,1,-1}, 1.0/8.0 }, { {-1,-1,1}, 1.0/8.0 } };

  //Prepare the data
  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > A0_all_j(10);
  std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > A0_all_dj(10);
  printMem("Prior to getData");
    
  if(cmdline.checkpoint_and_exit){
    checkpointRawOnly(args.tsep_k_pi, 
		      bubble_file_fmt_sloppy, bubble_file_fmt_exact, bubble_pimom_proj,
		      data_file_fmt_sloppy, data_file_fmt_exact, type1_pimom_proj,
		      args.data_dir_S, args.traj_start_S, args.traj_lessthan_S,
		      args.data_dir_C, args.traj_start_C, args.traj_lessthan_C,
		      args.traj_inc, args.bin_size, args.Lt, args.tsep_pipi, data_opt.read_opts);
    
    return 0;
  }

  // double M0_fake[10] = {-0.00302, 0.00273, 0.00459, 0.01036,
  // 		   -0.01271, -0.03366, 0.02259, 0.07907,
  // 		   -0.00678, -0.00090 };
  // double M1_fake[10] = {0.01084, 0.00373, -0.04050, -0.04755,
  // 		   0.05435, 0.11971, 0.00206, -0.05042,
  // 		   0.03636, 0.02769};

  std::vector<double> M0_fake(10,1.);
  std::vector<double> M1_fake(10,fake_exc_neg ? -1. : 1.);

  if(!fake_data){
    getDataSampleAMA(A0_all_j, A0_all_dj, args.tsep_k_pi, 
		     bubble_file_fmt_sloppy, bubble_file_fmt_exact, bubble_pimom_proj,
		     data_file_fmt_sloppy, data_file_fmt_exact, type1_pimom_proj,
		     args.data_dir_S, args.traj_start_S, args.traj_lessthan_S,
		     args.data_dir_C, args.traj_start_C, args.traj_lessthan_C,
		     args.traj_inc, args.bin_size, args.Lt, args.tsep_pipi, data_opt);
  }else{ //Generate fake data for testing/exploration
    int tseps[] = {10,12,14,16,18};

    jackknifeDistributionD Ethird(nsample, 1.);
    std::vector<double> M3(10,1.);

    for(int q=0;q<10;q++){
      for(int tsepidx=0;tsepidx<5;tsepidx++){
	int tsep = tseps[tsepidx];
	std::cout << "Fake data for Q" << q+1 << " and tsep=" << tsep << std::endl;
	for(int t=0;t<tsep;t++){
	  //Inflate the errors as a function of tsep-t to make more realistic
	  double infl =  exp(0.2 * (tsep-t));

	  jackknifeDistributionD val_pre = M0_fake[q] * cK * cpipi * exp(-mK * t) * exp(-Epipi * (tsep - t) );
	  if(!fake_1exp) val_pre = val_pre + M1_fake[q] * cK * cexc * exp(-mK * t) * exp(-Eexc * (tsep - t) );
	  
	  if(fake_data_errinflate){
	    double cen = val_pre.best();
	    jackknifeDistributionD val(nsample, [&](const int j){ return cen + infl * (val_pre.sample(j) - cen); });	
	    std::cout << t << " " << val << " (val_pre " << val_pre << " infl " << infl << ")" << std::endl;	
	    A0_all_j[q].push_back(amplitudeDataCoord(t,tsep), val);
	  }else{
	    std::cout << t << " " << val_pre << std::endl;	
	    A0_all_j[q].push_back(amplitudeDataCoord(t,tsep), val_pre);	    
	  }

	}
      }
    }
  }

  printMem("Prior to fitting");

  for(int q=0;q<10;q++){
    std::cout << "Starting Q" << q+1 << std::endl;

    //Sum data over top
    std::map<int, std::vector<int> > data_map;
    for(int d=0;d<A0_all_j[q].size();d++)
      data_map[ A0_all_j[q].coord(d).tsep_k_pi ].push_back(d);
    
    std::vector<int> tseps;
    for(auto it = data_map.begin(); it != data_map.end(); it++) tseps.push_back(it->first);
    std::sort(tseps.begin(),tseps.end());

    correlationFunction<double, jackknifeDistributionD> summed_corr; 
    summed_corr = getw<jackknifeDistributionD>(tseps,data_map,A0_all_j,q,a,b);

    correlationFunction<double, doubleJackknifeDistributionD> summed_corr_dj;
    if(!fake_data) summed_corr_dj = getw<doubleJackknifeDistributionD>(tseps,data_map,A0_all_dj,q,a,b);

    typedef typename composeFitPolicy<FitSum, frozenFitFuncPolicy, correlatedFitPolicy>::type FitPolicy;
    typedef typename FitPolicy::FitParameterDistribution FitParamsJack;
    typedef FitSum::Params FitParams;
    FitSum fitfunc(a,b);

    if(fake_data){//Check fit function
      FitParamsJack exact(nsample, FitParams(8));
      distributionStructPoke(exact, jackknifeDistributionD(nsample, M0_fake[q]), 0);
      distributionStructPoke(exact, jackknifeDistributionD(nsample, M1_fake[q]), 1);
      distributionStructPoke(exact, cK, 2);
      distributionStructPoke(exact, mK, 3);
      distributionStructPoke(exact, cpipi, 4);
      distributionStructPoke(exact, Epipi, 5);
      distributionStructPoke(exact, cexc, 6);
      distributionStructPoke(exact, Eexc, 7);

      std::cout << "Checking fitfunc against fake data for Q=" << q+1 << std::endl;
      for(int tt=0;tt<5;tt++){
	int tsep = tseps[tt];
	jackknifeDistributionD ffval(nsample, [&](const int b){ return fitfunc.value(tsep, exact.sample(b)); });
	const jackknifeDistributionD & rlval = summed_corr.value(tt);
	std::cout << tsep << " " << rlval << " " << ffval << std::endl;
      }
    }

    typename fitter<FitPolicy>::minimizerParamsType mlparams;
    mlparams.verbose = true;

    fitter<FitPolicy> fit(mlparams);
    fit.importFitFunc(fitfunc);

    
    NumericVector<jackknifeDistributionD> sigma_fake(tseps.size()); //for fake data
    NumericSquareMatrix<jackknifeDistributionD> inv_corr_fake(tseps.size());

    //for real data
    importCostFunctionParameters<correlatedFitPolicy, FitPolicy> *importer;
    if(fake_data){
      for(int i=0;i<tseps.size();i++){
	double err = summed_corr.value(i).standardError();
	sigma_fake[i] = jackknifeDistributionD(nsample, err);
	for(int j=0;j<tseps.size();j++)
	  inv_corr_fake(i,j) = jackknifeDistributionD(nsample, i==j ? 1. : 0.);
      }
      
      fit.importCostFunctionParameters(inv_corr_fake, sigma_fake);      
    }else{
      importer = new importCostFunctionParameters<correlatedFitPolicy, FitPolicy>(fit,summed_corr_dj);
    }

    //Params   M0 M1 cK mK cpipi Epipi cexc Eexc
    FitParamsJack freeze(nsample, FitParams(8));
    distributionStructPoke(freeze, cK, 2);
    distributionStructPoke(freeze, mK, 3);
    distributionStructPoke(freeze, cpipi, 4);
    distributionStructPoke(freeze, Epipi, 5);
    distributionStructPoke(freeze, cexc, 6);
    distributionStructPoke(freeze, Eexc, 7);

    if(fit_1exp){
      distributionStructPoke(freeze, jackknifeDistributionD(nsample,0.), 1);
      fit.freeze({1,2,3,4,5,6,7},freeze);
    }else{
      fit.freeze({2,3,4,5,6,7},freeze);
    }

    FitParamsJack result(nsample, FitParams(8,.5));
    jackknifeDistributionD chisq(nsample), chisq_per_dof(nsample);
    fit.fit(result, chisq, chisq_per_dof, summed_corr);

    std::cout << "Q" << q+1 << " result:\n" << result << "\nchi^2: " << chisq << "\nchi^2/dof: " << chisq_per_dof << std::endl;

    writeParamsStandard(result, stringize("result_Q%d.hdf5",q+1));
    writeParamsStandard(chisq, stringize("chisq_Q%d.hdf5",q+1));
    writeParamsStandard(chisq_per_dof, stringize("chisq_per_dof_Q%d.hdf5",q+1));

    int npoints = 60;
    int tsep_start = tseps[0];
    int tsep_end = tseps[tseps.size()-1];

    double delta = double(tsep_end - tsep_start)/double(npoints - 1);
    correlationFunction<double, jackknifeDistributionD> fit_curve;
    for(int i=0;i<npoints;i++){
      double t = tseps.front() + delta*i;
      jackknifeDistributionD val(nsample,
				 [&](const int s){
				   return fitfunc.value(t, result.sample(s));
				 });
      fit_curve.push_back(t,val);
    }    

    //Divide out ground-state time dep
    std::cout << "Dividing out ground-state time dependence for plotting:\n";
    for(int i=0;i<summed_corr.size();i++){
      jackknifeDistributionD tdep = groundTimeDependence(cK,mK, cpipi,Epipi,a,b, summed_corr.coord(i));

      std::cout << summed_corr.coord(i) << " " << summed_corr.value(i) << " -> ";

      summed_corr.value(i) = summed_corr.value(i)/tdep;

      std::cout << summed_corr.value(i) << std::endl;
    }
    for(int i=0;i<fit_curve.size();i++){
      jackknifeDistributionD tdep = groundTimeDependence(cK,mK, cpipi,Epipi,a,b, fit_curve.coord(i));
      fit_curve.value(i) = fit_curve.value(i)/tdep;
    }

    typedef DataSeriesAccessor<correlationFunction<double, jackknifeDistributionD>, 
			       ScalarCoordinateAccessor<double>, 
			       DistributionPlotAccessor<jackknifeDistributionD> > accessor;

    MatPlotLibScriptGenerate plot;
    typename MatPlotLibScriptGenerate::kwargsType kwargs;
    kwargs["color"] = "r";
    plot.plotData(accessor(summed_corr),kwargs);
    kwargs["alpha"] = 0.4;
    plot.errorBand(accessor(fit_curve),kwargs);

    std::string stub = stringize("plot_norm_sum_Q%d",q+1);
    plot.write(stub + ".py", stub + ".pdf");

    if(!fake_data) delete importer;
  }

  std::cout << "Done" << std::endl;
  
  return 0;
}

#endif

