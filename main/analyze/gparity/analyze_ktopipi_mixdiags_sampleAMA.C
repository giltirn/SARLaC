#if 1

//Code needs to be updated and fixed

int main(void){
  return 0;
}


#else

#include<ktopipi_common/ktopipi_common.h>
#include<ktopipi_sampleAMA_common/ktopipi_sampleAMA_common.h>

//#define PRINT_CORRECTION

using namespace SARLaC;

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

template<typename FitFunc>
struct Plot{};

template<>
struct Plot<FitKtoPiPiTwoExp>{
  typedef typename FitKtoPiPiTwoExp::Params Params;

  static void plot(const FitKtoPiPiTwoExp &fitfunc, const jackknifeDistribution<Params> &params, const allInputs &inputs,
		   const correlationFunction<amplitudeDataCoord, jackknifeDistributionD> &data_j){
    const SampleAMAargs &args = inputs.args;
    const SampleAMAcmdLine &cmdline = inputs.cmdline;

    //Plot results
    jackknifeDistributionD AK = distributionStructPeek(params, 0);
    jackknifeDistributionD mK = distributionStructPeek(params, 1);
    jackknifeDistributionD Cpipi = distributionStructPeek(params, 2);
    jackknifeDistributionD Epipi = distributionStructPeek(params, 3);

    const int nsample = AK.size();

    MatPlotLibScriptGenerate plot;
  
    typedef DataSeriesAccessor< correlationFunction<double, jackknifeDistributionD>, 
				ScalarCoordinateAccessor<double>, 
				DistributionPlotAccessor<jackknifeDistributionD> > accessor;  

    for(int tsep_k_pi_idx=0;tsep_k_pi_idx<args.tsep_k_pi.size();tsep_k_pi_idx++){
      int tsep_k_pi = args.tsep_k_pi[tsep_k_pi_idx];

      correlationFunction<double, jackknifeDistributionD> data_nrm;

      //Divide out ground-state time dependence and store as function of t_op_pi
      for(int i=0;i<data_j.size();i++) 
	if(data_j.coord(i).tsep_k_pi == tsep_k_pi && data_j.coord(i).t <= tsep_k_pi && data_j.coord(i).t >= args.tmin_k_op){
	  int tk_op = data_j.coord(i).t;
	  int top_pi = tsep_k_pi - tk_op;
	  jackknifeDistributionD tdep = AK*Cpipi*exp(-mK*tk_op)*exp(-Epipi*top_pi)/sqrt(2.);
	
	  data_nrm.push_back(top_pi, jackknifeDistributionD(data_j.value(i)/tdep));
	}

      plot.plotData(accessor(data_nrm), stringize("tsep_k_pi_%d",tsep_k_pi));
    }

    int tsep_k_pi_lrg = args.tsep_k_pi.back();

    correlationFunction<double, jackknifeDistributionD> fit_curve;
    int nplot = 60;
    double delta = double(tsep_k_pi_lrg - args.tmin_op_pi - args.tmin_k_op)/(nplot-1);

    for(int i=0;i<nplot;i++){
      double top_pi = args.tmin_op_pi + i*delta;
      double tk_op = tsep_k_pi_lrg - top_pi;

      jackknifeDistributionD val(nsample, [&](const int s){ return fitfunc.value(amplitudeDataCoord(tk_op,tsep_k_pi_lrg), params.sample(s)); });
      jackknifeDistributionD tdep = AK*Cpipi*exp(-mK*tk_op)*exp(-Epipi*top_pi)/sqrt(2.);
      val = val/tdep;
      fit_curve.push_back(top_pi, val);
    }

    plot.errorBand(accessor(fit_curve), "fit");

    plot.write("plot.py","plot.pdf");
  }
};

template<>
struct Plot<FitKtoPiPiTwoExpKaon>{
  typedef typename FitKtoPiPiTwoExpKaon::Params Params;

  static void plot(const FitKtoPiPiTwoExpKaon &fitfunc, const jackknifeDistribution<Params> &params, const allInputs &inputs,
		   const correlationFunction<amplitudeDataCoord, jackknifeDistributionD> &data_j){
    const SampleAMAargs &args = inputs.args;
    const SampleAMAcmdLine &cmdline = inputs.cmdline;

    //Plot results
    jackknifeDistributionD AK = distributionStructPeek(params, 0);
    jackknifeDistributionD mK = distributionStructPeek(params, 1);
    jackknifeDistributionD Cpipi = distributionStructPeek(params, 4);
    jackknifeDistributionD Epipi = distributionStructPeek(params, 5);

    const int nsample = AK.size();

    MatPlotLibScriptGenerate plot;
  
    typedef DataSeriesAccessor< correlationFunction<double, jackknifeDistributionD>, 
				ScalarCoordinateAccessor<double>, 
				DistributionPlotAccessor<jackknifeDistributionD> > accessor;  

    for(int tsep_k_pi_idx=0;tsep_k_pi_idx<args.tsep_k_pi.size();tsep_k_pi_idx++){
      int tsep_k_pi = args.tsep_k_pi[tsep_k_pi_idx];

      correlationFunction<double, jackknifeDistributionD> data_nrm;

      //Divide out ground-state time dependence and store as function of t_K_op
      for(int i=0;i<data_j.size();i++) 
	if(data_j.coord(i).tsep_k_pi == tsep_k_pi && data_j.coord(i).t <= tsep_k_pi - args.tmin_op_pi){
	  int tk_op = data_j.coord(i).t;
	  int top_pi = tsep_k_pi - tk_op;
	  jackknifeDistributionD tdep = AK*Cpipi*exp(-mK*tk_op)*exp(-Epipi*top_pi)/sqrt(2.);
	
	  data_nrm.push_back(tk_op, jackknifeDistributionD(data_j.value(i)/tdep));
	}

      plot.plotData(accessor(data_nrm), stringize("tsep_k_pi_%d",tsep_k_pi));
    }

    int tsep_k_pi_lrg = args.tsep_k_pi.back();

    correlationFunction<double, jackknifeDistributionD> fit_curve;
    int nplot = 60;
    double delta = double(tsep_k_pi_lrg - args.tmin_op_pi - args.tmin_k_op)/(nplot-1);

    for(int i=0;i<nplot;i++){
      double tk_op = args.tmin_k_op + i*delta;
      double top_pi = tsep_k_pi_lrg - tk_op;

      jackknifeDistributionD val(nsample, [&](const int s){ return fitfunc.value(amplitudeDataCoord(tk_op,tsep_k_pi_lrg), params.sample(s)); });
      jackknifeDistributionD tdep = AK*Cpipi*exp(-mK*tk_op)*exp(-Epipi*top_pi)/sqrt(2.);
      val = val/tdep;
      fit_curve.push_back(tk_op, val);
    }

    plot.errorBand(accessor(fit_curve), "fit");

    plot.write("plot.py","plot.pdf");
  }
};



template<>
struct Plot<FitKtoPiPiExcPiK>{
  typedef typename FitKtoPiPiExcPiK::Params Params;

  static void plot(const FitKtoPiPiExcPiK &fitfunc, const jackknifeDistribution<Params> &params, const allInputs &inputs,
		   const correlationFunction<amplitudeDataCoord, jackknifeDistributionD> &data_j){
    const SampleAMAargs &args = inputs.args;
    const SampleAMAcmdLine &cmdline = inputs.cmdline;

    //Plot results
    jackknifeDistributionD AK = distributionStructPeek(params, 0);
    jackknifeDistributionD mK = distributionStructPeek(params, 1);
    jackknifeDistributionD Cpipi = distributionStructPeek(params, 4);
    jackknifeDistributionD Epipi = distributionStructPeek(params, 5);

    const int nsample = AK.size();

    MatPlotLibScriptGenerate plot;
  
    typedef DataSeriesAccessor< correlationFunction<double, jackknifeDistributionD>, 
				ScalarCoordinateAccessor<double>, 
				DistributionPlotAccessor<jackknifeDistributionD> > accessor;  

    for(int tsep_k_pi_idx=0;tsep_k_pi_idx<args.tsep_k_pi.size();tsep_k_pi_idx++){
      int tsep_k_pi = args.tsep_k_pi[tsep_k_pi_idx];

      correlationFunction<double, jackknifeDistributionD> data_nrm;

      //Divide out ground-state time dependence and store as function of t_K_op
      for(int i=0;i<data_j.size();i++) 
	if(data_j.coord(i).tsep_k_pi == tsep_k_pi && data_j.coord(i).t <= tsep_k_pi - args.tmin_op_pi){
	  int tk_op = data_j.coord(i).t;
	  int top_pi = tsep_k_pi - tk_op;
	  jackknifeDistributionD tdep = AK*Cpipi*exp(-mK*tk_op)*exp(-Epipi*top_pi)/sqrt(2.);
	
	  data_nrm.push_back(tk_op, jackknifeDistributionD(data_j.value(i)/tdep));
	}

      plot.plotData(accessor(data_nrm), stringize("tsep_k_pi_%d",tsep_k_pi));
    }

    int tsep_k_pi_lrg = args.tsep_k_pi.back();

    correlationFunction<double, jackknifeDistributionD> fit_curve;
    int nplot = 60;
    double delta = double(tsep_k_pi_lrg - args.tmin_op_pi - args.tmin_k_op)/(nplot-1);

    for(int i=0;i<nplot;i++){
      double tk_op = args.tmin_k_op + i*delta;
      double top_pi = tsep_k_pi_lrg - tk_op;

      jackknifeDistributionD val(nsample, [&](const int s){ return fitfunc.value(amplitudeDataCoord(tk_op,tsep_k_pi_lrg), params.sample(s)); });
      jackknifeDistributionD tdep = AK*Cpipi*exp(-mK*tk_op)*exp(-Epipi*top_pi)/sqrt(2.);
      val = val/tdep;
      fit_curve.push_back(tk_op, val);
    }

    plot.errorBand(accessor(fit_curve), "fit");

    plot.write("plot.py","plot.pdf");
  }
};



template<typename FitFunc>
struct Fit{
  typedef typename FitFunc::Params Params;
  typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, correlatedFitPolicy>::type FitPolicies;
  
  static void fit(const allInputs &inputs,
		  const correlationFunction<amplitudeDataCoord, jackknifeDistributionD> &data_j,
		  const correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> &data_dj){
    const SampleAMAargs &args = inputs.args;
    const SampleAMAcmdLine &cmdline = inputs.cmdline;

    const int nsample = data_j.value(0).size();
    
    //Perform the fit
    correlationFunction<amplitudeDataCoord, jackknifeDistributionD> fit_data_j;
    correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> fit_data_dj;
    std::cout << "Data included in fit:\n";
    for(int i=0;i<data_j.size();i++){
      int t_k_op = data_j.coord(i).t;
      int t_op_pi = data_j.coord(i).tsep_k_pi - t_k_op;
      if(t_k_op >= args.tmin_k_op && t_op_pi >= args.tmin_op_pi){
	fit_data_j.push_back(data_j.coord(i), data_j.value(i));
	fit_data_dj.push_back(data_j.coord(i), data_dj.value(i));

	std::cout << "t_k_op=" << t_k_op << " t_op_pi=" << t_op_pi << " " << data_j.value(i) << std::endl;
      }
    }

    FitFunc fitfunc;

    fitter<FitPolicies> fit;
    fit.importFitFunc(fitfunc);
    std::unique_ptr<importCostFunctionParameters<correlatedFitPolicy,FitPolicies> > cost_corr;
    NumericSquareMatrix<jackknifeDistribution<double> > uncorr_invcorrmat;
    NumericVector<jackknifeDistribution<double> > uncorr_sigma;

    if(args.correlated) cost_corr.reset(new importCostFunctionParameters<correlatedFitPolicy,FitPolicies>(fit, fit_data_dj));
    else{
      uncorr_sigma = NumericVector<jackknifeDistribution<double> >(fit_data_dj.size(), 
								   [&](const int i){ return sqrt(doubleJackknifeDistributionD::covariance(fit_data_dj.value(i),fit_data_dj.value(i))); }
								   );
      uncorr_invcorrmat = NumericSquareMatrix<jackknifeDistribution<double> >(fit_data_dj.size(), 
									      [&](const int i, const int j){ return jackknifeDistributionD(nsample, i==j ? 1. : 0.); }
									      );
      fit.importCostFunctionParameters(uncorr_invcorrmat,uncorr_sigma);
    }


    std::cout << "Reading frozen fit params" << std::endl;
    assert(inputs.cmdline.load_freeze_data);
    readFrozenParams(fit, inputs.cmdline.freeze_data, nsample);

    jackknifeDistributionD chisq(nsample), chisq_per_dof(nsample);
    jackknifeDistribution<Params> params(nsample);
  
    std::cout << "Performing fit" << std::endl;
    fit.fit(params,chisq,chisq_per_dof,fit_data_j);
  
    std::cout << "Result: " << params << std::endl;
    std::cout << "Chisq: " << chisq << std::endl;
    std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;

#ifdef HAVE_HDF5
    writeParamsStandard(params, "params.hdf5");
    writeParamsStandard(chisq, "chisq.hdf5");
    writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");
#endif

    Plot<FitFunc>::plot(fitfunc, params, inputs, data_j);
  }
};


void fit(const allInputs &inputs,
	 const correlationFunction<amplitudeDataCoord, jackknifeDistributionD> &data_j,
	 const correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> &data_dj){
  switch(inputs.args.fitfunc){
  case KtoPiPiFitFunc::FitSeparateTwoExp:
    return Fit<FitKtoPiPiTwoExp>::fit(inputs, data_j, data_dj);
  case KtoPiPiFitFunc::FitSeparateTwoExpKaon:
    return Fit<FitKtoPiPiTwoExpKaon>::fit(inputs, data_j, data_dj);
  case KtoPiPiFitFunc::FitSeparateExcPiK:
    return Fit<FitKtoPiPiExcPiK>::fit(inputs, data_j, data_dj);
  default:
    error_exit(std::cout << "Unsupported fit func " << inputs.args.fitfunc << std::endl);
  }
};


class getDataFixedTsepKpi{
  const allInputs &inputs;
  const allBubbleData &bubble_data;
  const int tsep_k_pi;
  const int Lt;
  const int nsample;

  typedef NumericTensor<jackknifeDistributionD,1> CorrTypeJ;
  typedef NumericTensor<doubleJackknifeDistributionD,1> CorrTypeDJ;
  
  void computeVacSubs(std::vector<CorrTypeJ> &type4_vacsub_j, CorrTypeJ &mix4_vacsub_j, CorrTypeDJ &mix4_vacsub_dj, const allRawData &raw){
    //Compute mix4 vacuum subtraction. Also compute type-4 vacuum subtraction, not used for fit but for analysis
    const SampleAMAargs &args = inputs.args;
    const SampleAMAcmdLine &cmdline = inputs.cmdline;

    type4_vacsub_j = std::vector<CorrTypeJ>(10, NumericTensor<jackknifeDistributionD,1>({Lt}, jackknifeDistributionD(nsample,0.))); //[q]{t}
    mix4_vacsub_j = CorrTypeJ({Lt}, jackknifeDistributionD(nsample,0.));
    mix4_vacsub_dj = CorrTypeDJ({Lt}, doubleJackknifeDistributionD(nsample,0.));      

    std::vector<int> type4_nonzerotK = raw.raw_sloppy_S.nonzerotK(4);
    double nrm = type4_nonzerotK.size();
#pragma omp parallel for
    for(int t=0;t<Lt;t++){
      for(int i=0;i<type4_nonzerotK.size();i++){
	int tK = type4_nonzerotK[i];
	int tB = (tK + tsep_k_pi) % args.Lt;
		  	
#define RC(DIST,ARG) sampleAMAresampleCorrect<DIST>(raw.raw_sloppy_S. ARG, \
						    raw.raw_sloppy_C. ARG, \
						    raw.raw_exact_C. ARG, \
						    inputs.resamplers.resampler_S,inputs.resamplers.resampler_C);
	
  
	jackknifeDistributionD val_j = RC(jackknifeDistributionD, mix4_alltK_nobub({tK,t}));
	doubleJackknifeDistributionD val_dj = RC(doubleJackknifeDistributionD, mix4_alltK_nobub({tK,t}));

	val_j = val_j * bubble_data.bubble_j(&tB)/nrm;
	val_dj = val_dj * bubble_data.bubble_dj(&tB)/nrm;

	mix4_vacsub_j(&t) = mix4_vacsub_j(&t) + val_j;
	mix4_vacsub_dj(&t) = mix4_vacsub_dj(&t) + val_dj;

	for(int q=0;q<10;q++){
	  jackknifeDistributionD type4_j = RC(jackknifeDistributionD, A0_type4_alltK_nobub({q,tK,t}));
	  type4_j = type4_j * bubble_data.bubble_j(&tB)/nrm;
	  type4_vacsub_j[q](&t) = type4_vacsub_j[q](&t) + type4_j;
	}

#undef RC

      }
    }
  }
 
  struct mixNaccessor{
    int i;
    mixNaccessor(const int i): i(i){}
    const rawDataDistributionD & operator()(const int tK, const int t, const RawKtoPiPiData &raw) const{ return raw.mix_alltK(i)({tK,t}); }
    std::string descr(const int t) const{ return stringize("mix%d (tK avg) (t=%d)",i,t); }
    const std::vector<int> & nonzerotK(const RawKtoPiPiData &raw) const{ return raw.nonzerotK(i); }
  };

  template<typename DistributionType>
  void computeMatrixElem(correlationFunction<amplitudeDataCoord, DistributionType> &data, 
			 NumericTensor<DistributionType,1> &mix3, NumericTensor<DistributionType,1> &mix4,
			 const NumericTensor<DistributionType,1> &mix4_vacsub, const allRawData &raw){

    //Compute < <K|P|pipi> > via mix3 and mix4 diagrams
    mix3 = resampleAverageSampleAMA<DistributionType,mixNaccessor>(raw, Lt, inputs.resamplers, mixNaccessor(3));
    mix4 = resampleAverageSampleAMA<DistributionType,mixNaccessor>(raw, Lt, inputs.resamplers, mixNaccessor(4));

    //Perform vacuum subtraction and sum contributions
    mix4 = mix4 - mix4_vacsub;

    NumericTensor<DistributionType,1> kPpipi_full = mix3 + mix4;

    for(int t=0;t<Lt;t++){
      amplitudeDataCoord coor(t,tsep_k_pi);
      data.push_back(coor, kPpipi_full(&t));
    }
  }
  template<typename DistributionType>
  void computeMatrixElem(correlationFunction<amplitudeDataCoord, DistributionType> &data, 			 
		    const NumericTensor<DistributionType,1> &mix4_vacsub, const allRawData &raw){
    NumericTensor<DistributionType,1> mix3, mix4;
    computeMatrixElem(data,mix3,mix4,mix4_vacsub,raw);
  }

  static void write(const std::string &filename, 
		    const std::vector<CorrTypeJ> &what){
    std::vector<std::vector<jackknifeDistributionD> > tmp(what.size());
    for(int i=0;i<what.size();i++){
      tmp[i].resize(what[i].size(0));
      for(int j=0;j<tmp[i].size();j++)
	tmp[i][j] = what[i](&j);
    }
#ifdef HAVE_HDF5
    writeParamsStandard(tmp, filename);
#endif
  }
  
  struct alphaCptAccessor{
    int q, num_or_den; //num_or_den : 0 = type4,q,nobub   1 = mix4,nobub
    alphaCptAccessor(const int num_or_den, const int q = -1): q(q), num_or_den(num_or_den){}
    const rawDataDistributionD & operator()(const int tK, const int t, const RawKtoPiPiData &raw) const{ 
      return num_or_den ? raw.mix4_alltK_nobub({tK,t}) : raw.A0_type4_alltK_nobub({q,tK,t});
    }
    std::string descr(const int t) const{ return num_or_den ? stringize("mix4 nobub (tK avg) (t=%d)",t) : stringize("type4 Q%d nobub (tK avg) (t=%d)",q+1,t) ; }
    const std::vector<int> & nonzerotK(const RawKtoPiPiData &raw) const{ return raw.nonzerotK(4); }
  };
  struct typeDataAccessor{
    int q, type;
    typeDataAccessor(int q, int type): q(q), type(type){}
    const rawDataDistributionD & operator()(const int tK, const int t, const RawKtoPiPiData &raw) const{ 
      return raw.A0_alltK(type)({q,tK,t});
    }
    std::string descr(const int t) const{ return stringize("Q%d type %d (tK avg) (t=%d)",q,type,t); }
    const std::vector<int> & nonzerotK(const RawKtoPiPiData &raw) const{ return raw.nonzerotK(type); }
  };

  void analyzeSubtraction(const CorrTypeJ &mix3_j, const CorrTypeJ &mix4_j, const std::vector<CorrTypeJ> &type4_vacsub_j, const allRawData &raw){
    //While not used for the fit, it is worthwhile also computing and writing out the coefficients alpha and looking more closely at the size of the contraction
    std::vector<CorrTypeJ> alpha_j(10, CorrTypeJ({Lt}));

    CorrTypeJ mix4_nobub_srcavg = resampleAverageSampleAMA<jackknifeDistributionD,alphaCptAccessor>(raw, Lt, inputs.resamplers, alphaCptAccessor(1));
    for(int q=0;q<10;q++){
      std::cout << "alpha_" << q+1 << ":\n";
      CorrTypeJ A0_type4_nobub_srcavg = resampleAverageSampleAMA<jackknifeDistributionD,alphaCptAccessor>(raw, Lt, inputs.resamplers, alphaCptAccessor(0,q));
      for(int t=0;t<Lt;t++){
	alpha_j[q](&t) = A0_type4_nobub_srcavg(&t)/mix4_nobub_srcavg(&t);      
	std::cout << t << " " << alpha_j[q](&t) << std::endl;
      }
    }
    write(stringize("alpha_tsep_k_pi%d.hdf5",tsep_k_pi), alpha_j);

    std::vector<CorrTypeJ> A0_pre(10), A0_post(10), correction(10);

    for(int q=0;q<10;q++){
#define RA(type) resampleAverageSampleAMA<jackknifeDistributionD,typeDataAccessor>(raw, Lt, inputs.resamplers, typeDataAccessor(q,type))	  
      A0_pre[q] = RA(1) + RA(2) + RA(3) + RA(4) - type4_vacsub_j[q];
#undef RA
      correction[q] = CorrTypeJ({Lt},[&](const int *t){ return alpha_j[q](t)*(mix3_j(t) + mix4_j(t)); });
      A0_post[q] = A0_pre[q] - correction[q];
    }
      
    write(stringize("A0_presub_tsep_k_pi%d.hdf5",tsep_k_pi), A0_pre);
    write(stringize("A0_postsub_tsep_k_pi%d.hdf5",tsep_k_pi), A0_post);
    write(stringize("sub_tsep_k_pi%d.hdf5",tsep_k_pi), correction);
  }


public:

  getDataFixedTsepKpi(correlationFunction<amplitudeDataCoord, jackknifeDistributionD> &data_j,
		      correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> &data_dj,
		      const allInputs &inputs, const allBubbleData &bubble_data, const int tsep_k_pi): inputs(inputs), bubble_data(bubble_data), 
												       tsep_k_pi(tsep_k_pi), Lt(inputs.args.Lt),
												       nsample(inputs.resamplers.nS + inputs.resamplers.nC){

    std::cout << "Getting data for tsep_k_pi = " <<  tsep_k_pi << std::endl;
    readKtoPiPiDataSampleAMAoptions opt = inputs.cmdline.getSampleAMAreadOptions();

    std::vector<std::string> data_file_fmt_sloppy =
      { "traj_<TRAJ>_type1_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>_mom<MOM>",
	"traj_<TRAJ>_type2_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>",
	"traj_<TRAJ>_type3_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>",
	"traj_<TRAJ>_type4" };

    std::vector<std::string> data_file_fmt_exact(data_file_fmt_sloppy);
    for(int i=0;i<data_file_fmt_exact.size();i++) data_file_fmt_exact[i] = data_file_fmt_exact[i] + inputs.cmdline.symmetric_quark_momenta_figure_file_extension;

    std::vector<std::pair<threeMomentum, double> > type1_pimom_proj = {  { {1,1,1}, 1.0/8.0 }, { {-1,-1,-1}, 1.0/8.0 },  { {-1,1,1}, 3.0/8.0 }, { {1,-1,-1}, 3.0/8.0 }  };

    allRawData raw(bubble_data, tsep_k_pi, 
		   data_file_fmt_sloppy, data_file_fmt_exact, type1_pimom_proj,
		   inputs.args.data_dir_S, inputs.args.traj_start_S, inputs.args.traj_lessthan_S,
		   inputs.args.data_dir_C, inputs.args.traj_start_C, inputs.args.traj_lessthan_C,
		   inputs.args.traj_inc, inputs.args.bin_size, Lt, inputs.args.tsep_pipi, opt);

    std::vector<CorrTypeJ> type4_vacsub_j; //[q]{t}
    CorrTypeJ mix4_vacsub_j;
    CorrTypeDJ mix4_vacsub_dj;

    computeVacSubs(type4_vacsub_j,mix4_vacsub_j,mix4_vacsub_dj, raw);

    CorrTypeJ mix3_j, mix4_j;
    computeMatrixElem<jackknifeDistributionD>(data_j, mix3_j, mix4_j, mix4_vacsub_j, raw);
    computeMatrixElem<doubleJackknifeDistributionD>(data_dj, mix4_vacsub_dj, raw);
  
    analyzeSubtraction(mix3_j, mix4_j, type4_vacsub_j, raw);
  }

};



int main(const int argc, const char* argv[]){
  SampleAMAargs args;
  if(argc < 2){
    std::ofstream of("template.args");
    of << args;
    std::cout << "Wrote template argument file to template.args\n";
    return 0;
  }
  const std::string arg_file = argv[1];
  parse(args, arg_file);
  
  SampleAMAcmdLine cmdline(argc,argv,2);
  freezeCheck<FreezeParams>(cmdline);

  allInputs inputs(args,cmdline);

  readKtoPiPiDataSampleAMAoptions opt = cmdline.getSampleAMAreadOptions();

  std::string bubble_file_fmt_sloppy = "traj_<TRAJ>_FigureVdis_sep<TSEP_PIPI>_mom<PB>";
  std::string bubble_file_fmt_exact =  bubble_file_fmt_sloppy + "_symm";

  std::vector<std::pair<threeMomentum, double> > bubble_pimom_proj =  {  { {1,1,1}, 1.0/8.0 }, { {-1,-1,-1}, 1.0/8.0 },  
								         { {-1,1,1}, 1.0/8.0 }, { {1,-1,-1}, 1.0/8.0 }, 
								         { {1,-1,1}, 1.0/8.0 }, { {-1,1,-1}, 1.0/8.0 }, 
								         { {1,1,-1}, 1.0/8.0 }, { {-1,-1,1}, 1.0/8.0 } };

  allBubbleData bubble_data(bubble_file_fmt_sloppy, bubble_file_fmt_exact, bubble_pimom_proj,
			    args.data_dir_S, args.traj_start_S, args.traj_lessthan_S,
			    args.data_dir_C, args.traj_start_C, args.traj_lessthan_C,
			    args.traj_inc, args.bin_size, args.Lt, args.tsep_pipi, inputs.resamplers, opt);

  correlationFunction<amplitudeDataCoord, jackknifeDistributionD> data_j;
  correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> data_dj;

  if(cmdline.load_amplitude_data){
#ifdef HAVE_HDF5
    HDF5reader reader(cmdline.load_amplitude_data_file);
    read(reader, data_j, "data_j");
    read(reader, data_dj, "data_dj");
#else
    error_exit("Reading amplitude data requires HDF5\n");
#endif
  }else{    
    for(int tsep_k_pi_idx=0;tsep_k_pi_idx<inputs.args.tsep_k_pi.size();tsep_k_pi_idx++){
      int tsep_k_pi = inputs.args.tsep_k_pi[tsep_k_pi_idx];
      getDataFixedTsepKpi getdata(data_j, data_dj, inputs, bubble_data, tsep_k_pi);
    }//tsep_k_pi loop
  }

  if(cmdline.save_amplitude_data){
#ifdef HAVE_HDF5
    HDF5writer writer(cmdline.save_amplitude_data_file);
    write(writer, data_j, "data_j");
    write(writer, data_dj, "data_dj");
#else
    error_exit("Saving amplitude data requires HDF5\n");
#endif
  }

  fit(inputs, data_j, data_dj);
  
  std::cout << "Done" << std::endl;

  return 0;
}

#endif
