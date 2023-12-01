#ifndef _FIT_KTOSIGMA_KTOPIPI_GPARITY_BOOTSTRAP_ALPHA_VARY_PLOT_H_
#define _FIT_KTOSIGMA_KTOPIPI_GPARITY_BOOTSTRAP_ALPHA_VARY_PLOT_H_

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

//Generate a plot a each matrix element on a particular timeslice and tsep_K_snk as a function of varying the coeffcient of alpha away from 1
template<typename DistributionType, typename BinResampler>
DistributionType computeQalphaCoeff(const RawData &raw, const int q, const PiPiOperator op, const int t,
				    const int tsep_k_snk, const int tsep_k_snk_idx, const int Lt, const BinResampler &bin_resampler,
				    const double alpha_coeff){
  computeQamplitudeOpts opt;
  opt.alpha_scale = alpha_coeff;
  
  switch(op){
  case PiPiOperator::PiPiGnd:
    return computeQamplitude<DistributionType>(q, tsep_k_snk, *raw.raw_ktopipi_gnd[tsep_k_snk_idx], *raw.bubble_data_gnd, Lt, "", bin_resampler, opt)(&t);
  case PiPiOperator::PiPiExc:
    return computeQamplitude<DistributionType>(q, tsep_k_snk, *raw.raw_ktopipi_exc[tsep_k_snk_idx], *raw.bubble_data_exc, Lt, "", bin_resampler, opt)(&t);
  case PiPiOperator::Sigma:
    return computeQamplitude<DistributionType>(q, tsep_k_snk, *raw.raw_ktosigma[tsep_k_snk_idx], *raw.bubble_data_sigma, Lt, "", bin_resampler, opt)(&t);
  }
}

template<typename DistributionType, typename ArgsType, typename CMDlineType, typename BinResampler>
void alphaVaryPlot(const RawData &raw, const PiPiOperator op, const int q, const int tsep_k_snk_idx, const int tsep_op_snk,		    
		   const ArgsType &args, const CMDlineType &cmdline, const BinResampler &bin_resampler, const std::string &file_stub,
		   const int nstep_each_side = 5, const double coeff_step = 0.01){
  
  int tsep_k_snk;
  if(op == PiPiOperator::PiPiGnd || op == PiPiOperator::PiPiExc) tsep_k_snk =  args.tsep_k_pi[tsep_k_snk_idx];
  else tsep_k_snk =  args.tsep_k_sigma[tsep_k_snk_idx];

  int t = tsep_k_snk - tsep_op_snk;
      
  int ndata = 2*nstep_each_side + 1;
  std::vector<DistributionType> yvals(ndata);
  std::vector<double> xvals(ndata);
  int ii=0;
  for(int i=nstep_each_side-1;i>=0;i--){
    double x = 1. - (i+1)*coeff_step;
    xvals[ii] = x;
    yvals[ii] = computeQalphaCoeff<DistributionType>(raw, q, op, t, tsep_k_snk, tsep_k_snk_idx, args.Lt, bin_resampler, x);
    ii++;
  }
  xvals[ii] = 1.0;
  yvals[ii] = computeQalphaCoeff<DistributionType>(raw, q, op, t, tsep_k_snk, tsep_k_snk_idx, args.Lt, bin_resampler, 1.);
  ii++;

  for(int i=0; i<nstep_each_side;i++){
    double x = 1. + (i+1)*coeff_step;
    xvals[ii] = x;
    yvals[ii] = computeQalphaCoeff<DistributionType>(raw, q, op, t, tsep_k_snk, tsep_k_snk_idx, args.Lt, bin_resampler, x);
    ii++;
  }

  MatPlotLibScriptGenerate plot;
  struct acc{
    const std::vector<DistributionType> &yvals;
    const std::vector<double> &xvals;
    acc(const std::vector<DistributionType> &yvals, const std::vector<double> &xvals): yvals(yvals), xvals(xvals){}
    int size() const{ return yvals.size(); }

    inline double x(const int i) const{ return xvals[i]; }
    inline double dxp(const int i) const{ return 0; }
    inline double dxm(const int i) const{ return 0; }  

    inline double y(const int i) const{ return yvals[i].best(); }
    inline double dyp(const int i) const{ return yvals[i].standardError(); }
    inline double dym(const int i) const{ return yvals[i].standardError(); }  
  };
  plot.plotData(acc(yvals,xvals));
  plot.write(file_stub + ".py", file_stub + ".pdf");
}

typedef std::pair<threeMomentum, double> momMultiplicityPair;
typedef std::vector<std::string> typeFileFormat;

#define ALPHA_VARY_PLOT_ARGS_MEMBERS							\
  (int, tsep_k_snk_idx)							\
  (int, tsep_op_snk)							\
  (int, nstep_each_side)						\
  (double, coeff_step)
  
struct AlphaVaryPlotArgs{
  GENERATE_MEMBERS(ALPHA_VARY_PLOT_ARGS_MEMBERS);

  AlphaVaryPlotArgs(): tsep_k_snk_idx(0), tsep_op_snk(5), nstep_each_side(5), coeff_step(0.01){}
};

GENERATE_PARSER(AlphaVaryPlotArgs, ALPHA_VARY_PLOT_ARGS_MEMBERS);


CPSFIT_END_NAMESPACE

#endif
