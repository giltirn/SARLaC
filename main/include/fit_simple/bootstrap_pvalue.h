#ifndef _FIT_SIMPLE_BOOT_PVALUE_
#define _FIT_SIMPLE_BOOT_PVALUE_

#include<fit/bootstrap_pvalue.h>

struct BootResampler{
  void operator()(std::vector<rawDataCorrelationFunctionD> &raw_data, const std::vector<int> &map) const{
    for(int c=0;c<raw_data.size();c++)
      for(int i=0;i<raw_data[c].size();i++)
	bootstrapResampleRaw(raw_data[c].value(i), map);
  }
};

struct BootFitter{
  const Args &args;
  const CMDline &cmdline;
  parameterVectorD guess;

  MarquardtLevenbergParameters<double> min_params;

  BootFitter(const Args &args, const CMDline &cmdline, const parameterVectorD &guess): args(args), cmdline(cmdline), guess(guess){
    min_params.verbose = false;
    min_params.max_iter = 1e6;
    min_params.lambda_max = 1e10;
    min_params.lambda_factor = 1.5;
  }

  double operator()(const std::vector<rawDataCorrelationFunctionD> &raw_data, const correlationFunction<double, double> &corrections, const int b) const{
    bool do_j_b, do_j_ub;
    getJtypes(do_j_b, do_j_ub, args.covariance_strategy);
    
    jackknifeCorrelationFunctionD data_j_unbinned, data_j_binned;
    if(do_j_b){
      data_j_binned = resampleAndCombine<jackknifeDistributionD>(raw_data, args.Lt, args.bin_size, args.combination, args.outer_time_dep);
      recenter(data_j_binned, corrections);
    }
    if(do_j_ub){
      data_j_unbinned = resampleAndCombine<jackknifeDistributionD>(raw_data, args.Lt, 1, args.combination, args.outer_time_dep);
      recenter(data_j_unbinned, corrections);      
    }

    parameterVectorD cen_params = guess;
    double cen_chisq;
    int cen_dof;
    
    fitCentral(cen_params, cen_chisq, cen_dof, data_j_binned, data_j_unbinned, args, cmdline, &min_params);

    return cen_chisq;
  }
};

double computeBootstrapPvalue(const parameterVectorD &base_params,
			      const double base_chisq, int dof,
			      const std::vector<rawDataCorrelationFunctionD> &channels_raw,
			      const jackknifeCorrelationFunctionD &data_j_binned,
			      const Args &args, const CMDline &cmdline){
  if(!RNG.isInitialized()) RNG.initialize(1234);  

  int nsample = channels_raw[0].value(0).size();
  int ndata = data_j_binned.size();

  std::unique_ptr< FitFuncManagerBase<Args,CMDline> > fitfunc_manager = getFitFuncManager(args, cmdline);
  genericFitFuncBase const* fitfunc = fitfunc_manager->getFitFunc();

  correlationFunction<double, double> fit_data_cen(ndata);
  correlationFunction<double, double> fit_vals(ndata);

  std::cout << "For mean parameters, comparison of fit value and data:" << std::endl;
  for(int t=0;t<ndata;t++){
    fit_data_cen.coord(t) = fit_vals.coord(t) = data_j_binned.coord(t);
    fit_data_cen.value(t) = data_j_binned.value(t).mean();
    fit_vals.value(t) = fitfunc->value( genericFitFuncBase::getWrappedCoord(data_j_binned.coord(t)), base_params );
    
    std::cout << data_j_binned.coord(t) << " fit " << fit_vals.value(t) << " data " << fit_data_cen.value(t) << " diff " << fit_vals.value(t)-fit_data_cen.value(t) <<  std::endl;
  }
  
  BootResampler resampler;
  BootFitter fitter(args, cmdline, base_params);
  
  {
    std::cout << "Test fit functor on original data: this should agree with the fit to the central values" << std::endl;
    correlationFunction<double, double> corrections(ndata, [&](const int i){ return typename correlationFunction<double, double>::ElementType(fit_data_cen.coord(i), 0.); });

    fitter(channels_raw, corrections, 0);

    std::cout << "Test fit functor on recentered original data: this should have a chisq that is close to zero because the data have been adjusted to match the fit" << std::endl;
    for(int i=0;i<ndata;i++) corrections.value(i) = fit_vals.value(i) - fit_data_cen.value(i);
    fitter(channels_raw, corrections, 0);

    std::vector<int> map(nsample);
    for(int i=0;i<nsample;i++) map[i] = nsample - 1 - i;
    
    std::vector<rawDataCorrelationFunctionD> channels_raw_r(channels_raw);
    resampler(channels_raw_r, map);
    assert(channels_raw_r[0].value(0).sample(0) = channels_raw[0].value(0).sample(nsample-1));
	   
    std::cout << "Test the resampler by applying a reflection in sample index. The result should be the same as the fit to the recentered original data" << std::endl;
    fitter(channels_raw_r, corrections, 0);

    std::cout << "Get a couple of maps from a non-overlapping block bootstrap with block size " << args.bin_size << " and get fit results. This should give you an idea of the typical values of chisq in the true distribution. If these are larger than chisq of your fit you can expect a good p-value." << std::endl;

    std::vector<std::vector<int> > bmap = nonoverlappingBlockResampleTable(RNG, nsample, args.bin_size, 5);
    for(int i=0;i<5;i++){
      channels_raw_r = channels_raw;
      resampler(channels_raw_r, bmap[i]);
      fitter(channels_raw_r, corrections, i);
    }
  }
  std::cout << "Computing bootstrap p-value" << std::endl;

  int block_size = args.bin_size;
  std::vector<double> q2_boot_p;
  double boot_p = bootstrapPvalue(base_chisq, channels_raw, nsample, fit_data_cen, fit_vals, resampler, fitter, 1000, BootResampleTableType::NonOverlappingBlock, block_size, -1, RNG, &q2_boot_p);

  rawDataDistributionD boot_q2(q2_boot_p.size(), [&](const int s){ return q2_boot_p[s]; });
  
  std::cout << "Bootstrap distribution of q^2 has mean " << boot_q2.mean() << " and std.dev " << boot_q2.standardDeviation() << std::endl;
  
  writeParamsStandard(boot_q2, "boot_q2.hdf5");

  return boot_p;
}

#endif
