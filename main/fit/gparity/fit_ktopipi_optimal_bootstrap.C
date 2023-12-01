#include<ktopipi_common/ktopipi_common.h>

#include<fit_ktopipi_gnd_exc_ktosigma_combined_gparity/resampled_data.h>

#include<fit_ktopipi_optimal_bootstrap/args.h>
#include<fit_ktopipi_optimal_bootstrap/cmdline.h>

using namespace CPSfit;

template<typename DistributionType>
std::vector<correlationFunction<amplitudeDataCoord, DistributionType> > projectData(const std::vector<PiPiOperator> &ops, const std::vector<bootstrapDistributionD> &r,
										    const ResampledData<DistributionType> &data){
  assert(r.size() == ops.size());
  
  DistributionType zero(data(ops[0])[0].value(0)); zeroit(zero);

  std::vector<correlationFunction<amplitudeDataCoord, DistributionType> > out(data(ops[0]));

  for(int q=0;q<out.size();q++){
    for(int tt=0;tt<out[q].size();tt++){
      DistributionType &v = out[q].value(tt);
      v = zero;
      
      for(int op=0;op<ops.size();op++){
	assert(data(ops[op])[q].coord(tt) == out[q].coord(tt));
	DistributionType d = data(ops[op])[q].value(tt);
	d.best() = d.best() * r[op].best();
	for(int s=0;s<d.size();s++) d.sample(s) = d.sample(s) * r[op].sample(s);
	
	v = v+d;
      }
    }
  }

  return out;
}

template<typename DistributionType>
std::vector<correlationFunction<amplitudeDataCoord, DistributionType> > getDataInRange(const std::vector<correlationFunction<amplitudeDataCoord, DistributionType> > &data,
										       const int tmin_k_op, const int tmin_op_pi){
  int nq = data.size();
  std::vector<correlationFunction<amplitudeDataCoord, DistributionType> > out(nq);
  for(int q=0;q<nq;q++){
    for(int i=0;i<data[q].size();i++){
      auto &c = data[q].coord(i);
      int t = (int)c.t;
      int top_pi = c.tsep_k_pi - t;
      if(t>= tmin_k_op && top_pi >= tmin_op_pi)
	out[q].push_back(data[q][i]);
    }
  }
  return out;
}


int main(const int argc, const char* argv[]){
  printMem("Beginning of execution");
  
  Args args;
  if(argc < 2){
    std::ofstream of("template.args");
    of << args;
    std::cout << "Wrote template argument file to template.args\n";
    return 0;
  }
  const std::string arg_file = argv[1];
  parse(args, arg_file);

  CMDline cmdline(argc,argv,2); 

  if(cmdline.load_boot_resample_table) std::cout << "Loading boot resample table" << std::endl;
  else std::cout << "*Not* loading boot resample table" << std::endl;

  //Compute r
  std::vector<bootstrapDistributionD> r = computeR(args.op_amplitudes, args.Ascale); 

  //Read AK, mK and Epipi
  bootstrapDistributionD cK, mK, E0;
  { 
    std::vector<bootstrapDistributionD> p;
    readParamsStandard(p,  args.input_params.kaon2pt_fit_result);
    mK = p[args.input_params.idx_mK];
    cK = sqrt( p[args.input_params.idx_cK] );
  }
  {
    std::vector<bootstrapDistributionD> p;
    readParamsStandard(p, args.input_params.pipi_sigma_sim_fit_result);
    E0 = p[args.input_params.idx_E0];
  }
  std::cout << "Read inputs:" << std::endl;
  std::cout << "cK = " << cK << std::endl;
  std::cout << "mK = " << mK << std::endl;
  std::cout << "E0 = " << E0 << std::endl;
   

  //Read raw data
  RawData raw;
  if(!cmdline.load_resampled_data_container_checkpoint){
    if(cmdline.load_raw_data_container_checkpoint){
      std::cout << "Reading raw data container from checkpoint file" << std::endl;
      HDF5reader rd(cmdline.load_raw_data_container_checkpoint_file);
      read(rd, raw, "raw_data_container");
    }else{
      std::cout << "Reading raw data" << std::endl;
      raw.read(args, cmdline);
    }  
    if(cmdline.save_raw_data_container_checkpoint){
      std::cout << "Saving raw data container to checkpoint file" << std::endl;
      HDF5writer wr(cmdline.save_raw_data_container_checkpoint_file);
      write(wr, raw, "raw_data_container");
    }
  }
    
  std::cout << "Computing resampled data" << std::endl;
  ResampledData<bootstrapDistributionD> data_b;
  ResampledData<bootJackknifeDistributionD> data_bj;
  
  if(cmdline.load_resampled_data_container_checkpoint){
    std::cout << "Reading resampled data container from checkpoint file" << std::endl;
    HDF5reader rd(cmdline.load_resampled_data_container_checkpoint_file);
    read(rd, data_b, "data_b");
    read(rd, data_bj, "data_bj");
  }else{
    int nsample = raw.nsample();

    //Generate the resample table
    resampleTableOptions ropt;
    ropt.read_from_file = cmdline.load_boot_resample_table;
    ropt.read_file = cmdline.load_boot_resample_table_file;
    ropt.write_to_file = cmdline.save_boot_resample_table;
    ropt.write_file = cmdline.save_boot_resample_table_file;
    
    std::vector<std::vector<int> > rtable = generateResampleTable(nsample, args.nboot, args.resample_table_type, args.block_size, RNG, ropt);
    if(rtable.size() != args.nboot) error_exit(std::cout << "Expected resample table of size " << args.nboot << ", got " << rtable.size() << std::endl);

    bootstrapBlockResampler resampler(rtable);

    data_b.resample(raw, args, cmdline, "bootstrap", resampler);
    data_bj.resample(raw, args, cmdline, "boot-jackknife", resampler);
  }
  
  if(cmdline.save_resampled_data_container_checkpoint){
    std::cout << "Writing resampled data container to checkpoint file" << std::endl;
    HDF5writer wr(cmdline.save_resampled_data_container_checkpoint_file);
    write(wr, data_b, "data_b");
    write(wr, data_bj, "data_bj");
  }

  if(args.basis == Basis::Basis7){
    std::cout << "Converting to 7-basis" << std::endl;
    data_b.convertBasis10to7(); 
    data_bj.convertBasis10to7();
  }

  //Constrain tsep_k_snk to those in input file (in case we read resampled data that contains more)
  std::map<PiPiOperator, std::vector<int> const*> op_tsep_list = {  {PiPiOperator::PiPiGnd, &args.tsep_k_pi}, {PiPiOperator::PiPiExc, &args.tsep_k_pi}, {PiPiOperator::Sigma, &args.tsep_k_sigma} };

  for(int o=0;o<args.operators.size();o++){
    data_b.constrainSourceSinkSep(args.operators[o], *op_tsep_list[args.operators[o]]);
    data_bj.constrainSourceSinkSep(args.operators[o], *op_tsep_list[args.operators[o]]);
  }

  typedef std::vector<correlationFunction<amplitudeDataCoord, bootstrapDistributionD> > CorrFuncAllQB;
  typedef std::vector<correlationFunction<amplitudeDataCoord, bootJackknifeDistributionD> > CorrFuncAllQBJ;

  CorrFuncAllQB proj_data_b = projectData(args.operators, r, data_b);
  CorrFuncAllQBJ proj_data_bj = projectData(args.operators, r, data_bj);
  
  CorrFuncAllQB proj_data_inrange_b = getDataInRange(proj_data_b, args.tmin_k_op, args.tmin_op_snk);
  CorrFuncAllQBJ proj_data_inrange_bj = getDataInRange(proj_data_bj, args.tmin_k_op, args.tmin_op_snk);

  typedef FitKtoPiPi FitFunc;
  typedef FitFunc::Params Params;

  FitFunc fitfunc;
  Params guess;  
  genericFitFuncWrapper<FitFunc> fwrap(fitfunc,guess);
  
  bootstrapDistributionD one(proj_data_b[0].value(0).getInitializer(), 1.);

  int nq = proj_data_b.size();
  
  int nboot = one.size();

  std::vector<bootstrapDistributionD> chisq_all(nq);
  std::vector<bootstrapDistributionD> chisq_per_dof_all(nq); 
  std::vector<bootstrapDistribution<Params> > params_all(nq);
  std::vector<double> p_boot_all(nq);
  std::vector<rawDataDistributionD> q2_dist_r_all(nq);

  for(int q=0;q<nq;q++){
    std::cout << "Fitting Q" << q+1 << std::endl;

    simpleFitWrapper<bootstrapDistributionD> fitter(fwrap, MinimizerType::MarquardtLevenberg);
    fitter.generateCovarianceMatrix(proj_data_inrange_bj[q], args.correlated ? CostType::Correlated : CostType::Uncorrelated);

    //Freeze amplitude to 1
    fitter.freeze(2, one);

    //Freeze AK, mK and Epipi
    fitter.freeze(0, cK);
    fitter.freeze(1, mK);
    fitter.freeze(3, E0);

    bootstrapDistribution<Params> params(one.getInitializer(),guess);
    bootstrapDistributionD chisq(one), chisq_per_dof(one);
    int dof;

    fitter.fit(params, chisq, chisq_per_dof, dof, proj_data_inrange_b[q]);
    
    std::cout << "Q" << q+1 << " params:\n";
    {
      for(int p=0;p<params.sample(0).size();p++){
	std::string tag = params.sample(0).memberName(p);	
	bootstrapDistributionD tmp;
	standardIOhelper<bootstrapDistributionD, bootstrapDistribution<Params> >::extractStructEntry(tmp, params, p);
	std::cout << tag << " = " << tmp << std::endl;
      }
    }
    
    std::cout << std::endl;
    std::cout << "Chisq: " << chisq << std::endl;
    std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;

    //Compute bootstrap p-value
    std::cout << "Computing bootstrap p-value" << std::endl;

    correlationFunction<amplitudeDataCoord, bootstrapDistributionD> proj_data_inrange_b_rc(proj_data_inrange_b[q]);

    //recenter
    for(int tt=0;tt<proj_data_inrange_b_rc.size();tt++){
      double fval = fitfunc.value(proj_data_inrange_b_rc.coord(tt), params.best());
      double dval = proj_data_inrange_b[q].value(tt).best();

      proj_data_inrange_b_rc.value(tt) = proj_data_inrange_b_rc.value(tt) + (fval - dval);
    }
    bootstrapDistribution<Params> params_rc(one.getInitializer(),guess);
    bootstrapDistributionD chisq_rc(one), chisq_per_dof_rc(one);

    fitter.fit(params_rc, chisq_rc, chisq_per_dof_rc, dof, proj_data_inrange_b_rc);
    
    //Note we shouldn't use the "best" value here because it is not obtained using a bootstrap resampling; instead use the mean
    std::cout << "Bootstrap distribution of q^2 has mean " << chisq_rc.mean() << " and std.dev " << chisq_rc.standardDeviation() << std::endl;

    std::vector<double> q2_dist(nboot); for(int i=0;i<nboot;i++) q2_dist[i] = chisq_rc.sample(i);
    std::sort(q2_dist.begin(), q2_dist.end(), [&](const double a, const double b){ return a<b; });

    double p_boot = computePvalue(chisq.best(), q2_dist);
    std::cout << "Bootstrap p-value " << p_boot << std::endl;  

    //Save sorted bootstrap q2 distribution as rawDataDistribution
    rawDataDistributionD q2_dist_r(nboot, [&](const int i){ return q2_dist[i]; });

    p_boot_all[q] = p_boot;
    q2_dist_r_all[q] = std::move(q2_dist_r);
    chisq_all[q] = std::move(chisq);
    chisq_per_dof_all[q] = std::move(chisq_per_dof);
    params_all[q] = std::move(params);
  }
  
#ifdef HAVE_HDF5
  writeParamsStandard(params_all, "params.hdf5");
  writeParamsStandard(chisq_all, "chisq.hdf5");
  writeParamsStandard(chisq_per_dof_all, "chisq_per_dof.hdf5");
  writeParamsStandard(q2_dist_r_all, "boot_q2.hdf5");
#endif
  {
    std::ofstream of("p_boot.dat");
    for(int q=0;q<nq;q++) of << p_boot_all[q] << std::endl;
  }

  extractMdata<FitKtoPiPi, bootstrapDistribution<Params> > extract(params_all);
  plotErrorWeightedData(proj_data_b, extract, args.tmin_k_op, args.tmin_op_snk);

  std::cout << "Done" << std::endl;

  return 0;
}
