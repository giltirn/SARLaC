#include<ktopipi_common/ktopipi_common.h>

#include<fit_ktopipi_ktosigma_combined_gparity/simfit_generic.h>
#include<fit_ktopipi_gnd_exc_ktosigma_combined_gparity/resampled_data.h>
#include<fit_ktopipi_gnd_exc_ktosigma_combined_gparity/simfit_plot.h>
#include<fit_ktopipi_gnd_exc_ktosigma_combined_gparity/simfit.h>
#include<fit_ktopipi_gnd_exc_ktosigma_combined_gparity/fronthalf_backhalf.h>
#include<fit_ktopipi_gnd_exc_ktosigma_combined_gparity/weighted_avg_consistency.h>
#include<fit_ktopipi_gnd_exc_ktosigma_combined_gparity/plot_C_fixedtsepKop.h>

#include<fit_ktopipi_gnd_exc_ktosigma_combined_gparity_bootstrap/alpha_vary_plot.h>
#include<fit_ktopipi_gnd_exc_ktosigma_combined_gparity_bootstrap/args.h>
#include<fit_ktopipi_gnd_exc_ktosigma_combined_gparity_bootstrap/cmdline.h>
#include<fit_ktopipi_gnd_exc_ktosigma_combined_gparity_bootstrap/bootstrap_pvalue.h>


using namespace CPSfit;

//Make sure you delete the rtable when you are done
bootstrapBlockResampler getResampler(const RawData &raw, const Args &args, const CMDline &cmdline){
  int nsample = raw.nsample();

  //Generate the resample table
  resampleTableOptions ropt;
  ropt.read_from_file = cmdline.load_boot_resample_table;
  ropt.read_file = cmdline.load_boot_resample_table_file;
  ropt.write_to_file = cmdline.save_boot_resample_table;
  ropt.write_file = cmdline.save_boot_resample_table_file;

  std::vector<std::vector<int> >* rtable = new std::vector<std::vector<int> >(generateResampleTable(nsample, args.nboot, args.resample_table_type, args.block_size, RNG, ropt));
  if(rtable->size() != args.nboot) error_exit(std::cout << "Expected resample table of size " << args.nboot << ", got " << rtable->size() << std::endl);

  bootstrapBlockResampler resampler(*rtable);
  return resampler;
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

  simultaneousFitBase<bootstrapDistribution>* fitter = getFitter<bootstrapDistribution>(args.fitfunc, args.nstate, args.operators);

  fitter->load2ptFitParams(args.operators, args.input_params, args.nboot); 

  if(cmdline.load_frozen_fit_params) //Freeze matrix elements (not amplitudes/energies!)
    fitter->loadFrozenMatrixElements(cmdline.load_frozen_fit_params_file, args.nboot);
  
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

    if(cmdline.fronthalf_backhalf){
      if(!RNG.isInitialized()) RNG.initialize(1234);
      int nsample = raw.nsample();
      std::vector<std::vector<int> > rtable_fh = generateResampleTable(nsample/2, args.nboot, args.resample_table_type, args.block_size, RNG);
      std::vector<std::vector<int> > rtable_bh = generateResampleTable(nsample - nsample/2, args.nboot, args.resample_table_type, args.block_size, RNG);
	
      frontHalfBackHalfAnalysis<bootstrapDistributionD>(raw, bootstrapBlockResampler(rtable_fh), bootstrapBlockResampler(rtable_bh), args, cmdline);
      exit(0);
    }

    if(cmdline.remove_samples_in_range)
      raw.removeSamplesInRange(cmdline.remove_samples_in_range_start, cmdline.remove_samples_in_range_lessthan);      

    if(cmdline.write_alpha_and_pseudoscalar_matrix_elem){
      bootstrapBlockResampler resampler = getResampler(raw, args, cmdline);
      std::vector<std::vector<std::vector<std::vector<bootstrapDistributionD> > > > alpha(args.operators.size()); //[op_idx][tsep_k_snk_idx][q][t]
      std::vector<std::vector<correlationFunction<double, bootstrapDistributionD> > > op_P_K(args.operators.size()); //[op_idx][tsep_k_snk_idx][t]
      for(int o=0;o<args.operators.size();o++)
	raw.computeAlphaAndPseudoscalarMatrixElem(alpha[o], op_P_K[o], args.operators[o], args, cmdline, "alpha/op_P_K write", resampler);
      {
	HDF5writer wr("alpha.hdf5");
	write(wr, alpha, "alpha");
      }
      {
	HDF5writer wr("op_P_K.hdf5");
	write(wr, op_P_K, "op_P_K");
      }
      delete &resampler.rtable;
    }

    if(cmdline.alpha_vary_plot){      
      AlphaVaryPlotArgs pargs;
      parse(pargs, cmdline.alpha_vary_plot_args);

      bootstrapBlockResampler resampler = getResampler(raw, args, cmdline);
      for(int q=0;q<10;q++){
	for(int o=0;o<args.operators.size();o++){
	  PiPiOperator op = args.operators[o];
	  int tsep_k_snk;
	  if(op == PiPiOperator::PiPiGnd || op == PiPiOperator::PiPiExc) tsep_k_snk =  args.tsep_k_pi[pargs.tsep_k_snk_idx];
	  else tsep_k_snk =  args.tsep_k_sigma[pargs.tsep_k_snk_idx];

	  std::string file_stub = stringize("alpha_vary_Q%d_%s_tsepksnk%d_tsepopsnk%d", q+1,
					    opDescrFile(args.operators[o]).c_str(),
					    tsep_k_snk, pargs.tsep_op_snk);
	  alphaVaryPlot<bootstrapDistributionD>(raw, args.operators[o], q, pargs.tsep_k_snk_idx, pargs.tsep_op_snk,		    
						args, cmdline, resampler, file_stub, pargs.nstep_each_side, pargs.coeff_step);
	    
	}
      }
      delete &resampler.rtable;
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
    bootstrapBlockResampler resampler = getResampler(raw, args, cmdline);

    double alpha_coeff = cmdline.set_alpha_coeff ? cmdline.alpha_coeff : 1.;
    std::cout << "Generating resampled data with alpha coefficient " << alpha_coeff << " (default 1.)" << std::endl;

    data_b.resample(raw, args, cmdline, "bootstrap", resampler, alpha_coeff);
    data_bj.resample(raw, args, cmdline, "boot-jackknife", resampler, alpha_coeff);
  }
  
  if(cmdline.save_resampled_data_container_checkpoint){
    std::cout << "Writing resampled data container to checkpoint file" << std::endl;
    HDF5writer wr(cmdline.save_resampled_data_container_checkpoint_file);
    write(wr, data_b, "data_b");
    write(wr, data_bj, "data_bj");
  }

  if(cmdline.prune_data){
    pruneArgs pargs;
    parse(pargs, cmdline.prune_data_file);
    data_b.prune(pargs);
    data_bj.prune(pargs);
  }

  if(args.basis == Basis::Basis7){
    std::cout << "Converting to 7-basis" << std::endl;
    data_b.convertBasis10to7(); 
    data_bj.convertBasis10to7();
  }

  checkWeightedAvgConsistency(data_b, args.input_params, args.operators, args.tmin_k_op);

  //Constrain tsep_k_snk to those in input file (in case we read resampled data that contains more)
  std::map<PiPiOperator, std::vector<int> const*> op_tsep_list = {  {PiPiOperator::PiPiGnd, &args.tsep_k_pi}, {PiPiOperator::PiPiExc, &args.tsep_k_pi}, {PiPiOperator::Sigma, &args.tsep_k_sigma} };

  for(int o=0;o<args.operators.size();o++){
    data_b.constrainSourceSinkSep(args.operators[o], *op_tsep_list[args.operators[o]]);
    data_bj.constrainSourceSinkSep(args.operators[o], *op_tsep_list[args.operators[o]]);
  }

  //Plot correlation function for each operator and tsep_k_pi
  for(int o=0;o<args.operators.size();o++){
    for(int q=0; q < (args.basis == Basis::Basis7 ? 7 : 10); q++){
      plotCfixedTsepKop(data_b, args.operators[o], q+1, *op_tsep_list[args.operators[o]]);
    }
  }

  //Do fits
  std::cout << "Starting fits" << std::endl;
  typedef taggedValueContainer<double,std::string> Params;

  ResampledDataContainers<bootstrapDistribution> rdata(data_b, data_bj);

  std::vector<bootstrapDistribution<Params> > params;
  std::vector<bootstrapDistributionD> chisq;
    
  fitter->fit(params, chisq, rdata, args.operators,
	      args.Lt, args.tmin_k_op, args.tmin_op_snk, args.correlated, args.covariance_matrix);

  if(args.basis == Basis::Basis7){
    std::cout << "Converting 7 basis results to 10 basis" << std::endl;
    std::vector<std::vector<bootstrapDistributionD> > params_10 = convert7basisTo10basis(args.nstate, params);
    writeParamsStandard(params_10, "matrix_elems_10basis.hdf5");
    for(int q=0;q<10;q++){
      std::cout << "Q" << q+1 << std::endl;
      for(int m=0;m<args.nstate;m++)
	std::cout << "M" << m << " = " << params_10[q][m] << std::endl;
    }    
  }

  std::vector<double> q2_best(chisq.size());
  for(int i=0;i<chisq.size();i++) q2_best[i] = chisq[i].best();

  bootstrapPvalue(q2_best, data_b, data_bj, params, fitter, args.operators, args.Lt, args.tmin_k_op, args.tmin_op_snk, args.correlated, args.covariance_matrix);

  std::cout << "Done" << std::endl;
  
  delete fitter;

  return 0;
}
