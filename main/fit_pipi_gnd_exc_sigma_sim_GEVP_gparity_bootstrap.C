#include <utils.h>
#include <random.h>
#include <pipi_common/pipi_common.h>
using namespace CPSfit;

#include<fit_pipi_gnd_exc_sim_gparity/enums.h>
#include<fit_pipi_gnd_exc_sim_gparity/fit.h>

#include<fit_pipi_gnd_exc_sigma_sim_gparity/fitfunc.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/filters.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/raw_data.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/resampled_data.h>

#include<fit_pipi_multiop_comoving_GEVP_gparity/analyze.h>

#include<fit_pipi_gnd_exc_sigma_sim_GEVP_gparity_bootstrap/args.h>
#include<fit_pipi_gnd_exc_sigma_sim_GEVP_gparity_bootstrap/cmdline.h>


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

  MarquardtLevenbergParameters<double> mlp;
  mlp.verbose = true;
  if(cmdline.load_minimizer_params ){
    parse(mlp, cmdline.minimizer_params_file);
  }
  
  const std::vector<Operator> &ops = args.operators;
  int nop = ops.size();

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
		    ops, args.isospin);
    }
    if(cmdline.save_raw_data){
      std::cout << "Saving raw data to " << cmdline.save_raw_data_file << std::endl;
      HDF5writer wr(cmdline.save_raw_data_file);  raw_data.write(wr, "raw_data");
    }

    if(cmdline.remove_samples_in_range)
      raw_data.removeSamplesInRange(cmdline.remove_samples_in_range_start, cmdline.remove_samples_in_range_lessthan);
    if(cmdline.scramble_raw_data)
      raw_data.scrambleSamples();
  }

  ResampledData<bootstrapCorrelationFunctionD> data_b;
  ResampledData<bootJackknifeCorrelationFunctionD> data_bj;
  if(cmdline.load_combined_data){
    std::cout << "Reading previously-generated resampled data" << std::endl;
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

    if(args.isospin != 0 && args.do_vacuum_subtraction) error_exit(std::cout << "Vacuum subtraction is not appropriate to I!=0\n");
    bootstrapBlockResampler resampler(rtable);
    data_b.generatedResampledData(raw_data, resampler, args.Lt, args.tsep_pipi, args.do_vacuum_subtraction, args.timeslice_avg_vac_sub);
    data_bj.generatedResampledData(raw_data, resampler, args.Lt, args.tsep_pipi, args.do_vacuum_subtraction, args.timeslice_avg_vac_sub);
  }  

  if(cmdline.save_combined_data){
    HDF5writer wr(cmdline.save_combined_data_file);
    write(wr, data_b, "data_b");
    write(wr, data_bj, "data_bj");
  }

  std::cout << "Performing GEVP" << std::endl;

  //Put resampled data into matrices
  correlationFunction<double, NumericSquareMatrix<bootstrapDistributionD> > Cb(args.Lt);
  correlationFunction<double, NumericSquareMatrix<bootJackknifeDistributionD> > Cbj(args.Lt);
  for(int t=0;t<args.Lt;t++){
    Cb.coord(t) = Cbj.coord(t) = t;
    Cb.value(t).resize(nop);
    Cbj.value(t).resize(nop);

    for(int i=0;i<nop;i++)
      for(int j=i;j<nop;j++){
	Cb.value(t)(i,j) = Cb.value(t)(j,i) = data_b.correlator(ops[i],ops[j]).value(t);
	Cbj.value(t)(i,j) = Cbj.value(t)(j,i) = data_bj.correlator(ops[i],ops[j]).value(t);
      }
  }
  GEVPsolver<bootstrapDistributionD> gevpb(true);
  GEVPsolver<bootJackknifeDistributionD> gevpbj(true);

  if(cmdline.load_gevp_solutions){
    std::cout << "Reading GEVP solutions from file "<< cmdline.load_gevp_solutions_file << std::endl;
    HDF5reader rd(cmdline.load_gevp_solutions_file);
    std::cout << "Reading bootstrap GEVP solutions" << std::endl;
    gevpb.read(rd, "gevpb_solutions");
    std::cout << "Reading boot-jackknife GEVP solutions" << std::endl;
    gevpbj.read(rd, "gevpj_solutions");
  }else{
    std::cout << "Solving bootstrap GEVP" << std::endl;
    gevpb.solve(Cb, args.t_max, args.fit_t0, args.fit_t0);

    std::cout << "Solving boot-jackknife GEVP" << std::endl;   
    gevpbj.solve(Cbj, args.t_max, args.fit_t0, args.fit_t0);
  }

  if(cmdline.save_gevp_solutions){
    std::cout << "Saving GEVP solutions to file " << cmdline.save_gevp_solutions_file << std::endl;
    HDF5writer wr(cmdline.save_gevp_solutions_file);
    std::cout << "Writing bootstrap GEVP solutions" << std::endl;
    gevpb.write(wr, "gevpb_solutions");
    std::cout << "Writing boot-jackknife GEVP solutions" << std::endl;
    gevpbj.write(wr, "gevpj_solutions");
  }

  //Write the data out for analysis
  {
    std::vector<std::vector<bootstrapDistributionD> > data(nop, std::vector<bootstrapDistributionD>(args.fit_t_max-args.fit_t_min+1));
    std::vector<std::vector<bootstrapDistributionD> > resids(nop, std::vector<bootstrapDistributionD>(args.fit_t_max-args.fit_t_min+1));
    for(int t=args.fit_t_min; t<=args.fit_t_max; t++){
      auto const* evals_j = gevpb.evals(args.fit_t0, t);
      auto const* resids_j = gevpb.residuals(args.fit_t0, t);
      if(evals_j != NULL){
	for(int o=0;o<nop;o++){
	  data[o][t-args.fit_t_min] = (*evals_j)[o];
	  resids[o][t-args.fit_t_min] = (*resids_j)[o];
	}
      }else{
	for(int o=0;o<nop;o++){
	  zeroit(data[o][t-args.fit_t_min]);
	  zeroit(resids[o][t-args.fit_t_min]);
	}
      }
    }
    writeParamsStandard(data, "data.hdf5");    
    writeParamsStandard(resids, "data_resids.hdf5");    
  }


  std::vector<bootstrapDistributionD> energies;
  std::vector<bootstrapDistributionD> chisq;
  std::vector<int> dof;
  
  std::cout << "Fitting data" << std::endl;
  fitEigenvaluesNExpFixedt0(energies, chisq, dof, gevpb, gevpbj, nop,
			    args.fit_t_min, args.fit_t_max, args.fit_t0, mlp);

#ifdef HAVE_HDF5
  writeParamsStandard(energies, "energies.hdf5");
  writeParamsStandard(chisq, "chisq.hdf5");
#endif

  //Recenter
  GEVPsolver<bootstrapDistributionD> gevpb_recentered(gevpb);
  FitNStateExp fitfunc(1);
  for(int t=args.fit_t_min; t<=args.fit_t_max; t++){
    auto * evals = gevpb_recentered.evals(args.fit_t0, t);
    assert(evals != NULL);
    for(int n=0;n<nop;n++){
      parameterVector<double> pcen({1., energies[n].best()});
      double correction = fitfunc.value(t-args.fit_t0, pcen) - (*evals)[n].best();
      typedef iterate<bootstrapDistributionD> it_t;
      for(int s=0;s<it_t::size((*evals)[n]);s++)
	it_t::at(s, (*evals)[n]) += correction;
    }
  }
  std::vector<bootstrapDistributionD> tmp;
  std::vector<bootstrapDistributionD> chisq_recentered;
  std::vector<int> tmp_i;

  std::cout << "Fitting recentered data" << std::endl;
  fitEigenvaluesNExpFixedt0(tmp, chisq_recentered, tmp_i, gevpb_recentered, gevpbj, nop,
			    args.fit_t_min, args.fit_t_max, args.fit_t0, mlp);     
  
  std::vector<double> p_boot(nop);
  for(int n=0;n<nop;n++){
    std::cout << "Computing bootstrap p-value for state " << n << std::endl;
    std::cout << "Bootstrap distribution of q^2 has mean " << chisq_recentered[n].mean() << " and std.dev " << chisq_recentered[n].standardDeviation() << std::endl;
    p_boot[n] = computePvalue(chisq[n].best(), chisq_recentered[n]);
    std::cout << "pvalue=" << p_boot[n] << std::endl;
  }
  
  {
    std::ofstream of("p_boot.dat");
    for(int n=0;n<nop;n++)
      of << n << " " << p_boot[n] << std::endl;
  }

  std::cout << "Done\n";
  return 0;
}


