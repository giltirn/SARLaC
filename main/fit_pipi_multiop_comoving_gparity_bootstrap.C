#include <utils.h>
#include <random.h>
#include <pipi_common/pipi_common.h>
using namespace CPSfit;

#include<fit_pipi_gnd_exc_sim_gparity/enums.h>
#include<fit_pipi_gnd_exc_sim_gparity/fit.h>

#include<fit_pipi_multiop_comoving_gparity/enums.h>
#include<fit_pipi_multiop_comoving_gparity/fitfunc.h>
#include<fit_pipi_multiop_comoving_gparity/raw_data.h>
#include<fit_pipi_multiop_comoving_gparity/resampled_data.h>

#include<fit_pipi_gnd_exc_sigma_sim_gparity_bootstrap/fit.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity_bootstrap/bootstrap_pvalue.h>

#include<fit_pipi_multiop_comoving_gparity_bootstrap/args.h>
#include<fit_pipi_multiop_comoving_gparity_bootstrap/cmdline.h>


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

  //Setup the fit func parameter maps
  std::map< std::pair<Operator,Operator>, SubFitFuncParameterMap > subfit_pmaps;
  ParamTagIdxMap param_map;
  Params guess;
  setupParameterMaps(subfit_pmaps, param_map, guess, args.operators, args.fitfunc, args.Ascale, args.nstate);

  //Write guess template and exit if requested
  if(cmdline.save_guess_template){ saveGuess(guess, "guess_template.dat"); exit(0); }
  
  const std::vector<Operator> &ops = args.operators;

  //Read data
  std::map<threeMomentum, RawData> raw_data;
  if(!cmdline.load_combined_data){
    if(cmdline.load_raw_data){
      std::cout << "Reading raw data from " << cmdline.load_raw_data_file << std::endl;
      HDF5reader rd(cmdline.load_raw_data_file);  read(rd, raw_data, "raw_data");
      
      for(int p=0;p<args.p_tot.size();p++){
	auto it = raw_data.find(args.p_tot[p]);
	assert(it!=raw_data.end());
	for(int i=0;i<ops.size();i++)
	  for(int j=i;j<ops.size();j++)
	    assert(it->second.haveData(ops[i],ops[j]));
      }
    }else{
      for(int p=0;p<args.p_tot.size();p++){
	RawData &raw = raw_data[args.p_tot[p]];

	raw.read(args.isospin, args.Lt, args.data_dir, args.traj_start, args.traj_inc, args.traj_lessthan,
		 args.pipi_figure_file_format, args.pipi_bubble_file_format, args.tsep_pipi, args.tstep_pipi,
		 args.p_tot[p],
		 ops, cmdline.filemap_allow_ptot_parity);
      }
    }
    if(cmdline.save_raw_data){
      std::cout << "Saving raw data to " << cmdline.save_raw_data_file << std::endl;
      HDF5writer wr(cmdline.save_raw_data_file);  write(wr, raw_data, "raw_data");
    }
  }

      
    
  std::map<threeMomentum, ResampledData<bootstrapCorrelationFunctionD> > data_b;
  std::map<threeMomentum, ResampledData<bootJackknifeCorrelationFunctionD> > data_bj;
  if(cmdline.load_combined_data){
    loadCheckpoint(data_b, data_bj, cmdline.load_combined_data_file);
    for(int p=0;p<args.p_tot.size();p++){
      auto it = data_b.find(args.p_tot[p]);
      assert(it!=data_b.end());
      for(int i=0;i<ops.size();i++)
	for(int j=i;j<ops.size();j++)
	  assert(it->second.haveData(ops[i],ops[j]));
    }    

  }else{
    int nsample = raw_data[args.p_tot[0]].correlator(ops[0],ops[0]).value(0).size();

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

    for(int p=0;p<args.p_tot.size();p++){
      auto ptot = args.p_tot[p];
      data_b[ptot].generatedResampledData(ops,raw_data[ptot], resampler, args.isospin, args.Lt, args.tsep_pipi, ptot, args.do_vacuum_subtraction);
      data_bj[ptot].generatedResampledData(ops,raw_data[ptot], resampler, args.isospin, args.Lt, args.tsep_pipi, ptot, args.do_vacuum_subtraction);
    }
  }  
  
  if(cmdline.save_combined_data) saveCheckpoint(data_b, data_bj, cmdline.save_combined_data_file);

  //Add resampled data to full data set with generalized coordinate set appropriately
  correlationFunction<SimFitCoordGen,  bootstrapDistributionD> corr_comb_b;
  correlationFunction<SimFitCoordGen,  bootJackknifeDistributionD> corr_comb_bj;
  
  std::vector<std::string> vnm; //for printing

  std::map<std::unordered_map<std::string, std::string> const*, std::string> pmap_descr;
  for(int i=0;i<ops.size();i++){
    for(int j=i;j<ops.size();j++){
      std::ostringstream nm; nm << ops[i] << " " << ops[j];
      pmap_descr[&subfit_pmaps[{ops[i],ops[j]}]] = nm.str();

      for(int t=args.t_min;t<=args.t_max;t++){
	SimFitCoordGen coord(t, &subfit_pmaps[{ops[i],ops[j]}] , 2*args.tsep_pipi);

	//Average over the total momentum values (it is assumed these are equivalent)
	bootstrapDistributionD value_b = data_b[args.p_tot[0]].correlator(ops[i],ops[j]).value(t);
	bootJackknifeDistributionD value_bj = data_bj[args.p_tot[0]].correlator(ops[i],ops[j]).value(t);
	for(int p=1;p<args.p_tot.size();p++){
	  auto ptot = args.p_tot[p];
	  value_b = value_b + data_b[ptot].correlator(ops[i],ops[j]).value(t);
	  value_bj = value_bj + data_bj[ptot].correlator(ops[i],ops[j]).value(t);
	}
	value_b = value_b / double(args.p_tot.size());
	value_bj = value_bj / double(args.p_tot.size());
	
	corr_comb_b.push_back(coord, value_b);
	corr_comb_bj.push_back(coord, value_bj);
	vnm.push_back(nm.str());
      }
    }
  }

  std::cout << "Data in fit:" << std::endl;
  for(int i=0;i<corr_comb_b.size();i++){
    std::cout << vnm[i] << " " << corr_comb_b.coord(i).t << " " << corr_comb_b.value(i) << std::endl;
  }
  if(cmdline.print_data_sample){
    std::cout << "Data sample " << cmdline.print_data_sample_idx << std::endl;
    for(int i=0;i<corr_comb_b.size();i++){
      double val = cmdline.print_data_sample_idx == -1 ? corr_comb_b.value(i).best()  : corr_comb_b.value(i).sample(cmdline.print_data_sample_idx);
      printf("%s %d %.16e\n", vnm[i].c_str(), int(corr_comb_b.coord(i).t), val);
    }
  }

  std::cout << "Performing any data transformations required by the fit func" << std::endl;
  transformData(corr_comb_b, args.t_min, args.t_max, args.fitfunc);
  transformData(corr_comb_bj, args.t_min, args.t_max, args.fitfunc);


  std::cout << "Data post-transformation:" << std::endl;
  for(int i=0;i<corr_comb_b.size();i++){
    std::cout << pmap_descr[corr_comb_b.coord(i).param_map] << " " << corr_comb_b.coord(i).t << " " << corr_comb_b.value(i) << std::endl;
  }

  if(cmdline.load_guess) loadGuess(guess, cmdline.guess_file);
    
  std::cout << "Performing fit" << std::endl;
  bootstrapInitType binit(args.nboot);
  bootstrapDistribution<Params> params(binit, guess);
  bootstrapDistributionD chisq(binit), chisq_per_dof(binit);
  
  CovarianceMatrix cov_mat = CovarianceMatrix::Regular;
  
  fitOptions opt;
  cmdline.exportOptions(opt);

  fit(params, chisq, chisq_per_dof,
      corr_comb_b, corr_comb_bj, args.fitfunc, param_map,
      args.nstate, args.Lt, args.t_min, args.t_max, args.correlated, cov_mat, args.Ascale, args.Cscale, opt);

  std::cout << "Params:\n";
  {
    for(int p=0;p<params.sample(0).size();p++){
      std::string tag = params.sample(0).tag(p);
      bootstrapDistributionD tmp;
      standardIOhelper<bootstrapDistributionD, bootstrapDistribution<Params> >::extractStructEntry(tmp, params, p);
      std::cout << tag << " = " << tmp << std::endl;
    }
  }

  double dof = chisq.best()/chisq_per_dof.best();  
  bootstrapDistributionD pvalue_chisq(binit);
  for(int s=0;s<iterate<bootstrapDistributionD>::size(pvalue_chisq);s++) 
    iterate<bootstrapDistributionD>::at(s,pvalue_chisq) = chiSquareDistribution::pvalue(dof, iterate<bootstrapDistributionD>::at(s, chisq));

  std::cout << std::endl;
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
  std::cout << "P-value(chi^2): " << pvalue_chisq << std::endl;

#ifdef HAVE_HDF5
  writeParamsStandard(params, "params.hdf5");
  writeParamsStandard(chisq, "chisq.hdf5");
  writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");
  writeParamsStandard(pvalue_chisq, "pvalue_chisq.hdf5");
#endif

  bootstrapPvalue(chisq.best(), 
		  corr_comb_b, corr_comb_bj, params,
		  args.fitfunc, param_map,
		  args.nstate, args.Lt, args.t_min, args.t_max, 
		  args.correlated, cov_mat, args.Ascale, args.Cscale, opt, true);

  std::cout << "Done\n";
  return 0;
}
