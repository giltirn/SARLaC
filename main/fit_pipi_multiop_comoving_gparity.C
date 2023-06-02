#include <utils.h>
#include <random.h>
#include <pipi_common/pipi_common.h>
using namespace CPSfit;

#include<fit_pipi_gnd_exc_sim_gparity/enums.h>
#include<fit_pipi_gnd_exc_sim_gparity/fit.h>

#include<fit_pipi_multiop_comoving_gparity/args.h>
#include<fit_pipi_multiop_comoving_gparity/cmdline.h>
#include<fit_pipi_multiop_comoving_gparity/fitfunc.h>
#include<fit_pipi_multiop_comoving_gparity/raw_data.h>
#include<fit_pipi_multiop_comoving_gparity/resampled_data.h>

#include<fit_pipi_gnd_exc_sigma_sim_gparity/fit.h>

int main(const int argc, const char* argv[]){
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

      
    
  std::map<threeMomentum, ResampledData<jackknifeCorrelationFunctionD> > data_j;
  std::map<threeMomentum, ResampledData<doubleJackknifeCorrelationFunctionD> > data_dj;
  if(cmdline.load_combined_data){
    loadCheckpoint(data_j, data_dj, cmdline.load_combined_data_file);
    for(int p=0;p<args.p_tot.size();p++){
      auto it = data_j.find(args.p_tot[p]);
      assert(it!=data_j.end());
      for(int i=0;i<ops.size();i++)
	for(int j=i;j<ops.size();j++)
	  assert(it->second.haveData(ops[i],ops[j]));
    }    

  }else{
    for(int p=0;p<args.p_tot.size();p++){
      auto ptot = args.p_tot[p];
      data_j[ptot].generatedResampledData(ops,raw_data[ptot], args.bin_size, args.isospin, args.Lt, args.tsep_pipi, ptot, args.do_vacuum_subtraction);
      data_dj[ptot].generatedResampledData(ops,raw_data[ptot], args.bin_size, args.isospin, args.Lt, args.tsep_pipi, ptot, args.do_vacuum_subtraction);
    }
  }  
  
  if(cmdline.save_combined_data) saveCheckpoint(data_j, data_dj, cmdline.save_combined_data_file);

  int nsample = (args.traj_lessthan - args.traj_start)/args.traj_inc/args.bin_size;
  
  //Add resampled data to full data set with generalized coordinate set appropriately
  correlationFunction<SimFitCoordGen,  jackknifeDistributionD> corr_comb_j;
  correlationFunction<SimFitCoordGen,  doubleJackknifeDistributionD> corr_comb_dj;
  
  std::vector<std::string> vnm; //for printing

  std::map<std::unordered_map<std::string, std::string> const*, std::string> pmap_descr;
  for(int i=0;i<ops.size();i++){
    for(int j=i;j<ops.size();j++){
      std::ostringstream nm; nm << ops[i] << " " << ops[j];
      pmap_descr[&subfit_pmaps[{ops[i],ops[j]}]] = nm.str();

      for(int t=args.t_min;t<=args.t_max;t++){
	SimFitCoordGen coord(t, &subfit_pmaps[{ops[i],ops[j]}] , 2*args.tsep_pipi);

	//Average over the total momentum values (it is assumed these are equivalent)
	jackknifeDistributionD value_j(nsample,0.);
	doubleJackknifeDistributionD value_dj(nsample,0.);
	for(int p=0;p<args.p_tot.size();p++){
	  auto ptot = args.p_tot[p];
	  value_j = value_j + data_j[ptot].correlator(ops[i],ops[j]).value(t);
	  value_dj = value_dj + data_dj[ptot].correlator(ops[i],ops[j]).value(t);
	}
	value_j = value_j / double(args.p_tot.size());
	value_dj = value_dj / double(args.p_tot.size());
	
	corr_comb_j.push_back(coord, value_j);
	corr_comb_dj.push_back(coord, value_dj);
	vnm.push_back(nm.str());
      }
    }
  }

  std::cout << "Data in fit:" << std::endl;
  for(int i=0;i<corr_comb_j.size();i++){
    std::cout << vnm[i] << " " << corr_comb_j.coord(i).t << " " << corr_comb_j.value(i) << std::endl;
  }
  
  std::cout << "Performing any data transformations required by the fit func" << std::endl;
  transformData(corr_comb_j, args.t_min, args.t_max, args.fitfunc);
  transformData(corr_comb_dj, args.t_min, args.t_max, args.fitfunc);


  std::cout << "Data post-transformation:" << std::endl;
  for(int i=0;i<corr_comb_j.size();i++){
    std::cout << pmap_descr[corr_comb_j.coord(i).param_map] << " " << corr_comb_j.coord(i).t << " " << corr_comb_j.value(i) << std::endl;
  }

  if(cmdline.load_guess) loadGuess(guess, cmdline.guess_file);
    
  std::cout << "Performing fit" << std::endl;
  
  jackknifeDistribution<Params> params(nsample, guess);
  jackknifeDistributionD chisq(nsample), chisq_per_dof(nsample);
  
  //Does not support block double jackknife (this technique gives more reliable covariance matrices when binning)
  CovarianceMatrix cov_mat = CovarianceMatrix::Regular;
  correlationFunction<SimFitCoordGen,  blockDoubleJackknifeDistributionD> not_used;
  
  fitOptions opt;
  cmdline.exportOptions(opt);

  fit(params, chisq, chisq_per_dof,
      corr_comb_j, corr_comb_dj, not_used, args.fitfunc, param_map,
      args.nstate, args.Lt, args.t_min, args.t_max, args.correlated, cov_mat, args.Ascale, args.Cscale, opt);

  std::cout << "Params:\n";
  {
    for(int p=0;p<params.sample(0).size();p++){
      std::string tag = params.sample(0).tag(p);
      jackknifeDistributionD tmp;
      standardIOhelper<jackknifeDistributionD, jackknifeDistribution<Params> >::extractStructEntry(tmp, params, p);
      std::cout << tag << " = " << tmp << std::endl;
    }
  }

  double dof = chisq.sample(0)/chisq_per_dof.sample(0);  
  jackknifeDistributionD pvalue(nsample, [&](const int s){ return chiSquareDistribution::pvalue(dof, chisq.sample(s)); });

  std::cout << std::endl;
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
  std::cout << "P-value: " << pvalue << std::endl;

#ifdef HAVE_HDF5
  writeParamsStandard(params, "params.hdf5");
  writeParamsStandard(chisq, "chisq.hdf5");
  writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");
  writeParamsStandard(pvalue, "pvalue.hdf5");
#endif

  if(args.correlated) analyzeChisq(corr_comb_j, params, args.fitfunc, args.nstate, args.Lt, args.t_min, args.t_max, args.Ascale, args.Cscale, pmap_descr); 

  std::cout << "Done\n";
  return 0;
}
