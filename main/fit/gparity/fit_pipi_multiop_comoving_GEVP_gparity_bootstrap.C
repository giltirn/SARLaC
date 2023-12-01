#include <utils.h>
#include <pipi_common/pipi_common.h>
using namespace CPSfit;

#include<fit_pipi_gnd_exc_sim_gparity/enums.h>
#include<fit_pipi_multiop_comoving_gparity/enums.h>
#include<fit_pipi_multiop_comoving_gparity/raw_data.h>
#include<fit_pipi_multiop_comoving_gparity/resampled_data.h>

#include<fit_pipi_multiop_comoving_GEVP_gparity/corr_subtract.h>
#include<fit_pipi_multiop_comoving_GEVP_gparity/resampled_data.h>
#include<fit_pipi_multiop_comoving_GEVP_gparity/analyze.h>

#include<fit_pipi_multiop_comoving_GEVP_gparity_bootstrap/args.h>
#include<fit_pipi_multiop_comoving_GEVP_gparity_bootstrap/cmdline.h>

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
  const int nop = ops.size();
  
  //Read data
  std::map<threeMomentum, RawData> raw_data;
  if(!cmdline.load_combined_data){
    if(cmdline.load_raw_data){
      std::cout << "Reading raw data from " << cmdline.load_raw_data_file << std::endl;
      HDF5reader rd(cmdline.load_raw_data_file);  read(rd, raw_data, "raw_data");
      
      for(int p=0;p<args.p_tot.size();p++){
	auto it = raw_data.find(args.p_tot[p]);
	assert(it!=raw_data.end());
	for(int i=0;i<nop;i++)
	  for(int j=i;j<nop;j++)
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

  //Load/generate resampled data
  std::map<threeMomentum, ResampledData<bootstrapCorrelationFunctionD> > data_j;
  std::map<threeMomentum, ResampledData<bootJackknifeCorrelationFunctionD> > data_dj;

  if(cmdline.load_combined_data){
    loadCheckpoint(data_j,data_dj, cmdline.load_combined_data_file);

    for(int p=0;p<args.p_tot.size();p++){
      auto it = data_j.find(args.p_tot[p]);
      assert(it!=data_j.end());
      for(int i=0;i<nop;i++)
	for(int j=i;j<nop;j++)
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
      data_j[ptot].generatedResampledData(ops,raw_data[ptot], resampler, args.isospin, args.Lt, args.tsep_pipi, ptot, args.do_vacuum_subtraction);
      data_dj[ptot].generatedResampledData(ops,raw_data[ptot], resampler, args.isospin, args.Lt, args.tsep_pipi, ptot, args.do_vacuum_subtraction);
    }
  } 

  if(cmdline.save_combined_data) saveCheckpoint(data_j, data_dj, cmdline.save_combined_data_file);

  //Put resampled data into matrices
  correlationFunction<double, NumericSquareMatrix<bootstrapDistributionD> > C = createCorrelatorMatrix(data_j,args.Lt,ops,args.p_tot);
  correlationFunction<double, NumericSquareMatrix<bootJackknifeDistributionD> > Cdj = createCorrelatorMatrix(data_dj,args.Lt,ops,args.p_tot);

  //Perform transformations on data as desired
  if(cmdline.subtract_nbr_tslice){
    correlatorSubtractNeighbor(C);
    correlatorSubtractNeighbor(Cdj);
  }

  if(cmdline.fix_t_sub){
    correlatorSubtractFixedT(C, cmdline.fix_t_sub_time);
    correlatorSubtractFixedT(Cdj, cmdline.fix_t_sub_time);
  }

  //We don't need a derived solver as we use the eigenvalues/eigenvectors directly, not the effective energies
  GEVPsolverBase<bootstrapDistributionD> gevp(cmdline.verbose_solver);
  GEVPsolverBase<bootJackknifeDistributionD> gevp_dj(cmdline.verbose_solver);

  if(cmdline.load_gevp_solutions){
    std::cout << "Reading GEVP solutions from file "<< cmdline.load_gevp_solutions_file << std::endl;
    HDF5reader rd(cmdline.load_gevp_solutions_file);
    std::cout << "Reading bootstrap GEVP solutions" << std::endl;
    gevp.read(rd, "gevp_solutions");
    std::cout << "Reading boot-jackknife GEVP solutions" << std::endl;
    gevp_dj.read(rd, "gevp_dj_solutions");
  }else{
    std::cout << "Solving bootstrap GEVP" << std::endl;
    gevp.solve(C, args.t_max, args.fit_t0, args.fit_t0);

    std::cout << "Solving boot-jackknife GEVP" << std::endl;   
    gevp_dj.solve(Cdj, args.t_max, args.fit_t0, args.fit_t0);
  }

  if(cmdline.save_gevp_solutions){
    std::cout << "Saving GEVP solutions to file " << cmdline.save_gevp_solutions_file << std::endl;
    HDF5writer wr(cmdline.save_gevp_solutions_file);
    std::cout << "Writing bootstrap GEVP solutions" << std::endl;
    gevp.write(wr, "gevp_solutions");
    std::cout << "Writing boot-jackknife GEVP solutions" << std::endl;
    gevp_dj.write(wr, "gevp_dj_solutions");
  }

  std::vector<bootstrapDistributionD> energies;
  std::vector<bootstrapDistributionD> chisq;
  std::vector<int> dof;
  
  std::cout << "Fitting data" << std::endl;
  fitEigenvaluesNExpFixedt0(energies, chisq, dof, gevp, gevp_dj, nop,
			    args.fit_t_min, args.fit_t_max, args.fit_t0, mlp);

#ifdef HAVE_HDF5
  writeParamsStandard(energies, "energies.hdf5");
  writeParamsStandard(chisq, "chisq.hdf5");
#endif

  //Recenter
  GEVPsolverBase<bootstrapDistributionD> gevp_recentered(gevp);
  FitNStateExp fitfunc(1);
  for(int t=args.fit_t_min; t<=args.fit_t_max; t++){
    auto * evals = gevp_recentered.evals(args.fit_t0, t);
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
  fitEigenvaluesNExpFixedt0(tmp, chisq_recentered, tmp_i, gevp_recentered, gevp_dj, nop,
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
