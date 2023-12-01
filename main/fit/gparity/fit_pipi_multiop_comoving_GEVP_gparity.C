#if 1

int main(void){
  return 0;
}

#else

#include <utils.h>
#include <pipi_common/pipi_common.h>
using namespace CPSfit;

#include<fit_pipi_gnd_exc_sim_gparity/enums.h>
#include<fit_pipi_multiop_comoving_gparity/enums.h>
#include<fit_pipi_multiop_comoving_gparity/raw_data.h>
#include<fit_pipi_multiop_comoving_gparity/resampled_data.h>

#include<fit_pipi_multiop_comoving_GEVP_gparity/args.h>
#include<fit_pipi_multiop_comoving_GEVP_gparity/corr_subtract.h>
#include<fit_pipi_multiop_comoving_GEVP_gparity/cmdline.h>
#include<fit_pipi_multiop_comoving_GEVP_gparity/resampled_data.h>
#include<fit_pipi_multiop_comoving_GEVP_gparity/analyze.h>

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

  int nsample = (args.traj_lessthan - args.traj_start)/args.traj_inc/args.bin_size;      

  //Load/generate resampled data
  std::map<threeMomentum, ResampledData<jackknifeCorrelationFunctionD> > data_j;
  std::map<threeMomentum, ResampledData<doubleJackknifeCorrelationFunctionD> > data_dj;

  if(cmdline.load_combined_data){
    if(args.double_jackknife_cov)
      loadCheckpoint(data_j,data_dj, cmdline.load_combined_data_file);
    else
      loadCheckpoint(data_j, cmdline.load_combined_data_file);

    for(int p=0;p<args.p_tot.size();p++){
      auto it = data_j.find(args.p_tot[p]);
      assert(it!=data_j.end());
      for(int i=0;i<nop;i++)
	for(int j=i;j<nop;j++)
	  assert(it->second.haveData(ops[i],ops[j]));
    }    

  }else{
    for(int p=0;p<args.p_tot.size();p++){
      auto ptot = args.p_tot[p];
      data_j[ptot].generatedResampledData(ops, raw_data[ptot], args.bin_size, args.isospin, args.Lt, args.tsep_pipi, ptot, args.do_vacuum_subtraction);
      if(args.double_jackknife_cov)
	data_dj[ptot].generatedResampledData(ops, raw_data[ptot], args.bin_size, args.isospin, args.Lt, args.tsep_pipi, ptot, args.do_vacuum_subtraction);
    }
  } 

  if(cmdline.save_combined_data) saveCheckpoint(data_j, cmdline.save_combined_data_file);

  //Put resampled data into matrices
  correlationFunction<double, NumericSquareMatrix<jackknifeDistributionD> > C = createCorrelatorMatrix(data_j,args.Lt,ops,args.p_tot);
  correlationFunction<double, NumericSquareMatrix<doubleJackknifeDistributionD> > Cdj;
  if(args.double_jackknife_cov) Cdj = createCorrelatorMatrix(data_dj,args.Lt,ops,args.p_tot);

  //Perform transformations on data as desired
  if(cmdline.subtract_from_data){
    correlatorSubtract(C, cmdline.subtract_from_data_file);
    if(args.double_jackknife_cov) correlatorSubtract(Cdj, cmdline.subtract_from_data_file);
  }

  if(cmdline.subtract_nbr_tslice){
    correlatorSubtractNeighbor(C);
    if(args.double_jackknife_cov) correlatorSubtractNeighbor(Cdj);
  }

  if(cmdline.fix_t_sub){
    correlatorSubtractFixedT(C, cmdline.fix_t_sub_time);
    if(args.double_jackknife_cov) correlatorSubtractFixedT(Cdj, cmdline.fix_t_sub_time);
  }

  //Instantiate the appropriate solver and solve/analyze

  if(cmdline.subtract_nbr_tslice){
    GEVPsubNeighborTslice<jackknifeDistributionD> gevp(cmdline.verbose_solver);
    gevp.solve(C, args.t_max, args.fit_t0min, args.fit_t0max);

    if(args.double_jackknife_cov){
      GEVPsubNeighborTslice<doubleJackknifeDistributionD> gevp_dj(cmdline.verbose_solver);
      gevp_dj.solve(Cdj, args.t_max, args.fit_t0min, args.fit_t0max);
      analyze_GEVP(gevp, gevp_dj, C, args.t_max, args.fit_tmin, args.fit_tmax, args.fit_t0min, args.fit_t0max, args.Ascale);
    }else{
      analyze_GEVP(gevp, C, args.t_max, args.fit_tmin, args.fit_tmax, args.fit_t0min, args.fit_t0max, args.Ascale);
    }
  }else if(cmdline.fix_t_sub){
    if(args.double_jackknife_cov) error_exit(std::cout << "Fixed tsub does not support double-jackknife\n");
    GEVPsubFixedTslice gevp(cmdline.fix_t_sub_time, cmdline.verbose_solver);
    gevp.solve(C, args.t_max, args.fit_t0min, args.fit_t0max);
    analyze_GEVP(gevp, C, args.t_max, args.fit_tmin, args.fit_tmax, args.fit_t0min, args.fit_t0max, args.Ascale);
  }else{
    GEVPsolver<jackknifeDistributionD> gevp(cmdline.verbose_solver);
    gevp.solve(C, args.t_max, args.fit_t0min, args.fit_t0max);

    if(args.double_jackknife_cov){
      GEVPsolver<doubleJackknifeDistributionD> gevp_dj(cmdline.verbose_solver);
      gevp_dj.solve(C, args.t_max, args.fit_t0min, args.fit_t0max);
      analyze_GEVP(gevp, gevp_dj, C, args.t_max, args.fit_tmin, args.fit_tmax, args.fit_t0min, args.fit_t0max, args.Ascale);
    }else{
      analyze_GEVP(gevp, C, args.t_max, args.fit_tmin, args.fit_tmax, args.fit_t0min, args.fit_t0max, args.Ascale);
    }
  }

  std::cout << "Done\n";
  return 0;
}

#endif
