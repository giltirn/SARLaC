#include <utils.h>
#include <random.h>
#include <pipi_common/pipi_common.h>
using namespace CPSfit;

#include<fit_pipi_gnd_exc_sim_gparity/enums.h>

#include<fit_pipi_gnd_exc_sigma_sim_gparity/enums.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/raw_data.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/resampled_data.h>

#include<fit_pipi_multiop_comoving_GEVP_gparity/corr_subtract.h>
#include<fit_pipi_multiop_comoving_GEVP_gparity/analyze.h>

#include<fit_pipi_gnd_exc_sigma_sim_GEVP_gparity/args.h>
#include<fit_pipi_gnd_exc_sigma_sim_GEVP_gparity/cmdline.h>

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
  
  const std::vector<Operator> &ops = args.operators;
  const int nop = ops.size();

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
		    ops);
    }
    if(cmdline.save_raw_data){
      std::cout << "Saving raw data to " << cmdline.save_raw_data_file << std::endl;
      HDF5writer wr(cmdline.save_raw_data_file);  raw_data.write(wr, "raw_data");
    }
  }

  int nsample = (args.traj_lessthan - args.traj_start)/args.traj_inc/args.bin_size;      
  
  //Load/generate resampled data
  ResampledData<jackknifeCorrelationFunction> data_j;
  if(cmdline.load_combined_data){
    loadCheckpoint(data_j, cmdline.load_combined_data_file);
    for(int i=0;i<ops.size();i++)
      for(int j=i;j<ops.size();j++)
	assert(data_j.haveData(ops[i],ops[j]));
    
  }else{
    data_j.generatedResampledData(raw_data, args.bin_size, args.Lt, args.tsep_pipi, args.do_vacuum_subtraction, args.timeslice_avg_vac_sub);
  }  

  if(cmdline.save_combined_data) saveCheckpoint(data_j, cmdline.save_combined_data_file);

  //Put resampled data into matrices
  correlationFunction<double, NumericSquareMatrix<jackknifeDistributionD> > C(args.Lt);
  for(int t=0;t<args.Lt;t++){
    C.coord(t) = t;
    C.value(t).resize(nop);

    for(int i=0;i<nop;i++)
      for(int j=i;j<nop;j++)
	C.value(t)(i,j) = C.value(t)(j,i) = data_j.correlator(ops[i],ops[j]).value(t);
  }

  if(cmdline.subtract_from_data)
    correlatorSubtract(C, cmdline.subtract_from_data_file);

  analze_GEVP(C, args.t_max, args.fit_tmin, args.fit_tmax, args.Ascale, cmdline.verbose_solver);

  std::cout << "Done\n";
  return 0;
}


