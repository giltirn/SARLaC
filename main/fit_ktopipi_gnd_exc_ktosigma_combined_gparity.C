#include<ktopipi_common/ktopipi_common.h>

#include<fit_ktopipi_ktosigma_combined_gparity/simfit_generic.h>
#include<fit_ktopipi_gnd_exc_ktosigma_combined_gparity/args.h>
#include<fit_ktopipi_gnd_exc_ktosigma_combined_gparity/cmdline.h>
#include<fit_ktopipi_gnd_exc_ktosigma_combined_gparity/resampled_data.h>
#include<fit_ktopipi_gnd_exc_ktosigma_combined_gparity/simfit_plot.h>
#include<fit_ktopipi_gnd_exc_ktosigma_combined_gparity/simfit.h>

using namespace CPSfit;






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

  int nsample = (args.traj_lessthan - args.traj_start)/args.traj_inc/args.bin_size;

  simultaneousFitBase<jackknifeDistribution>* fitter = getFitter<jackknifeDistribution>(args.fitfunc, args.nstate);

  fitter->load2ptFitParams(args.operators, args.input_params, nsample); 

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
  ResampledData<jackknifeDistributionD> data_j;
  ResampledData<doubleJackknifeA0StorageType> data_dj;
  ResampledData<blockDoubleJackknifeA0StorageType> data_bdj;
  
  bool do_dj = args.covariance_matrix == CovarianceMatrix::Regular;
  bool do_bdj = args.covariance_matrix == CovarianceMatrix::Block;

  if(do_dj) std::cout << "Using double jackknife" << std::endl;
  if(do_bdj) std::cout << "Using block double jackknife" << std::endl;

  if(cmdline.load_resampled_data_container_checkpoint){
    std::cout << "Reading resampled data container from checkpoint file" << std::endl;
    HDF5reader rd(cmdline.load_resampled_data_container_checkpoint_file);
    read(rd, data_j, "data_j");
    if(do_dj) read(rd, data_dj, "data_dj");
    if(do_bdj) read(rd, data_bdj, "data_bdj");
  }else{
    data_j.resample(raw, args, cmdline, "single jackknife");
    if(do_dj) data_dj.resample(raw, args, cmdline, "double jackknife");
    if(do_bdj) data_bdj.resample(raw, args, cmdline, "block double jackknife");
  }
  
  if(cmdline.save_resampled_data_container_checkpoint){
    std::cout << "Writing resampled data container to checkpoint file" << std::endl;
    HDF5writer wr(cmdline.save_resampled_data_container_checkpoint_file);
    write(wr, data_j, "data_j");
    if(do_dj) write(wr, data_dj, "data_dj");
    if(do_bdj) write(wr, data_bdj, "data_bdj");
  }

  if(args.basis == Basis::Basis7){
    std::cout << "Converting to 7-basis" << std::endl;
    data_j.convertBasis10to7(); 
    if(do_dj) data_dj.convertBasis10to7();
    if(do_bdj) data_bdj.convertBasis10to7();
  }

  std::cout << "Starting fits" << std::endl;
  typedef taggedValueContainer<double,std::string> Params;

  ResampledDataContainers<jackknifeDistribution> rdata(data_j, data_dj, data_bdj);

  std::vector<jackknifeDistribution<Params> > params = fitter->fit(rdata, args.operators,
								   args.Lt, args.tmin_k_op, args.tmin_op_snk, args.correlated, args.covariance_matrix);

  if(args.basis == Basis::Basis7){
    std::cout << "Converting 7 basis results to 10 basis" << std::endl;
    std::vector<std::vector<jackknifeDistributionD> > params_10 = convert7basisTo10basis(args.nstate, params);
    writeParamsStandard(params_10, "matrix_elems_10basis.hdf5");
    for(int q=0;q<10;q++){
      std::cout << "Q" << q+1 << std::endl;
      for(int m=0;m<args.nstate;m++)
	std::cout << "M" << m << " = " << params_10[q][m] << std::endl;
    }    
  }

  std::cout << "Done" << std::endl;
  
  delete fitter;

  return 0;
}
