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

  simultaneousFitBase* fitter = getFitter(args.fitfunc, args.nstate);

  fitter->load2ptFitParams(args.operators, args.input_params, nsample); 

  RawData raw;
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
    
  std::cout << "Computing resampled data" << std::endl;
  ResampledData<jackknifeDistributionD> data_j;
  ResampledData<doubleJackknifeA0StorageType> data_dj;
  
  data_j.resample(raw, args, cmdline, "single jackknife");
  data_dj.resample(raw, args, cmdline, "double jackknife");

  if(args.basis == Basis::Basis7){
    std::cout << "Converting to 7-basis" << std::endl;
    data_j.convertBasis10to7(); data_dj.convertBasis10to7();
  }

  std::cout << "Starting fits" << std::endl;
  typedef taggedValueContainer<double,std::string> Params;

  std::vector<jackknifeDistribution<Params> > params = fitter->fit(data_j, data_dj, args.operators,
								   args.Lt, args.tmin_k_op, args.tmin_op_snk, args.correlated);

  if(args.basis == Basis::Basis7){
    std::cout << "Converting 7 basis results to 10 basis" << std::endl;
    std::vector<std::vector<jackknifeDistributionD> > params_10(10, std::vector<jackknifeDistributionD>(args.nstate, jackknifeDistributionD(nsample)));    
    
    struct InContainer{
      size_t s;
      size_t idx;
      const std::vector<jackknifeDistribution<Params> > &d;
      double operator[](const int q) const{ return d[q].sample(s)(idx); }
      InContainer(size_t s, size_t idx, const std::vector<jackknifeDistribution<Params> > &d): s(s), idx(idx), d(d){}
    };
    struct OutContainer{
      size_t s;
      size_t state;
      std::vector<std::vector<jackknifeDistributionD> > &d;
      double & operator[](const int q){ return d[q][state].sample(s); }
      OutContainer(size_t s, size_t state, std::vector<std::vector<jackknifeDistributionD> > &d): s(s), state(state), d(d){}
    };

    auto const* tag_map = params[0].sample(0).getTagMap();
    assert(tag_map != NULL);

    for(int i=0;i<args.nstate;i++){
      auto it = tag_map->find(stringize("M%d",i));
      assert(it != tag_map->end());
      int Midx = it->second;
      
      for(int s=0;s<nsample;s++){
	InContainer in(s, Midx, params);
	OutContainer out(s, i, params_10);
	convert7to10(out,in);
      }
    }

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
