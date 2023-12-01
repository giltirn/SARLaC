#include<common.h>
#include<fit.h>
#include<parser.h>
#include<plot.h>
#include<random.h>

using namespace SARLaC;

#include<fit_simple/read_data.h>
#include<fit_simple/main.h>

#include<fit_simple_bootstrap/cmdline.h>
#include<fit_simple_bootstrap/args.h>
#include<fit_simple_bootstrap/fit.h>
#include<fit_simple_bootstrap/bootstrap_pvalue.h>

//Basic fitting
int main(const int argc, const char** argv){
  RNG.initialize(1234);
  CMDline cmdline(argc,argv,2);

  Args args;
  if(argc < 2){
    std::ofstream of("template.args");
    (std::cout << "No parameter file provided: writing template to 'template.args' and exiting\n").flush();
    of << args;
    return 1;
  }    
  
  parse(args, argv[1]);

  //Get raw data
  const int nchannel = args.data.size();
  std::vector<rawDataCorrelationFunctionD> channels_raw(nchannel);

  if(!cmdline.load_combined_data){
    //Load from checkpoint if desired
    if(cmdline.load_raw_data){ 
      std::cout << "Reading raw data from " << cmdline.load_raw_data_file << std::endl;
      HDF5reader reader(cmdline.load_raw_data_file);
      read(reader, channels_raw, "channels_raw");
    }
    //Load from original files
    else{ 
      for(int i=0;i<nchannel;i++)
	readData(channels_raw[i], args.data[i], args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);
    }

    //Save checkpoint if desired
    if(cmdline.save_raw_data){
      std::cout << "Writing raw data to " << cmdline.save_raw_data_file << std::endl;
      HDF5writer writer(cmdline.save_raw_data_file);
      write(writer, channels_raw, "channels_raw");
    }

    //Optional, in-place transformations on raw data
    transformOptions opt;
#define CP(A) opt.A = cmdline.A
    CP(remove_samples_in_range); CP(remove_samples_in_range_start); CP(remove_samples_in_range_lessthan); CP(scramble_raw_data);
#undef CP

    transformRaw(channels_raw, opt); 
  }

  int nsample = channels_raw[0].value(0).size();

  //Generate the resample table
  resampleTableOptions ropt;
  ropt.read_from_file = cmdline.load_boot_resample_table;
  ropt.read_file = cmdline.load_boot_resample_table_file;
  ropt.write_to_file = cmdline.save_boot_resample_table;
  ropt.write_file = cmdline.save_boot_resample_table_file;
  
  std::vector<std::vector<int> > rtable = generateResampleTable(nsample, args.nboot, args.resample_table_type, args.block_size, RNG, ropt);

  //Resample the data
  bootJackknifeCorrelationFunctionD data_bj;
  bootstrapCorrelationFunctionD data_b;

  bootstrapBlockResampler resampler(rtable);

  if(cmdline.load_combined_data){
    std::cout << "Reading resampled data from " << cmdline.load_combined_data_file << std::endl;
    HDF5reader reader(cmdline.load_combined_data_file);
    read(reader, data_b, "data_b");
    read(reader, data_bj, "data_bj");
  }else{
    std::cout << "Resampling bootstrap data" << std::endl;
    data_b = resampleAndCombine<bootstrapDistributionD>(channels_raw, args.Lt, args.combination, args.outer_time_dep, resampler);
    std::cout << "Resampling boot-jackknife data" << std::endl;
    data_bj = resampleAndCombine<bootJackknifeDistributionD>(channels_raw, args.Lt, args.combination, args.outer_time_dep, resampler);
  }

  if(cmdline.save_combined_data){
    std::cout << "Writing resampled data to " << cmdline.save_combined_data_file << std::endl;
    HDF5writer writer(cmdline.save_combined_data_file);
    write(writer, data_b, "data_b");
    write(writer, data_bj, "data_bj");
  }

  //Perform the fit
  bootstrapDistribution<parameterVectorD> params;
  bootstrapDistributionD chisq;
  int dof;

  fit(params, chisq, dof, data_b,data_bj, args, cmdline, true);

  //Compute p-value
  bootstrapPvalue(data_b, data_bj, params, chisq, args, cmdline);
}

