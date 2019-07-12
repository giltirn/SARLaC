#include<common.h>
#include<fit.h>
#include<parser.h>
#include<plot.h>
#include<random.h>

using namespace CPSfit;

#include<fit_simple/cmdline.h>
#include<fit_simple/args.h>
#include<fit_simple/read_data.h>
#include<fit_simple/fit.h>
#include<fit_simple/main.h>
#include<fit_simple/bootstrap_pvalue.h>

//Basic fitting
int main(const int argc, const char** argv){
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
    transformRaw(channels_raw, args, cmdline); 
  }

  //Resample the data
  bool do_dj, do_bdj;
  getDJtypes(do_dj, do_bdj, args.covariance_strategy);

  if(do_dj) std::cout << "Using double jackknife" << std::endl;
  if(do_bdj) std::cout << "Using block double jackknife" << std::endl;

  blockDoubleJackknifeCorrelationFunctionD data_bdj;
  doubleJackknifeCorrelationFunctionD data_dj;
  jackknifeCorrelationFunctionD data_j;

  if(cmdline.load_combined_data){
    std::cout << "Reading resampled data from " << cmdline.load_combined_data_file << std::endl;
    HDF5reader reader(cmdline.load_combined_data_file);
    read(reader, data_j, "data_j");
    if(do_dj) read(reader, data_dj, "data_dj");
    if(do_bdj) read(reader, data_bdj, "data_bdj");
  }else{
    data_j = resampleAndCombine<jackknifeDistributionD>(channels_raw, args.Lt, args.bin_size, args.combination, args.outer_time_dep);
    if(do_dj) data_dj = resampleAndCombine<doubleJackknifeDistributionD>(channels_raw, args.Lt, args.bin_size, args.combination, args.outer_time_dep);
    if(do_bdj) data_bdj = resampleAndCombine<blockDoubleJackknifeDistributionD>(channels_raw, args.Lt, args.bin_size, args.combination, args.outer_time_dep);
  }

  if(cmdline.save_combined_data){
    std::cout << "Writing resampled data to " << cmdline.save_combined_data_file << std::endl;
    HDF5writer writer(cmdline.save_combined_data_file);
    write(writer, data_j, "data_j");
    if(do_dj) write(writer, data_dj, "data_dj");
    if(do_bdj) write(writer, data_bdj, "data_bdj");
  }

  //Perform the fit
  jackknifeDistribution<parameterVectorD> params;
  jackknifeDistributionD chisq;
  int dof;

  fit(params, chisq, dof, data_j,data_dj, data_bdj, args, cmdline);

  //Compute the bootstrap p-value
  {
    std::cout << "Computing bootstrap p-value" << std::endl;
    std::cout << "Performing fit to central value of data" << std::endl;

    assert(!cmdline.load_combined_data);
    
    bool do_j_b, do_j_ub;
    getJtypes(do_j_b, do_j_ub, args.covariance_strategy);
    
    jackknifeCorrelationFunctionD data_j_unbinned;
    if(do_j_ub) data_j_unbinned = resampleAndCombine<jackknifeDistributionD>(channels_raw, args.Lt, 1, args.combination, args.outer_time_dep);

    parameterVectorD cen_params = params.mean();
    double cen_chisq;
    int cen_dof;
    
    fitCentral(cen_params, cen_chisq, cen_dof, data_j, data_j_unbinned, args, cmdline);

    double boot_p = computeBootstrapPvalue(cen_params, cen_chisq, cen_dof, channels_raw, data_j, args, cmdline);

    {
      std::ofstream of("p_boot.dat");
      of << boot_p << std::endl;
    }
    {
      std::ofstream of("chisq_cen.dat");
      of << cen_chisq << std::endl;
    }

    std::cout << "Bootstrap p-value " << boot_p << std::endl;
  }

}

