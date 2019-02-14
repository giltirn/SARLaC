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

  doubleJackknifeCorrelationFunctionD data_dj;
  jackknifeCorrelationFunctionD data_j;

  if(cmdline.load_combined_data){
#ifdef HAVE_HDF5
    std::cout << "Reading resampled data from " << cmdline.load_combined_data_file << std::endl;
    HDF5reader reader(cmdline.load_combined_data_file);
    read(reader, data_j, "data_j");
    if(args.covariance_strategy != CovarianceStrategy::FrozenCorrelated) read(reader, data_dj, "data_dj");
#else
    error_exit("main: Loading amplitude data requires HDF5\n");
#endif
  }else{
    const int nchannel = args.data.size();
    std::vector<rawDataCorrelationFunctionD> channels_raw(nchannel);
    if(cmdline.load_raw_data){ //Load from checkpoint if desired
#ifdef HAVE_HDF5
      std::cout << "Reading raw data from " << cmdline.load_raw_data_file << std::endl;
      HDF5reader reader(cmdline.load_raw_data_file);
      read(reader, channels_raw, "channels_raw");
#else
      error_exit("main: Loading raw data requires HDF5\n");
#endif
    }else{ //Load from original files
      for(int i=0;i<nchannel;i++)
	readData(channels_raw[i], args.data[i], args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);
    }
    
    if(cmdline.save_raw_data){
#ifdef HAVE_HDF5
      std::cout << "Writing raw data to " << cmdline.save_raw_data_file << std::endl;
      HDF5writer writer(cmdline.save_raw_data_file);
      write(writer, channels_raw, "channels_raw");
#else
      error_exit("main: Writing raw data requires HDF5\n");
#endif
    }

    if(cmdline.remove_samples_in_range){
      std::cout << "Removing samples in range [" << cmdline.remove_samples_in_range_start << ", " <<  cmdline.remove_samples_in_range_lessthan << ")" << std::endl;
      for(int c=0;c<nchannel;c++)
	for(int t=0;t<channels_raw[c].size();t++)
	  channels_raw[c].value(t) = removeSamplesInRange(channels_raw[c].value(t), cmdline.remove_samples_in_range_start, cmdline.remove_samples_in_range_lessthan);
    }
    if(cmdline.scramble_raw_data){ //useful as a check to see if binning is actually doing anything more than reducing resolution on the covariance matrix
      int nsample = channels_raw[0].value(0).size();
      if(!RNG.isInitialized()) RNG.initialize(1234);
      std::vector<int> reord(nsample);
      std::list<int> rem; 
      for(int i=0;i<nsample;i++) rem.push_back(i);
      
      for(int i=0;i<nsample;i++){
	int off = (int)uniformRandom<float>(0,rem.size());
	auto it = std::next(rem.begin(), off);
	reord[i] = *it;
	rem.erase(it);
      }
      
      //Check indices are unique
      std::cout << "Reordered samples: ";
      std::set<int> con;
      for(int i=0;i<nsample;i++){
	std::cout << reord[i] << " ";
	con.insert(reord[i]);
      }
      std::cout << std::endl;
      assert(con.size() == nsample); 

      for(int c=0;c<nchannel;c++)
	for(int t=0;t<channels_raw[c].size();t++)
	  channels_raw[c].value(t) = rawDataDistributionD(nsample, [&](const int s){ return channels_raw[c].value(t).sample(reord[s]); });
    }

    data_j = resampleAndCombine<jackknifeDistributionD>(channels_raw, args, cmdline);
    if(args.covariance_strategy != CovarianceStrategy::FrozenCorrelated) data_dj = resampleAndCombine<doubleJackknifeDistributionD>(channels_raw, args, cmdline);
  }

  fit(data_j,data_dj, args, cmdline);
}

