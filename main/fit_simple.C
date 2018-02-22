#include<common.h>
#include<fit.h>
#include<parser.h>
#include<plot.h>

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

  const int nchannel = args.data.size();
  std::vector<doubleJackknifeCorrelationFunctionD> channels_dj(nchannel);
  std::vector<jackknifeCorrelationFunctionD> channels_j(nchannel);
  for(int i=0;i<nchannel;i++){
    rawDataCorrelationFunctionD channel_raw;
    readData(channel_raw, args.data[i], args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);
    bin(channel_raw, args.bin_size);    

    channels_dj[i] = doubleJackknifeCorrelationFunctionD(args.Lt, [&](const int t){
	return typename doubleJackknifeCorrelationFunctionD::ElementType(t,  doubleJackknifeDistributionD(channel_raw.value(t)));
      });
    channels_j[i] = jackknifeCorrelationFunctionD(args.Lt, [&](const int t){
	return typename jackknifeCorrelationFunctionD::ElementType(t,  jackknifeDistributionD(channel_raw.value(t)));
      });
  }
  doubleJackknifeCorrelationFunctionD data_dj;
  applyCombination(data_dj,channels_dj,args.combination);
  applyTimeDep(data_dj, args.outer_time_dep, args.Lt);

  jackknifeCorrelationFunctionD data_j;
  applyCombination(data_j,channels_j,args.combination);
  applyTimeDep(data_j, args.outer_time_dep, args.Lt);

  fit(data_j,data_dj, args, cmdline);
}

