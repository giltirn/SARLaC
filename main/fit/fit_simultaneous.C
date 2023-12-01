#include<common.h>
#include<fit.h>
#include<parser.h>
#include<plot.h>

using namespace SARLaC;

#include<fit_simple/data_info.h>
#include<fit_simple/read_data.h>

#include<fit_simultaneous/fitfunc.h>
#include<fit_simultaneous/args.h>
#include<fit_simultaneous/fit.h>

//Simultaneous fits of correlation functions

int main(const int argc, const char** argv){
  Args args;
  if(argc < 2){
    std::ofstream of("template.args");
    (std::cout << "No parameter file provided: writing template to 'template.args' and exiting\n").flush();
    of << args;
    return 1;
  }    
  
  parse(args, argv[1]);

  FitCosh fcosh(args.Lt);
  FitSinh fsinh(args.Lt);

  const int ncorr = args.correlators.size();
  
  FitFunc fitfunc;
  fitfunc.setNcorrelators(ncorr);

  jackknifeSimFitCorrelationFunction data_j;
  doubleJackknifeSimFitCorrelationFunction data_dj;

  for(int c=0;c<ncorr;c++){
    switch(args.correlators[c].fitfunc){
    case FitFuncType::FCosh:
      fitfunc.setCorrelatorFitFunc(c,&fcosh); break;
    case FitFuncType::FSinh:
      fitfunc.setCorrelatorFitFunc(c,&fsinh); break;
    default:
      assert(0);
    }

    doubleJackknifeCorrelationFunctionD corrdata_dj;
    jackknifeCorrelationFunctionD corrdata_j;

    const int nchannel = args.correlators[c].data.size();
    std::vector<doubleJackknifeCorrelationFunctionD> channels_dj(nchannel);
    std::vector<jackknifeCorrelationFunctionD> channels_j(nchannel);
    for(int i=0;i<nchannel;i++){
      rawDataCorrelationFunctionD channel_raw;
      readData(channel_raw, args.correlators[c].data[i], args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);
      bin(channel_raw, args.bin_size);    

      channels_dj[i] = doubleJackknifeCorrelationFunctionD(args.Lt, [&](const int t){
	  return typename doubleJackknifeCorrelationFunctionD::ElementType(t,  doubleJackknifeDistributionD(channel_raw.value(t)));
	});
      channels_j[i] = jackknifeCorrelationFunctionD(args.Lt, [&](const int t){
	  return typename jackknifeCorrelationFunctionD::ElementType(t,  jackknifeDistributionD(channel_raw.value(t)));
	});
    }

    applyCombination(corrdata_dj,channels_dj,args.correlators[c].combination);
    applyTimeDep(corrdata_dj, args.correlators[c].outer_time_dep, args.Lt);

    applyCombination(corrdata_j,channels_j,args.correlators[c].combination);
    applyTimeDep(corrdata_j, args.correlators[c].outer_time_dep, args.Lt);

    for(int i=0;i<corrdata_j.size();i++){
      data_j.push_back(FitSimCoord(c,corrdata_j.coord(i)), std::move(corrdata_j.value(i)) );
      data_dj.push_back(FitSimCoord(c,corrdata_dj.coord(i)), std::move(corrdata_dj.value(i)) );
    }
  }

  fit(data_j,data_dj,fitfunc,args);

  std::cout << "Done" << std::endl;
  return 0;
}
