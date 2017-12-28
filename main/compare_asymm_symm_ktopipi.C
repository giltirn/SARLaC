#include<array>
#include<vector>
#include<complex>
#include<iostream>
#include<fstream>

#include <boost/timer/timer.hpp>

#include<utils.h>
#include<distribution.h>
#include<common_defs.h>
#include<numeric_tensors.h>
#include<correlationfunction.h>
#include<fit_wrapper.h>
#include<parser.h>
#include<hdf5_serialize.h>

#include <fit_pipi_gparity/data_containers.h>
#include <fit_pipi_gparity/mom_data_containers.h>
#include <fit_pipi_gparity/read_data.h>

#include <fit_ktopipi_gparity/freeze.h>
#include <fit_ktopipi_gparity/cmdline.h>
#include <fit_ktopipi_gparity/args.h>
#include <fit_ktopipi_gparity/data_containers.h>
#include <fit_ktopipi_gparity/read_data.h>
#include <fit_ktopipi_gparity/compute_amplitude.h>
#include <fit_ktopipi_gparity/fitfunc.h>
#include <fit_ktopipi_gparity/plot.h>
#include <fit_ktopipi_gparity/fit.h>
#include <fit_ktopipi_gparity/main.h>

#include <compare_asymm_symm_ktopipi/cmdline.h>
#include <compare_asymm_symm_ktopipi/args.h>
#include <compare_simple_correlators/compare.h>

int main(const int argc, const char* argv[]){
  printMem("Beginning of execution");
  
  ComparisonArgs args;
  if(argc < 2){
    std::ofstream of("template.args");
    of << args;
    std::cout << "Wrote template argument file to template.args\n";
    return 0;
  }
  const std::string arg_file = argv[1];
  parse(args, arg_file);
  
  assert( (args.traj_lessthan - args.traj_start) % args.traj_inc == 0 );  
  const int nsample = (args.traj_lessthan - args.traj_start)/args.traj_inc;

  ComparisonCMDline cmdline(argc,argv,2);
  
  //Prepare the data
  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > A0_all_asymm_j(10), A0_all_symm_j(10);
  {
    CMDline c = cmdline.toCMDline(Asymmetric);
    Args a = args.toArgs(Asymmetric);
    std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > A0_all_dj(10);
    getData(A0_all_asymm_j, A0_all_dj,a,c);
  }
  {
    CMDline c = cmdline.toCMDline(Symmetric);
    Args a = args.toArgs(Symmetric);
    std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > A0_all_dj(10);
    getData(A0_all_symm_j, A0_all_dj,a,c);
  }

  assert(A0_all_asymm_j.size() == A0_all_symm_j.size());
  const int sz = A0_all_asymm_j.size();
  
  for(int i=0;i<args.tsep_k_pi.size();i++){
    const int tsep_k_pi = args.tsep_k_pi[i];
    for(int q=0;q<10;q++){
      correlationFunction<double, jackknifeDistributionD> A0_q_asymm_tsep_j, A0_q_symm_tsep_j;
      for(int a=0;a<sz;a++){
	if(A0_all_asymm_j[q].coord(a).tsep_k_pi == tsep_k_pi){
	  typedef correlationFunction<double, jackknifeDistributionD>::ElementType Elem;
	  assert(A0_all_asymm_j[q].coord(a) == A0_all_symm_j[q].coord(a));
	  double t = A0_all_asymm_j[q].coord(a).t;
	  
	  A0_q_asymm_tsep_j.push_back(Elem(t, A0_all_asymm_j[q].value(a)));
	  A0_q_symm_tsep_j.push_back(Elem(t, A0_all_symm_j[q].value(a)));
	}
      }
      std::ostringstream plot_stub; plot_stub << "reldiff_Q" << q+1 << "_tsep_k_pi_" << tsep_k_pi;
      compareRelativeDifferences(A0_q_asymm_tsep_j, A0_q_symm_tsep_j, plot_stub.str());
    }
  }

  std::cout << "Done" << std::endl;
  
  return 0;
}

