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


void errorWeightedAverage(std::vector<correlationFunction<double, jackknifeDistributionD> > &out_asymm,
			  std::vector<correlationFunction<double, jackknifeDistributionD> > &out_symm,
			  const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &in_asymm,
			  const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &in_symm,
			  const int tmin_k_op){
  assert(in_asymm[0].size() > 0);
  const int nsample = in_asymm[0].value(0).size();
  
  for(int q=0;q<10;q++){
    std::map<int, std::vector<int> > dmap;
    for(int d=0;d<in_asymm[q].size();d++){
      const int tsep_k_op = int(in_asymm[q].coord(d).t);
      const int tsep_k_pi = in_asymm[q].coord(d).tsep_k_pi;
      const int tsep_op_pi = tsep_k_pi - tsep_k_op;
      if(tsep_op_pi >= 0 && tsep_k_op >= tmin_k_op) dmap[tsep_op_pi].push_back(d);
    }
    for(auto it = dmap.begin(); it != dmap.end(); it++){
      const int tsep_op_pi = it->first;
      const std::vector<int> &include = it->second;

      double wsum = 0.;
      jackknifeDistributionD wavg_asymm(nsample,0.);
      jackknifeDistributionD wavg_symm(nsample,0.);

      for(int dd=0; dd<include.size();dd++){
	const int d = include[dd];
	const double stderr = in_asymm[q].value(d).standardError(); //use weights from asymm for both for consistent combination
	const double w = 1/stderr/stderr;
	wsum = wsum + w;
	wavg_asymm = wavg_asymm + w*in_asymm[q].value(d);
	wavg_symm = wavg_symm + w*in_symm[q].value(d);
      }
      wavg_asymm = wavg_asymm / wsum;
      wavg_symm = wavg_symm / wsum;

      out_asymm[q].push_back(tsep_op_pi, wavg_asymm);
      out_symm[q].push_back(tsep_op_pi, wavg_symm);
    }
  }
  
}

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

  const int sz = A0_all_asymm_j[0].size();
  for(int q=0;q<10;q++){
    assert(A0_all_asymm_j[q].size() == A0_all_symm_j[q].size());
    assert(A0_all_asymm_j[q].size() == sz);
  }
  
  for(int i=0;i<args.tsep_k_pi.size();i++){
    const int tsep_k_pi = args.tsep_k_pi[i];
    for(int q=0;q<10;q++){
      std::cout << "tsep_k_pi=" << tsep_k_pi << " Q=" << q+1 << std::endl;
      correlationFunction<double, jackknifeDistributionD> A0_q_asymm_tsep_j, A0_q_symm_tsep_j;
      for(int a=0;a<sz;a++){
	if(A0_all_asymm_j[q].coord(a).tsep_k_pi == tsep_k_pi){
	  typedef correlationFunction<double, jackknifeDistributionD>::ElementType Elem;
	  assert(A0_all_asymm_j[q].coord(a) == A0_all_symm_j[q].coord(a));
	  double tsep_op_pi = tsep_k_pi - A0_all_asymm_j[q].coord(a).t;
	  if(tsep_op_pi < 0) continue;

	  A0_q_asymm_tsep_j.push_back(Elem(tsep_op_pi, A0_all_asymm_j[q].value(a)));
	  A0_q_symm_tsep_j.push_back(Elem(tsep_op_pi, A0_all_symm_j[q].value(a)));
	}
      }
      int fsz = A0_q_symm_tsep_j.size();
      std::cout << "Found " << fsz << " data with tsep_k_pi="<<tsep_k_pi<<" Q="<<q+1<<std::endl;
      if(fsz == 0) error_exit(std::cout << "Error: Couldn't find any data for tsep_k_pi="<<tsep_k_pi<<" Q="<<q+1<<std::endl);

      std::ostringstream plot_stub; plot_stub << "reldiff_Q" << q+1 << "_tsep_k_pi_" << tsep_k_pi;
      compareRelativeDifferences(A0_q_asymm_tsep_j, A0_q_symm_tsep_j, plot_stub.str());
    }
  }

  //Divide out the kaon mass exponential term so the data is dependent only on tsep_op_pi (in the kaon plateaux region)
  jackknifeDistributionD mK_asymm, mK_symm;
  {
    std::vector<jackknifeDistributionD> tmp;
    readParamsStandard(tmp,args.mK_file_asymm);
    mK_asymm = std::move(tmp[args.mK_param_asymm]);
  }
  {
    std::vector<jackknifeDistributionD> tmp;
    readParamsStandard(tmp,args.mK_file_symm);
    mK_symm = std::move(tmp[args.mK_param_symm]);
  }
  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > A0_noKexp_all_asymm_j = A0_all_asymm_j;
  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > A0_noKexp_all_symm_j = A0_all_symm_j;

  for(int q=0;q<10;q++){
    for(int a=0;a<sz;a++){
      double t = A0_all_asymm_j[q].coord(a).t;
      A0_noKexp_all_asymm_j[q].value(a) = A0_noKexp_all_asymm_j[q].value(a)/exp(-mK_asymm * t);
      A0_noKexp_all_symm_j[q].value(a) = A0_noKexp_all_symm_j[q].value(a)/exp(-mK_symm * t);
    }
  }

  //Do the error weighted averages
  std::vector<correlationFunction<double, jackknifeDistributionD> > errW_asymm(10), errW_symm(10);
  errorWeightedAverage(errW_asymm, errW_symm, A0_noKexp_all_asymm_j, A0_noKexp_all_symm_j, args.weighted_avg_tmin_k_op);

  for(int q=0;q<10;q++){
    std::cout << "Weighted average of Q="<<q+1 <<" data with kaon time dependence divided out\n";
    std::ostringstream plot_stub; plot_stub << "reldiff_Q" << q+1 << "_errW";
    compareRelativeDifferences(errW_asymm[q], errW_symm[q], plot_stub.str());
  }

  std::cout << "Done" << std::endl;
  
  return 0;
}

