#include<common.h>
#include<fit.h>
#include<parser.h>
#include<plot.h>
#include<random.h>

using namespace SARLaC;

#include<fit_simple/cmdline.h>
#include<fit_simple/args.h>
#include<fit_simple/read_data.h>
#include<fit_simple/fit.h>
#include<fit_simple/main.h>

#include<fit_simple_sampleAMA/args.h>

//Basic fitting with 'sample AMA' technique
int main(const int argc, const char** argv){
  CMDline cmdline(argc,argv,2);
  
  ArgsSampleAMA args;
  if(argc < 2){
    std::ofstream of("template.args");
    (std::cout << "No parameter file provided: writing template to 'template.args' and exiting\n").flush();
    of << args;
    return 1;
  }    
  
  parse(args, argv[1]);

  jackknifeCorrelationFunctionD corrected_j(args.Lt);
  doubleJackknifeCorrelationFunctionD corrected_dj(args.Lt);

  if(cmdline.load_combined_data){
#ifdef HAVE_HDF5
    std::cout << "Reading resampled data from " << cmdline.load_combined_data_file << std::endl;
    HDF5reader reader(cmdline.load_combined_data_file);
    read(reader, corrected_j, "data_j");
    read(reader, corrected_dj, "data_dj");
#else
    error_exit("main: Loading amplitude data requires HDF5\n");
#endif
  }else{
    //Label data on the 'sloppy' configurations with _S and those on the 'correction' configurations as _C  
    const int nChannel = args.data.size();
    std::vector<rawDataCorrelationFunctionD> sloppy_channels_S(nChannel);
    std::vector<rawDataCorrelationFunctionD> sloppy_channels_C(nChannel);
    std::vector<rawDataCorrelationFunctionD> exact_channels_C(nChannel);

    for(int i=0;i<nChannel;i++){
      readData(sloppy_channels_S[i], args.data[i].toDataInfo(Sloppy,'S'), args.Lt, args.sloppy_traj_start, args.traj_inc, args.sloppy_traj_lessthan);

      readData(sloppy_channels_C[i], args.data[i].toDataInfo(Sloppy,'C'), args.Lt, args.correction_traj_start, args.traj_inc, args.correction_traj_lessthan);
      readData(exact_channels_C[i], args.data[i].toDataInfo(Exact,'C'), args.Lt, args.correction_traj_start, args.traj_inc, args.correction_traj_lessthan);
    }

    //Apply time-deps/combinations/binning/resampling to each set of data separately
    rawDataCorrelationFunctionD sloppy_data_S, sloppy_data_C, exact_data_C;  
    std::vector<rawDataCorrelationFunctionD>* cc[3] = {&sloppy_channels_S, &sloppy_channels_C, &exact_channels_C};
    rawDataCorrelationFunctionD* dd[3] = {&sloppy_data_S, &sloppy_data_C, &exact_data_C};
  
    for(int i=0;i<3;i++){
      applyCombination(*dd[i],*cc[i],args.combination);
      applyTimeDep(*dd[i], args.outer_time_dep, args.Lt);
      bin(*dd[i], args.bin_size);
    }
  
    //Apply a superjackknife-esque procedure to correct the sloppy data
    const int nS=sloppy_data_S.value(0).size();
    const int nC=sloppy_data_C.value(0).size();
  
    for(int t=0;t<args.Lt;t++){
      corrected_j.coord(t) = corrected_dj.coord(t) = sloppy_data_S.coord(t);

      jackknifeDistributionD sloppy_S_sj_j = superJackknifeResampleS(sloppy_data_S.value(t),nS,nC);
      doubleJackknifeDistributionD sloppy_S_sj_dj = superDoubleJackknifeResampleS(sloppy_data_S.value(t),nS,nC);    

      jackknifeDistributionD sloppy_C_sj_j = superJackknifeResampleC(sloppy_data_C.value(t),nS,nC);
      doubleJackknifeDistributionD sloppy_C_sj_dj = superDoubleJackknifeResampleC(sloppy_data_C.value(t),nS,nC);    

      jackknifeDistributionD exact_C_sj_j = superJackknifeResampleC(exact_data_C.value(t),nS,nC);
      doubleJackknifeDistributionD exact_C_sj_dj = superDoubleJackknifeResampleC(exact_data_C.value(t),nS,nC);    

      corrected_j.value(t) = sloppy_S_sj_j + exact_C_sj_j - sloppy_C_sj_j;
      corrected_dj.value(t) = sloppy_S_sj_dj + exact_C_sj_dj - sloppy_C_sj_dj;

      jackknifeDistributionD sigma_sloppy_j = doubleJackknifeDistributionD::covariance(sloppy_S_sj_dj, sloppy_S_sj_dj);
      jackknifeDistributionD sigma_corrected_j = doubleJackknifeDistributionD::covariance(corrected_dj.value(t), corrected_dj.value(t));
      std::cout << t << " " << "sloppy: value=" << sloppy_S_sj_j << " sigma=" << sigma_sloppy_j << " corrected: value=" << corrected_j.value(t) << " sigma=" << sigma_corrected_j << std::endl;
    }
  }  

  blockDoubleJackknifeCorrelationFunctionD data_bj;
  assert(args.covariance_strategy != CovarianceStrategy::CorrelatedBlockHybrid &&
	 args.covariance_strategy != CovarianceStrategy::CorrelatedBlock);

  jackknifeDistribution<parameterVectorD> params;
  jackknifeDistributionD chisq;
  int dof;
  fit(params, chisq, dof, corrected_j, corrected_dj, data_bj, args.toArgs(), cmdline);
}

