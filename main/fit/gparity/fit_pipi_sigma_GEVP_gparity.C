#if 1

//Fixme
int main(void){
  return 0;
}

#else

#include <pipi_common/pipi_common.h>
#include <pipi_common/analyze_chisq.h>

#include<fit/GEVP.h>

using namespace SARLaC;

#include <fit_pipi_gparity/args.h>
#include <fit_pipi_gparity/cmdline.h>

#include <fit_pipi_sigma_sim_gparity/fitfunc.h>
#include <fit_pipi_sigma_sim_gparity/fit.h>
#include <fit_pipi_sigma_sim_gparity/read_data.h>

#include <fit_pipi_sigma_GEVP_gparity/args.h>
#include <fit_pipi_sigma_GEVP_gparity/cmdline.h>

void getResampledData(correlationFunction<double, NumericSquareMatrix<jackknifeDistributionD> > &C,
		      const rawData &raw_data,
		      const int Lt, const int tsep_pipi, const int bin_size, bool do_vacuum_subtraction){

  //Get double-jack data
  auto pipi_to_sigma = binResample<jackknifeCorrelationFunction>(raw_data.pipi_to_sigma_raw, bin_size);
  auto sigma2pt = binResample<jackknifeCorrelationFunction>(raw_data.sigma2pt_raw, bin_size);
  auto pipi = binResample<jackknifeCorrelationFunction>(raw_data.pipi_raw, bin_size);
  
  //Compute vacuum subtractions
  if(do_vacuum_subtraction){
    { //Pipi->sigma
      PiPiProjectorA1Basis111 proj_pipi;
      bubbleData pipi_self_proj = projectSourcePiPiBubble(raw_data.pipi_self_data, proj_pipi);
      pipi_to_sigma = pipi_to_sigma - computePiPiToSigmaVacSub<jackknifeCorrelationFunction>(raw_data.sigma_self_data, pipi_self_proj, bin_size);
    }
    { //sigma->sigma
      sigma2pt = sigma2pt - computeSigmaVacSub<jackknifeCorrelationFunction>(raw_data.sigma_self_data, bin_size);
    }
    { //Pipi->pipi
      pipi = pipi - computePiPi2ptVacSub<jackknifeCorrelationFunction>(raw_data.pipi_self_data, bin_size, tsep_pipi);
    }
  }
  
  //Fold data
  pipi_to_sigma = fold(pipi_to_sigma, tsep_pipi);
  sigma2pt = fold(sigma2pt, 0);
  pipi = fold(pipi, 2*tsep_pipi);

  C.resize(Lt);
  for(int t=0;t<Lt;t++){
    C.coord(t) = t;
    C.value(t) = { pipi.value(t), pipi_to_sigma.value(t),
		   pipi_to_sigma.value(t), sigma2pt.value(t) };
  }
}		      
		      
int main(const int argc, const char* argv[]){
  PiPiSigmaGEVPargs args;
  if(argc < 2){
    std::ofstream of("template.args");
    of << args;
    std::cout << "Wrote template argument file to template.args\n";
    return 0;
  }
  {
    const std::string arg_file = argv[1];
    parse(args,arg_file);
    std::cout << "Read arguments: \n" << args << std::endl;
  }

  PiPiSigmaGEVPcmdLine cmdline(argc,argv,2);

  correlationFunction<double, NumericSquareMatrix<jackknifeDistributionD> > C;

  if(cmdline.load_resampled_data){
    HDF5reader rd(cmdline.load_resampled_data_file);
    read(rd,C,"C");
  }else{
    rawData raw_data;
    if(cmdline.load_checkpoint){
      raw_data.readDataFromCheckpoint(cmdline.load_checkpoint_file);
    }else{
      raw_data.readDataFromOrigFiles(args.data_dir, 
				     args.pipi2pt_figure_file_fmt, args.sigma2pt_file_fmt, args.pipitosigma_file_fmt,
				     args.pipi_bubble_file_fmt, args.sigma_bubble_file_fmt,
				     args.tsep_pipi, args.tstep_pipi2pt, args.tstep_pipitosigma,
				     args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan, !cmdline.use_pipitosigma_disconn_complex_prod);
    }
    if(cmdline.save_checkpoint){
      raw_data.saveCheckpoint(cmdline.save_checkpoint_file);
    }
    getResampledData(C, raw_data, args.Lt, args.tsep_pipi, args.bin_size, args.do_vacuum_subtraction);
  }
  if(cmdline.save_resampled_data){
    HDF5writer wr(cmdline.save_resampled_data_file);
    write(wr, C, "C");
  }

  GEVPsolver<jackknifeDistributionD> gevp;
  gevp.solve(C, args.t_max);

  std::vector<std::vector<jackknifeDistributionD> > E_all(2);

  for(int t=2;t<args.t_max;t++){
    int t0=t/2;
    auto E = gevp.effectiveEnergy(t0,t);
    if(E.size() > 0){
      std::cout << t0 << " " << t;      
      for(int n=0;n<2;n++) std::cout << " " << E[n];
      std::cout << std::endl;
    }
    E_all.push_back(E);
  }
  
  writeParamsStandard(E_all,"effective_energies.hdf5");
  
  std::cout << "Amplitudes:\n";
  for(int t0=0; t0<args.t_max; t0++){
    for(int t=t0+1; t<args.t_max; t++){
      std::vector<std::vector<jackknifeDistributionD> > Coeffs_all = gevp. effectiveAmplitude(t0,t,C);     
      std::cout << t0 << " " << t << std::endl;
      for(int op=0;op<2;op++){
	for(int state=0;state<2;state++)
	  if(Coeffs_all.size() != 0){
	    Coeffs_all[op][state] = Coeffs_all[op][state]/sqrt(args.Ascale);
	    std::cout << Coeffs_all[op][state] << " ";
	  }
	  else std::cout << "- ";
	std::cout << std::endl;
      }
    }
  }

  std::cout << "Done\n";
  return 0;
}

#endif
