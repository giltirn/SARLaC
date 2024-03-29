#include <pipi_common/pipi_common.h>

using namespace SARLaC;

#include <fit_pipitosigma_gparity/fit.h>
#include <fit_pipitosigma_gparity/args.h>
#include <fit_pipitosigma_gparity/cmdline.h>


int main(const int argc, const char* argv[]){
  PiPiToSigmaArgs args;
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
  
  PiPiToSigmaCMDline cmdline(argc,argv,2);

  //Read sigma bubble as complex
  sigmaSelfContractionZ sigma_self_data_Z;
  readSigmaSelf(sigma_self_data_Z, args.data_dir, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);

  sigmaSelfContraction sigma_self_data = reIm(sigma_self_data_Z, 0); //copy real part

  //Read pipi bubble as complex
  PiPiProjectorA1Basis111 proj_pipi;
  bubbleDataZ pipi_self_data_Z;
  getProjectedSourcePiPiBubble(pipi_self_data_Z, args.data_dir, args.traj_start, args.traj_inc, args.traj_lessthan, args.Lt, args.tsep_pipi, proj_pipi);

  bubbleData pipi_self_data = reIm(pipi_self_data_Z, 0); //real part

  readReconstructPiPiToSigmaWithDisconnAllTsrcOptions opt;
  opt.force_disconn_tstep_src = cmdline.force_disconn_tstep_src;
  opt.disconn_tstep_src = cmdline.disconn_tstep_src;
  opt.compute_disconn_ReRe = !cmdline.use_disconn_complex_prod;

  rawDataCorrelationFunctionD correlator_raw = readReconstructPiPiToSigmaWithDisconnAllTsrc(args.data_dir, args.Lt, args.tstep_src, proj_pipi, args.traj_start, args.traj_inc, args.traj_lessthan,
										       pipi_self_data_Z, sigma_self_data_Z,
										       opt);

  //Resample data
  auto correlator_j = binResample<jackknifeCorrelationFunctionD>(correlator_raw, args.bin_size);
  auto correlator_dj = binResample<doubleJackknifeCorrelationFunctionD>(correlator_raw, args.bin_size);

  const int nsample = correlator_dj.value(0).size();

  if(args.do_vacuum_subtraction){
    correlator_j = correlator_j - computePiPiToSigmaVacSub<jackknifeCorrelationFunctionD>(sigma_self_data, pipi_self_data, args.bin_size);
    correlator_dj = correlator_dj - computePiPiToSigmaVacSub<doubleJackknifeCorrelationFunctionD>(sigma_self_data, pipi_self_data, args.bin_size);
  }

  correlator_j = fold(correlator_j, args.tsep_pipi);
  correlator_dj = fold(correlator_dj, args.tsep_pipi);

  std::cout << "Resampled data:\n";
  for(int t=0;t<args.Lt;t++) std::cout << t << " " << correlator_j.value(t) << std::endl;

  //Filter out the data that is to be fitted
  filterXrange<double> trange(args.t_min,args.t_max);
  
  doubleJackknifeCorrelationFunctionD correlator_dj_inrange;
  jackknifeCorrelationFunctionD correlator_j_inrange;
  for(int d=0;d<correlator_dj.size();d++)
    if(trange.accept(correlator_dj.coord(d),correlator_dj.value(d) )){
      correlator_dj_inrange.push_back(correlator_dj[d]);
      correlator_j_inrange.push_back(correlator_j[d]);
    }
  
  //Perform the fit
  PiPiToSigmaFitArgs fargs;
  args.transfer(fargs);
  cmdline.transfer(fargs);
  std::pair<jackknifeDistributionD,jackknifeDistributionD> Epipi_and_const = fit_pipi_to_sigma(correlator_j_inrange, correlator_dj_inrange, fargs);
  
  std::cout << "Done\n";
  return 0;
}

