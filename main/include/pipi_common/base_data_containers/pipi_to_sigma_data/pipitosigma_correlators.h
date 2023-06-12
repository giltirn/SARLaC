#include<config.h>
#include<utils/macros.h>

#include "pipitosigma_reconstruct_disconn.h"
#include "../../mom_data_containers/pipi_bubble_mom_data_container.h"
#include "../../mom_data_containers/pipi_mom_data/pipi_projected_bubble.h"

CPSFIT_START_NAMESPACE

struct readReconstructPiPiToSigmaWithDisconnAllTsrcOptions{
  //Choose a tstep for the vacuum diagram that is not 1 for testing
  bool force_disconn_tstep_src;
  int disconn_tstep_src;
  
  //Choose whether to compute the disconnected part using Re ( pipi_bubble * sigma_bubble ) or Re ( pipi_bubble ) * Re ( sigma_bubble )
  //The former is how it was computed in the parallel code, but the latter (default) has very slightly better stat error
  bool compute_disconn_ReRe;

  bool include_V_diagram;

  readReconstructPiPiToSigmaWithDisconnAllTsrcOptions(): force_disconn_tstep_src(false), compute_disconn_ReRe(true), include_V_diagram(true){}
};

//Combine above, reading the pipi->sigma data and reconstructing the disconnected component for all tsrc
//pipi_self_data_Z should be pre-projected
template<typename ReadPolicy>
rawDataCorrelationFunctionD readReconstructPiPiToSigmaWithDisconnAllTsrc(const ReadPolicy &rp,
								    const int Lt, const int tstep_src, const PiPiProjectorBase &proj_pipi,
								    const bubbleDataZ &pipi_self_data_Z, const sigmaSelfContractionZ &sigma_self_data_Z,
								    const readReconstructPiPiToSigmaWithDisconnAllTsrcOptions &opt = readReconstructPiPiToSigmaWithDisconnAllTsrcOptions()){

  bubbleData pipi_self_data = reIm(pipi_self_data_Z, 0); //real part
  sigmaSelfContraction sigma_self_data = reIm(sigma_self_data_Z, 0); //copy real part

  //Get pipi->sigma data
  figureData pipitosigma_data;
  readPiPiToSigma(pipitosigma_data, Lt, proj_pipi, rp);
  
  //Reconstruct disconnected and connected part
  figureData pipitosigma_disconn_data_ReZZ;
  reconstructPiPiToSigmaDisconnected(pipitosigma_disconn_data_ReZZ, pipi_self_data_Z, sigma_self_data_Z); // Re ( pipi_bubble * sigma_bubble )
 
  figureData pipitosigma_disconn_data_ReZReZ;
  reconstructPiPiToSigmaDisconnected(pipitosigma_disconn_data_ReZReZ, pipi_self_data, sigma_self_data); // Re ( pipi_bubble ) * Re ( sigma_bubble )

  figureData pipitosigma_conn_data;
  reconstructPiPiToSigmaConnected(pipitosigma_conn_data, pipitosigma_data, pipitosigma_disconn_data_ReZZ, tstep_src);

  //(Very slightly) better statistics if we use the Re ( pipi_bubble ) * Re ( sigma_bubble ) for the disconnected part, taking advantage of the fact that the bubbles are real under the ensemble avg
  figureData &pipitosigma_disconn_data = opt.compute_disconn_ReRe ? pipitosigma_disconn_data_ReZReZ : pipitosigma_disconn_data_ReZZ;

  //The code computes the disconnected component for all tsrc, but this option can be used to constrain the number of source timeslices to observe the effect
  if(opt.force_disconn_tstep_src){
    for(int tsrc=0;tsrc<Lt;tsrc++){
      if(tsrc % opt.disconn_tstep_src != 0)
	for(int t=0;t<Lt;t++) pipitosigma_disconn_data(tsrc, t).zero();
    }
  }

  //Source avg connected and disconnected parts and sum the contributions
  rawDataCorrelationFunctionD correlator_raw_conn = sourceAverage(pipitosigma_conn_data);
  rawDataCorrelationFunctionD correlator_raw_disconn = sourceAverage(pipitosigma_disconn_data);

  rawDataCorrelationFunctionD correlator_raw = correlator_raw_conn;
  if(opt.include_V_diagram)
    for(int t=0;t<Lt;t++) correlator_raw.value(t) = correlator_raw.value(t) + correlator_raw_disconn.value(t);

  std::cout << "Pipi->sigma raw data connected/disconnected parts:\n";
  for(int t=0;t<Lt;t++) std::cout << t << " " << correlator_raw_conn.value(t) << " " << correlator_raw_disconn.value(t) << std::endl;

  return correlator_raw;
}

rawDataCorrelationFunctionD readReconstructPiPiToSigmaWithDisconnAllTsrc(const std::string &data_dir, const int Lt, const int tstep_src, const PiPiProjectorBase &proj_pipi,
								    const int traj_start, const int traj_inc, const int traj_lessthan,
								    const bubbleDataZ &pipi_self_data_Z, const sigmaSelfContractionZ &sigma_self_data_Z,
								    const readReconstructPiPiToSigmaWithDisconnAllTsrcOptions &opt = readReconstructPiPiToSigmaWithDisconnAllTsrcOptions()){
  PiPiToSigmaBasicReadPolicy rd(data_dir, traj_start, traj_inc, traj_lessthan);
  return readReconstructPiPiToSigmaWithDisconnAllTsrc(rd, Lt, tstep_src,  proj_pipi, pipi_self_data_Z, sigma_self_data_Z, opt);
}
rawDataCorrelationFunctionD readReconstructPiPiToSigmaWithDisconnAllTsrc(const std::string &file_fmt, const std::string &data_dir, const int Lt, const int tstep_src, const PiPiProjectorBase &proj_pipi,
								    const int traj_start, const int traj_inc, const int traj_lessthan,
								    const bubbleDataZ &pipi_self_data_Z, const sigmaSelfContractionZ &sigma_self_data_Z,
								    const readReconstructPiPiToSigmaWithDisconnAllTsrcOptions &opt = readReconstructPiPiToSigmaWithDisconnAllTsrcOptions()){
  PiPiToSigmaGenericReadPolicy rd(file_fmt, data_dir, traj_start, traj_inc, traj_lessthan);
  return readReconstructPiPiToSigmaWithDisconnAllTsrc(rd, Lt, tstep_src,  proj_pipi, pipi_self_data_Z, sigma_self_data_Z, opt);
}

//This version takes in the all-momentum bubble containers and does the momentum-projection of the pipi bubble diagram
rawDataCorrelationFunctionD readReconstructPiPiToSigmaWithDisconnAllTsrc(const std::string &data_dir, const int Lt, const int tstep_src, const PiPiProjectorBase &proj_pipi,
								    const int traj_start, const int traj_inc, const int traj_lessthan,								    
								    const bubbleDataAllMomentaZ &pipi_self_data_Z, const sigmaSelfContractionZ &sigma_self_data_Z,
								    const readReconstructPiPiToSigmaWithDisconnAllTsrcOptions &opt = readReconstructPiPiToSigmaWithDisconnAllTsrcOptions()){
  bubbleDataZ pipi_self_proj_Z = projectSourcePiPiBubble(pipi_self_data_Z, proj_pipi);
  return readReconstructPiPiToSigmaWithDisconnAllTsrc(data_dir, Lt, tstep_src, proj_pipi, traj_start, traj_inc, traj_lessthan, pipi_self_proj_Z, sigma_self_data_Z, opt);
}
//Same as above but with user-specified file format
rawDataCorrelationFunctionD readReconstructPiPiToSigmaWithDisconnAllTsrc(const std::string &file_fmt, const std::string &data_dir, const int Lt, const int tstep_src, const PiPiProjectorBase &proj_pipi,
								    const int traj_start, const int traj_inc, const int traj_lessthan,								    
								    const bubbleDataAllMomentaZ &pipi_self_data_Z, const sigmaSelfContractionZ &sigma_self_data_Z,
								    const readReconstructPiPiToSigmaWithDisconnAllTsrcOptions &opt = readReconstructPiPiToSigmaWithDisconnAllTsrcOptions()){
  bubbleDataZ pipi_self_proj_Z = projectSourcePiPiBubble(pipi_self_data_Z, proj_pipi);
  return readReconstructPiPiToSigmaWithDisconnAllTsrc(file_fmt, data_dir, Lt, tstep_src, proj_pipi, traj_start, traj_inc, traj_lessthan, pipi_self_proj_Z, sigma_self_data_Z, opt);
}



CPSFIT_END_NAMESPACE
