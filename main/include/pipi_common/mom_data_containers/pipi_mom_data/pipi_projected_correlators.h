#pragma once

#include<config.h>
#include<utils/macros.h>

#include "pipi_all_mom_data.h"

CPSFIT_START_NAMESPACE

//Perform the rotational state projection
template<typename DataAllMomentumType>
typename DataAllMomentumType::ContainerType project(const char fig, const DataAllMomentumType &raw_data, 
						    const PiPiProjectorBase &proj_src, const PiPiProjectorBase &proj_snk){
  std::cout << "Computing projection of figure " << fig << " with DataAllMomentumType = " << printType<DataAllMomentumType>() << "\n"; 
  boost::timer::auto_cpu_timer t(std::string("Report: Computed projection of figure ") + fig + " with DataAllMomentumType = " + printType<DataAllMomentumType>() + " in %w s\n");

  typename DataAllMomentumType::ContainerType out, tmp;

  for(int psnki=0; psnki<proj_snk.nMomenta();psnki++){
    threeMomentum psnk = proj_snk.momentum(psnki);

    for(int psrci=0; psrci<proj_src.nMomenta();psrci++){
      threeMomentum psrc = proj_src.momentum(psrci);

      tmp = std::real(proj_snk.coefficient(psnki)*proj_src.coefficient(psrci)) * raw_data(fig, momComb(psnk,psrc));
      if(psnki == 0 && psrci == 0) out = std::move(tmp);
      else out = out + tmp;
    }
  }

  return out;
}

//Given the parsed, raw data, compute the raw , unbinned, unresampled pipi correlation function from the underlying contraction data. This includes projecting the pipi states onto
//a user-selected linear combination (for example projecting onto the A1 cubic representation)
void combineRawPiPiContractions(rawDataCorrelationFunctionD &pipi_raw, const figureDataAllMomenta &raw_data,
			const PiPiProjectorBase &proj_src, const PiPiProjectorBase &proj_snk, const int isospin, 
			const int bin_size, const std::string &extra_descr = "", bool output_raw_data = true){
  if(isospin != 0 && isospin != 2) error_exit(std::cout << "combineRawPiPiContractions only supports isospin 0,2\n");
  const char figs[4] = {'C','D','R','V'};
  const std::vector<double> coeffs = isospin == 0 ? std::vector<double>({1., 2., -6., 3.}) : std::vector<double>({-2., 2., 0., 0.});
  std::string ee = extra_descr != "" ? "_" + extra_descr : "";  

  rawDataCorrelationFunctionD fig_corr[4];  

  for(int f=0;f<4;f++){
    figureData proj_data = project(figs[f], raw_data, proj_src, proj_snk);
    fig_corr[f] = sourceAverage(proj_data);

    if(output_raw_data){
      //These data are saved after binning
      rawDataCorrelationFunctionD realavg_b(fig_corr[f]);
      bin(realavg_b, bin_size);
      outputRawCorrelator(stringize("raw_data_%cpart%s.dat",figs[f],ee.c_str()), realavg_b, coeffs[f]);
    }
  }    
  
  pipi_raw = coeffs[0]*fig_corr[0] + coeffs[1]*fig_corr[1] + coeffs[2]*fig_corr[2] + coeffs[3]*fig_corr[3];
  std::cout << "Raw data " << extra_descr << ":\n" << pipi_raw << std::endl;
}

//Read pipi 2pt data into correlation function
void readPiPi2pt(rawDataCorrelationFunctionD &pipi_raw, bubbleDataAllMomentaZ &raw_bubble_data,
		 const std::string &data_dir, 
		 const std::string &figure_file_fmt, const std::string &bubble_file_fmt, 
		 const int tsep_pipi, const int tstep_pipi, const int Lt,
		 const int traj_start, const int traj_inc, const int traj_lessthan, 
		 const threeMomentum &ptot,
		 const PiPiProjector proj_src_t = PiPiProjector::A1momSet111, const PiPiProjector proj_snk_t = PiPiProjector::A1momSet111, const int isospin = 0,
		 bool filemap_allow_ptot_parity = false){
  std::unique_ptr<PiPiProjectorBase> proj_src( getProjector(proj_src_t, ptot) );
  std::unique_ptr<PiPiProjectorBase> proj_snk( getProjector(proj_snk_t, -ptot) );
  
  PiPiSymmetrySubsetFigureFileMapping ffn(data_dir, figure_file_fmt, traj_start, tsep_pipi, getSrcSnkMomentumSet(*proj_src, *proj_snk), 
					  ptot, MomentumUnit::PiOverL, filemap_allow_ptot_parity);

  bubbleFilenamePolicyGeneric bfn_src(bubble_file_fmt, ptot, Source);
  bubbleFilenamePolicyGeneric bfn_snk(bubble_file_fmt, -ptot, Sink);
  
  figureDataAllMomenta raw_data;
  readRawPiPi2ptData(raw_data, raw_bubble_data, ffn, bfn_src, bfn_snk, data_dir, traj_start, traj_inc, traj_lessthan, Lt, tstep_pipi, tsep_pipi, *proj_src, *proj_snk);

  //Combine diagrams to construct raw correlator
  combineRawPiPiContractions(pipi_raw, raw_data, *proj_src, *proj_snk, isospin, 1, "", false);
}
//Same as above but outputting just the real part of the bubble data
inline void readPiPi2pt(rawDataCorrelationFunctionD &pipi_raw, bubbleDataAllMomenta &raw_bubble_data,
			const std::string &data_dir, 
			const std::string &figure_file_fmt, const std::string &bubble_file_fmt, 
			const int tsep_pipi, const int tstep_pipi, const int Lt,
			const int traj_start, const int traj_inc, const int traj_lessthan,
			const threeMomentum &ptot,
			const PiPiProjector proj_src = PiPiProjector::A1momSet111, const PiPiProjector proj_snk = PiPiProjector::A1momSet111, const int isospin = 0,
			bool filemap_allow_ptot_parity = false){
  bubbleDataAllMomentaZ raw_bubble_data_Z;
  readPiPi2pt(pipi_raw, raw_bubble_data_Z, data_dir, figure_file_fmt,  bubble_file_fmt, tsep_pipi, tstep_pipi, Lt,
	      traj_start, traj_inc, traj_lessthan, ptot, proj_src, proj_snk, isospin, filemap_allow_ptot_parity);
  raw_bubble_data = reIm(raw_bubble_data_Z, 0);
}


//Read zero total momentum pipi 2pt data 
inline void readPiPi2pt(rawDataCorrelationFunctionD &pipi_raw, bubbleDataAllMomentaZ &raw_bubble_data,
		 const std::string &data_dir, 
		 const std::string &figure_file_fmt, const std::string &bubble_file_fmt, 
		 const int tsep_pipi, const int tstep_pipi, const int Lt,
		 const int traj_start, const int traj_inc, const int traj_lessthan, 
		 const PiPiProjector proj_src_t = PiPiProjector::A1momSet111, const PiPiProjector proj_snk_t = PiPiProjector::A1momSet111, const int isospin = 0){
  readPiPi2pt(pipi_raw, raw_bubble_data, data_dir, figure_file_fmt, bubble_file_fmt, tsep_pipi, tstep_pipi, Lt, 
	      traj_start, traj_inc, traj_lessthan, {0,0,0}, proj_src_t, proj_snk_t, isospin);
}

inline void readPiPi2pt(rawDataCorrelationFunctionD &pipi_raw, bubbleDataAllMomenta &raw_bubble_data,
			const std::string &data_dir, 
			const std::string &figure_file_fmt, const std::string &bubble_file_fmt, 
			const int tsep_pipi, const int tstep_pipi, const int Lt,
			const int traj_start, const int traj_inc, const int traj_lessthan,
			const PiPiProjector proj_src = PiPiProjector::A1momSet111, const PiPiProjector proj_snk = PiPiProjector::A1momSet111, const int isospin = 0){
  bubbleDataAllMomentaZ raw_bubble_data_Z;
  readPiPi2pt(pipi_raw, raw_bubble_data_Z, data_dir, figure_file_fmt,  bubble_file_fmt, tsep_pipi, tstep_pipi, Lt,
	      traj_start, traj_inc, traj_lessthan, proj_src, proj_snk, isospin);
  raw_bubble_data = reIm(raw_bubble_data_Z, 0);
}

CPSFIT_END_NAMESPACE
