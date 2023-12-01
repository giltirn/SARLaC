#pragma once

#include "../figure_data_container.h"
#include "../../correlator_utils/threemomentum.h"
#include "../../correlator_utils/mom_project.h"

SARLAC_START_NAMESPACE

///////////////////
//FILENAME POLICIES
//////////////////

struct PiPiToSigmaGenericReadPolicy{
  int traj_start, traj_inc, traj_lessthan;
  std::string dir;
  subStringReplace fmt;
  
  PiPiToSigmaGenericReadPolicy(const std::string &file_fmt, const std::string &dir, const int traj_start, const int traj_inc, const int traj_lessthan):
    traj_start(traj_start), traj_inc(traj_inc), traj_lessthan(traj_lessthan), dir(dir){
    fmt.chunkString(file_fmt, { subStringSpecify("<CONF>"), subStringSpecify("<MOM_QUARK_SIGMA>"), subStringSpecify("<MOM_PI>") });
  }

  int nsample() const{ return (traj_lessthan - traj_start)/traj_inc; }
  
  std::string filename(const int sample, const threeMomentum &mom_quark_sigma, const threeMomentum &mom_pi) const{
    std::ostringstream os;
    os << dir << '/';
    fmt.replace(os, { anyToStr(traj_start + sample * traj_inc), momStr(mom_quark_sigma), momStr(mom_pi) } );
    return os.str();
  }
};
  
struct PiPiToSigmaBasicReadPolicy: public PiPiToSigmaGenericReadPolicy{
  PiPiToSigmaBasicReadPolicy(const std::string &dir, const int traj_start, const int traj_inc, const int traj_lessthan):
    PiPiToSigmaGenericReadPolicy("traj_<CONF>_pipitosigma_sigmawdagmom<MOM_QUARK_SIGMA>_pionmom<MOM_PI>_v2", dir, traj_start, traj_inc, traj_lessthan)
  {}
};
  
////////////////////////
//READ PIPI->SIGMA DATA
//This includes the disconnected component. 
//In practise we only measure on a subset of pipi src timeslices. To improve the statistics we can reconstruct the disconnected contribution from the missing timeslices
//using the pipi and sigma bibble data via code in this directory
//////////////////////

template<typename ReadPolicy>
void readPiPiToSigma(figureData &raw_data, const int Lt, const PiPiProjectorBase &proj_pipi, const ReadPolicy &rd){
  std::cout << "Reading pipi->sigma data\n"; boost::timer::auto_cpu_timer t("Read pipi->sigma in %w s\n");
  int nsample = rd.nsample();

  raw_data.setup(Lt,rawDataDistributionD(nsample));
  raw_data.zero();

  std::vector<threeMomentum> quark_mom = { {1,1,1}, {-1,-1,-1},
					   {-3,1,1}, {3,-1,-1},
					   {1,-3,1}, {-1,3,-1},
					   {1,1,-3}, {-1,-1,3} };

  
  figureData tmp_raw_data(Lt,rawDataDistributionD(nsample));

  for(int ppiidx = 0 ; ppiidx < proj_pipi.nMomenta() ; ppiidx++){
    threeMomentum ppi = proj_pipi.momentum(ppiidx) * 2; //Tianle's conventions for the pion energy are in units of pi/2L     
    double pipi_coeff = std::real(proj_pipi.coefficient(ppiidx));    

    for(int psigqidx = 0 ; psigqidx < 8 ; psigqidx++){
#pragma omp parallel for
      for(int sample=0; sample < nsample; sample++){
	std::string filename = rd.filename(sample, quark_mom[psigqidx], ppi);
	std::cout << "Parsing " << filename << std::endl;
	tmp_raw_data.parseCDR(filename, sample);
      }

      raw_data = raw_data + 1./8 * pipi_coeff*tmp_raw_data;
    }
  }
}
//Call the above with the default file format
void readPiPiToSigma(figureData &raw_data, const std::string &data_dir, const int Lt, const PiPiProjectorBase &proj_pipi,
		     const int traj_start, const int traj_inc, const int traj_lessthan){ 
  PiPiToSigmaBasicReadPolicy rd(data_dir, traj_start, traj_inc, traj_lessthan);
  readPiPiToSigma(raw_data, Lt, proj_pipi, rd);
}
//Call the above with a user-specified file format
void readPiPiToSigma(figureData &raw_data, const std::string &file_fmt, const std::string &data_dir, const int Lt, const PiPiProjectorBase &proj_pipi,
		     const int traj_start, const int traj_inc, const int traj_lessthan){ 
  PiPiToSigmaGenericReadPolicy rd(file_fmt, data_dir, traj_start, traj_inc, traj_lessthan);
  readPiPiToSigma(raw_data, Lt, proj_pipi, rd);
}



SARLAC_END_NAMESPACE
