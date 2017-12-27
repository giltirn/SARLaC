#ifndef _PIPI_READ_DATA_H_
#define _PIPI_READ_DATA_H_

std::string figureFile(const std::string &data_dir, const char fig, const int traj, const threeMomentum &psnk, const threeMomentum &psrc, const int tsep_pipi, const bool use_symmetric_quark_momenta){
  std::ostringstream os;
  os << data_dir << "/traj_" << traj << "_Figure" << fig << "_sep" << tsep_pipi << "_mom" << momStr(psrc) << "_mom" << momStr(psnk);
  if(use_symmetric_quark_momenta) os << "_symm";
  return os.str();
}

void readFigure(figureDataAllMomenta &raw_data, const char fig, const std::string &data_dir, const int tsep_pipi, const int Lt,
		const int traj_start, const int traj_inc, const int traj_lessthan, const bool use_symmetric_quark_momenta){
  std::cout << "Reading figure " << fig << "\n"; boost::timer::auto_cpu_timer t(std::string("Report: Read figure ") + fig + " in %w s\n");
  int nsample = (traj_lessthan - traj_start)/traj_inc;

  raw_data.setup(Lt,nsample);
  
  threeMomentum R[8] = { {1,1,1}, {-1,-1,-1},
			 {1,1,-1}, {-1,-1,1},
			 {1,-1,1}, {-1,1,-1},
			 {-1,1,1}, {1,-1,-1} };
  for(int psnk=0;psnk<8;psnk++)
    for(int psrc=0;psrc<8;psrc++){
      figureData &into = raw_data(fig, momComb(R[psnk], R[psrc]));      
      
#pragma omp parallel for
      for(int sample=0; sample < nsample; sample++){
	int traj = traj_start + sample * traj_inc;
	std::string filename = figureFile(data_dir, fig, traj, R[psnk], R[psrc], tsep_pipi, use_symmetric_quark_momenta);
	std::cout << "Parsing " << filename << std::endl;
	into.parseCDR(filename, sample);
      }
    }
}

std::string bubbleFile(const std::string &data_dir, const int traj, const threeMomentum &p, const int tsep_pipi, const bool use_symmetric_quark_momenta){
  std::ostringstream os;
  os << data_dir << "/traj_" << traj << "_FigureVdis_sep" << tsep_pipi << "_mom" << momStr(p);
  if(use_symmetric_quark_momenta) os << "_symm";
  return os.str();
}

void readBubble(bubbleDataAllMomenta &raw_data, const std::string &data_dir, const int tsep_pipi, const int Lt,
		const int traj_start, const int traj_inc, const int traj_lessthan, const bool use_symmetric_quark_momenta){
  std::cout << "Reading bubble\n"; boost::timer::auto_cpu_timer t("Report: Read bubble in %w s\n");
  int nsample = (traj_lessthan - traj_start)/traj_inc;

  raw_data.setup(Lt,nsample);
  
  threeMomentum R[8] = { {1,1,1}, {-1,-1,-1},
			 {1,1,-1}, {-1,-1,1},
			 {1,-1,1}, {-1,1,-1},
			 {-1,1,1}, {1,-1,-1} };
  for(int p=0;p<8;p++){
    bubbleData &into = raw_data(R[p]);
    
#pragma omp parallel for
    for(int sample=0; sample < nsample; sample++){
      int traj = traj_start + sample * traj_inc;
      std::string filename = bubbleFile(data_dir, traj, R[p], tsep_pipi, use_symmetric_quark_momenta);
      std::cout << "Parsing " << filename << std::endl;
      into.parse(filename, sample);
    }
  }
}

#endif
