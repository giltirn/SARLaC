#ifndef _SIGMA_READ_DATA_H_
#define _SIGMA_READ_DATA_H_

std::string sigmaFile(const std::string &data_dir, const int traj, const threeMomentum &psnk_quark, const threeMomentum &psrc_quark){
  std::ostringstream os;
  os << data_dir << "/traj_" << traj << "_sigmacorr_mompsrc"<< momStr(psrc_quark) << "psnk" << momStr(psnk_quark) << "_v2";
  return os.str();
}


void readSigmaSigma(figureData &raw_data, const std::string &data_dir, const int Lt,
		    const int traj_start, const int traj_inc, const int traj_lessthan){
  std::cout << "Reading sigma 2pt data\n"; boost::timer::auto_cpu_timer t("Read sigma 2pt in %w s\n");
  int nsample = (traj_lessthan - traj_start)/traj_inc;

  raw_data.setup(Lt,nsample);
  raw_data.zero();

  std::vector<threeMomentum> quark_mom = { {1,1,1}, {-1,-1,-1},
					   {-3,1,1}, {3,-1,-1},
					   {1,-3,1}, {-1,3,-1},
					   {1,1,-3}, {-1,-1,3} };
  figureData tmp_raw_data(Lt,nsample);

  const int nmom = quark_mom.size();
  
  for(int psnk=0;psnk<nmom;psnk++){
    for(int psrc=0;psrc<nmom;psrc++){
            
#pragma omp parallel for
      for(int sample=0; sample < nsample; sample++){
	int traj = traj_start + sample * traj_inc;
	std::string filename = sigmaFile(data_dir, traj, quark_mom[psnk], quark_mom[psrc]);
	std::cout << "Parsing " << filename << std::endl;
	tmp_raw_data.parseCDR(filename, sample);
      }

      raw_data = raw_data + tmp_raw_data;
    }
  }
  
  raw_data = raw_data/double(nmom*nmom);
}


class sigmaSelfContraction{
  typedef rawDataDistribution<double> DistributionType;
  NumericVector<DistributionType> d; //(t).sample(cfg)
  int Lt;

public:
  sigmaSelfContraction(){}
  sigmaSelfContraction(const int _Lt, const int _Nsample): Lt(_Lt), d(_Lt, DistributionType(_Nsample)){}

  void setup(const int _Lt, const int _Nsample){ 
    Lt = _Lt; 
    d.resize(_Lt, DistributionType(_Nsample)); 
  }
  
  inline int getLt() const{ return Lt; }
  inline int getNsample() const{ return d[0].size(); }

  inline DistributionType & operator()(const int t){ return d[t]; } //time coordinate is for pi1 (earlier pion for sink, later for src)
  inline const DistributionType & operator()(const int t) const { return d[t]; }

  inline DistributionType & at(const int t){ return d[t]; }
  inline const DistributionType & at(const int t) const { return d[t]; }

  void zero(){
    for(int t=0;t<Lt;t++) zeroit(d(t));
  }

  void parse(std::istream &in, const int sample){
    int t;
    for(int t_expect=0;t_expect<Lt;t_expect++){
      if(!(in >> t)) error_exit(std::cout << "sigmaSelfContraction::parse failed to read t for config " << sample << "\n");
      if(t != t_expect) error_exit(std::cout << "sigmaSelfContraction::parse t doesn't match expectations: " << t << ":" << t_expect << " for config " << sample << "\n");

      double &re = d(t).sample(sample);
      double im; //discard because it is zero
      if(!(in >> re >> im)) error_exit(std::cout << "sigmaSelfContraction::parse failed to real values for config " << sample << "\n");
    }
  }
  void parse(const std::string &filename, const int sample){
    std::ifstream is(filename.c_str());
    if(!is.good()){ std::cout << "Could not open file \"" << filename << "\"\n"; std::cout.flush(); exit(-1); }
    parse(is,sample);
    if(is.fail() || is.bad()){ std::cout << "Error reading file \"" << filename << "\"\n"; std::cout.flush(); exit(-1); }
    is.close();
  }

  void bin(const int bin_size){
#pragma omp parallel for
    for(int t=0;t<Lt;t++){
      d(t) = d(t).bin(bin_size);
    }
  }

  GENERATE_HDF5_SERIALIZE_METHOD((Lt)(d));
};
#ifdef HAVE_HDF5
inline void write(HDF5writer &writer, const sigmaSelfContraction &d, const std::string &tag){ d.write(writer,tag); }
inline void read(HDF5reader &reader, sigmaSelfContraction &d, const std::string &tag){ d.read(reader,tag); }
#endif







std::string sigmaSelfFile(const std::string &data_dir, const int traj, const threeMomentum &p_quark){
  std::ostringstream os;
  os << data_dir << "/traj_" << traj << "_sigmaself_mom" << momStr(p_quark) << "_v2";
  return os.str();
}

void readSigmaSelf(sigmaSelfContraction &raw_data, const std::string &data_dir, const int Lt,
		const int traj_start, const int traj_inc, const int traj_lessthan){    
  const int nsample = (traj_lessthan - traj_start)/traj_inc;

  std::vector<threeMomentum> quark_mom = { {1,1,1}, {-1,-1,-1},
					   {-3,1,1}, {3,-1,-1},
					   {1,-3,1}, {-1,3,-1},
					   {1,1,-3}, {-1,-1,3} };
  const int nmom = quark_mom.size();

  raw_data.setup(Lt,nsample);
  raw_data.zero();

  sigmaSelfContraction tmp_data(Lt,nsample);

  for(int p=0;p<nmom;p++){    
#pragma omp parallel for
    for(int sample=0; sample < nsample; sample++){
      int traj = traj_start + sample * traj_inc;
      std::string filename = sigmaSelfFile(data_dir, traj, quark_mom[p]);
      std::cout << "Parsing " << filename << std::endl;
      tmp_data.parse(filename, sample);
    }
    
    for(int t=0;t<Lt;t++) raw_data(t) = raw_data(t) + tmp_data(t) / double(nmom);
  }
}


#endif
