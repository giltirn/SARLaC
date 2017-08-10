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

#include <fit_pipi_gparity/data_containers.h>
#include <fit_pipi_gparity/mom_data_containers.h>
#include <fit_pipi_gparity/read_data.h>

typedef rawDataDistribution<std::complex<double> > rawDataDistributionZD;

struct contractionData{
  std::vector<rawDataDistributionZD> d;

  contractionData(const int nsample): d(8,rawDataDistributionZD(nsample,std::complex<double>(0.))){}

  void parse(std::istream &in, const int sample){
    for(int i=0;i<8;i++){
      double* cd = (double*) &d[i].sample(sample);	    
      if(!(in >> cd[0] >> cd[1])) error_exit(std::cout << "type1234Data::contractions::contractionData::parse failed to parse contraction data for config " << sample << "\n");
    }
  }
    
  inline int map_off(const int i, const char A, const int j, const char B) const{
    int iA = A == 'V' ? 0 : 1;
    int iB = B == 'V' ? 0 : 1;
    return i + 2*(iA + 2*(j + 2*iB));      
  }
      
  int mapping(const int i, const char A, const int j, const char B) const{
    static int map[8];
    static bool initted = false;
    if(!initted){
      map[ map_off(0,'V',0,'A') ] = 0;
      map[ map_off(0,'A',0,'V') ] = 1;
      map[ map_off(0,'V',1,'A') ] = 2;
      map[ map_off(0,'A',1,'V') ] = 3;
      map[ map_off(1,'V',0,'A') ] = 4;
      map[ map_off(1,'A',0,'V') ] = 5;
      map[ map_off(1,'V',1,'A') ] = 6;
      map[ map_off(1,'A',1,'V') ] = 7;
      initted = true;
    }
    return map[ map_off(i,A,j,B) ];
  }
    
  inline rawDataDistributionZD & operator()(const int i, const char A, const int j, const char B){
    return d[mapping(i,A,j,B)];
  }
  inline const rawDataDistributionZD & operator()(const int i, const char A, const int j, const char B) const{
    return d[mapping(i,A,j,B)];
  }
  inline contractionData operator+(const contractionData &r) const{
    contractionData out(*this);
    for(int i=0;i<8;i++) out.d[i] = out.d[i] + r.d[i];
    return out;
  }
  inline contractionData operator*(const double d) const{
    contractionData out(*this);
    for(int i=0;i<8;i++) out.d[i] = out.d[i] * d;
    return out;
  }

  bool isZero() const{
    std::complex<double> zro(0.,0.);
    for(int i=0;i<8;i++)
      for(int s=0;s<d[i].size();s++)
	if(d[i].sample(s) != zro) return false;
    return true;
  }
    
};
    
struct contractions{
  std::vector<contractionData> con;
  std::vector<rawDataDistributionZD> extra; //F0 g5, -F1 g5

  contractions(){}
  contractions(const int nsample, const int ncontract, const int nextra = 0): con(ncontract, contractionData(nsample)), extra(nextra,rawDataDistributionZD(nsample,std::complex<double>(0.))){}

  //Parse the ncontract * 8 * 2 double entries for a given tK, t, sample
  void parse(std::istream &in, const int sample){
    for(int c=0;c<con.size();c++) con[c].parse(in,sample);

    for(int e=0;e<extra.size();e++){
      double* cd = (double*) &extra[e].sample(sample);
      if(!(in >> cd[0] >> cd[1]) ) error_exit(std::cout << "type1234Data::contractions::parse failed to parse extra idx " << e << " data for config " << sample << "\n");
    }
  }

  const contractionData & C(const int i) const{ return con[i]; }
  const rawDataDistributionZD & mix() const{ return extra[0]; } //only use the F0 g5 spin-flavor structure in practise 
  
  inline bool isZero() const{ 
    for(int i=0;i<con.size();i++) if(!con[i].isZero()) return false;

    std::complex<double> zro(0.,0.);
    for(int e=0;e<extra.size();e++)
      for(int s=0;s<extra[e].size();s++)
	if(extra[e].sample(s) != zro) return false;
    return true;    
  }

  contractions operator+(const contractions &r){
    if(r.con.size() != con.size() || r.extra.size() != extra.size()) error_exit(std::cout << "contractions::operator+ cannot sum two contractions of different sizes!\n");
    contractions out(*this);
    for(int i=0;i<con.size();i++) out.con[i] = out.con[i] + r.con[i];
    for(int i=0;i<extra.size();i++) out.extra[i] = out.extra[i] + r.extra[i];
    return out;
  }
  contractions operator*(const double d){
    contractions out(*this);
    for(int i=0;i<con.size();i++) out.con[i] = out.con[i] * d;
    for(int i=0;i<extra.size();i++) out.extra[i] = out.extra[i] * d;
    return out;
  }      
};



class type1234Data{
  typedef contractions tdata; //data associated with given t, i.e. the contractions
  typedef std::vector<tdata> tKdata;  //data associated with a given tK, i.e. Lt * tdata
  
  std::vector<tKdata> data; //[tK][t]

  int type;
  int Lt;
  int nsample;
  int ncontract;
  int nextra;
public:
  type1234Data(const int _type, const int _Lt, const int _nsample): type(_type),Lt(_Lt), nsample(_nsample), nextra(0){
    switch(type){
    case 1:
    case 2:
      ncontract = 6; break;
    case 3:
    case 4:
      ncontract = 10;
      nextra = 2; break;  
    default:
      error_exit(std::cout << "type1234Data::type1234Data invalid type " << type << std::endl);
    };

    data.resize(Lt, tKdata(Lt, tdata(nsample,ncontract,nextra)));
  }
  inline int getType() const{ return type; }
  inline int getLt() const{ return Lt; }
  inline int getNsample() const{ return nsample; }
  inline int getNcontract() const{ return ncontract; }
  inline int getNextra() const{ return nextra; }

  const contractions & operator()(const int tK, const int t) const{ return data[tK][t]; }
  
  void parse(std::istream &in, const int sample){
    int tK, t;
    for(int tK_expect=0;tK_expect<getLt();tK_expect++){      
      for(int t_expect=0;t_expect<getLt();t_expect++){
	if(!(in >> tK)) error_exit(std::cout << "type1234Data::parse failed to read tK for config " << sample << "\n");
	if(tK != tK_expect) error_exit(std::cout << "type1234Data::parse tK doesn't match expectations: " << tK << ":" << tK_expect << " for config " << sample << "\n");

	if(!(in >> t)) error_exit(std::cout << "type1234Data::parse failed to read t for config " << sample << "\n");
	if(t != t_expect) error_exit(std::cout << "type1234Data::parse t doesn't match expectations: " << t << ":" << t_expect << " for config " << sample << "\n");

	data[tK][t].parse(in,sample);
      }
    }
  }
  void parse(const std::string &filename, const int sample){
    std::cout << "type1234Data parsing " << filename << " to sample " << sample << std::endl;
    std::ifstream is(filename.c_str());
    if(!is.good()){ std::cout << "Could not open file \"" << filename << "\"\n"; std::cout.flush(); exit(-1); }
    parse(is,sample);
    if(is.fail() || is.bad()){ std::cout << "Error reading file \"" << filename << "\"\n"; std::cout.flush(); exit(-1); }
    is.close();
  }  
  inline type1234Data operator+(const type1234Data &r) const{
    type1234Data out(*this);
#pragma omp parallel for
    for(int tt=0;tt<Lt*Lt;tt++){
      int tK = tt / Lt,  t = tt % Lt;
      out.data[tK][t] = out.data[tK][t] + r.data[tK][t];
    }
    return out;
  }
  inline type1234Data operator*(const double r) const{
    type1234Data out(*this);
#pragma omp parallel for
    for(int tt=0;tt<Lt*Lt;tt++){
      int tK = tt / Lt,  t = tt % Lt;
      out.data[tK][t] = out.data[tK][t] * r;
    }
    return out;
  }  
    
};


class srcAvgType1234Data{
  typedef contractions tdata; //data associated with given t, i.e. the contractions
  std::vector<tdata> data; //[t]

  int type;
  int Lt;
  int nsample;
  int ncontract;
  int nextra;
public:
  srcAvgType1234Data(const type1234Data &in): type(in.getType()), Lt(in.getLt()), nsample(in.getNsample()), ncontract(in.getNcontract()), nextra(in.getNextra()),
					      data(in.getLt(), tdata(in.getNsample(),in.getNcontract(),in.getNextra())){    
    for(int t=0;t<Lt;t++){
      int n_nonzero = 0;
      for(int tK=0;tK<Lt;tK++){
	if(!in(tK,t).isZero()){
	  data[t] = data[t] + in(tK,t);
	  ++n_nonzero;
	}
      }
      data[t] = data[t] * (1./n_nonzero);
      if(t==0) std::cout << "srcAvgType1234Data type " << type << " averaged over " << n_nonzero << " values of tK\n";
    }
  }
  inline int getType() const{ return type; }
  inline int getLt() const{ return Lt; }
  inline int getNsample() const{ return nsample; }
  inline int getNcontract() const{ return ncontract; }
  inline int getNextra() const{ return nextra; }

  const contractions & operator()(const int t) const{ return data[t]; }
};




inline std::string typeFile(const int traj, const int type, const int tsep_k_pi, const int tsep_pipi, const std::string &data_dir, const threeMomentum &mom = {0,0,0}){
  std::ostringstream os;
  os << data_dir << "/traj_" << traj << "_type" << type;
  if(type != 4) os << "_deltat_" << tsep_k_pi << "_sep_" << tsep_pipi;
  if(type == 1) os << "_mom" << momStr(mom);
  return os.str();
}

type1234Data readType(const int type, const int traj_start, const int traj_inc, const int traj_lessthan, const int tsep_k_pi, const int tsep_pipi, const int Lt, const std::string &data_dir){
  int nsample = (traj_lessthan - traj_start)/traj_inc;
  if(type == 1){
    type1234Data type1_mom111(1,Lt,nsample);    
    type1234Data type1_mom_1_1_1(1,Lt,nsample);
    type1234Data type1_mom_111(1,Lt,nsample);
    type1234Data type1_mom1_1_1(1,Lt,nsample);

#pragma omp parallel for
    for(int i=0;i<nsample;i++){
      const int traj = traj_start + i*traj_inc;
      type1_mom111.parse(typeFile(traj,type,tsep_k_pi,tsep_pipi,data_dir,{1,1,1}),i);
      type1_mom_1_1_1.parse(typeFile(traj,type,tsep_k_pi,tsep_pipi,data_dir,{-1,-1,-1}),i);
      type1_mom_111.parse(typeFile(traj,type,tsep_k_pi,tsep_pipi,data_dir,{-1,1,1}),i);
      type1_mom1_1_1.parse(typeFile(traj,type,tsep_k_pi,tsep_pipi,data_dir,{1,-1,-1}),i);
    }
    return type1_mom111 * (1.0/8.0) + type1_mom_1_1_1 * (1.0/8.0) + type1_mom_111 * (3.0/8.0) + type1_mom1_1_1 * (3.0/8.0); //Project onto A2
  }else if(type == 2 || type == 3 || type == 4){
    type1234Data typedata(type,Lt,nsample);
#pragma omp parallel for
    for(int i=0;i<nsample;i++){
      const int traj = traj_start + i*traj_inc;
      typedata.parse(typeFile(traj,type,tsep_k_pi,tsep_pipi,data_dir),i);
    }
    if(type == 2 || type == 3) typedata = typedata * 0.5; //correct for missing coefficient
    return typedata;
  }else{
    error_exit(std::cout << "readType invalid type " << type << std::endl);
  }
}

enum LR{ VpA, VmA }; 
rawDataDistributionZD computeLRcontraction(const int cidx, const int i, const LR g1, const int j, const LR g2, const contractions &from){
  //(V+aA)(V+bA) = V bA + aA V
  double a = g1 == VpA ? 1 : -1;
  double b = g2 == VpA ? 1 : -1;
  return b * from.C(cidx)(i,'V',j,'A') + a * from.C(cidx)(i,'A',j,'V');
}

typedef dataSeries<int, rawDataDistributionZD> rawTimeSeriesZD; //time, dist(value)

NumericVector<rawTimeSeriesZD> computeAmplitudeType4(const srcAvgType1234Data &type4_srcavg){
  const int nsample = type4_srcavg.getNsample();
  const int Lt = type4_srcavg.getLt();
  NumericVector<rawTimeSeriesZD> out(10, rawTimeSeriesZD(type4_srcavg.getLt(), [nsample](const int t){ return typename rawTimeSeriesZD::ElementType(t,rawDataDistributionZD(nsample)); }) );
  
  for(int t=0;t<Lt;t++){
#define A(IDX) out[IDX-1].value(t)
#define C(IDX,I,G1,J,G2) computeLRcontraction(IDX-23, I, G1, J, G2, type4_srcavg(t));  
    A(1) = 3.*C(23, 0,VmA,1,VpA) -3.*C(26, 1,VpA,0,VmA);
    A(2) = 3.*C(24, 0,VmA,1,VpA) -3.*C(27, 1,VpA,0,VmA);
    A(3) = 3.*C(23, 0,VmA,1,VpA) +3.*C(23, 0,VmA,0,VmA) -3.*C(26, 1,VpA,0,VmA) - 3.*C(26, 0,VmA,0,VmA) +3.*C(29, 0,VmA,0,VmA) -3.*C(31, 0,VmA,0,VmA);
    A(4) = 3.*C(24, 0,VmA,1,VpA) +3.*C(25, 0,VmA,0,VmA) -3.*C(27, 1,VpA,0,VmA) -3.*C(28, 0,VmA,0,VmA) + 3.*C(30, 0,VmA,0,VmA) -3.*C(32, 0,VmA,0,VmA);
    A(5) = 3.*C(23, 0,VmA,1,VmA) +3.*C(23, 0,VmA,0,VpA) -3.*C(26, 1,VmA,0,VmA) -3.*C(26, 0,VpA,0,VmA) +3.*C(29, 0,VmA,0,VpA) -3.*C(31, 0,VmA,0,VpA);
    A(6) = 3.*C(24, 0,VmA,1,VmA) +3.*C(25, 0,VmA,0,VpA) -3.*C(27, 1,VmA,0,VmA) -3.*C(28, 0,VpA,0,VmA) +3.*C(30, 0,VmA,0,VpA) -3.*C(32, 0,VmA,0,VpA);
    A(7) = 3.*C(23, 0,VmA,1,VmA) -1.5*C(23, 0,VmA,0,VpA) -3.*C(26, 1,VmA,0,VmA) +1.5*C(26, 0,VpA,0,VmA) -1.5*C(29, 0,VmA,0,VpA) + 1.5*C(31, 0,VmA,0,VpA);
    A(8) = 3.*C(24, 0,VmA,1,VmA) -1.5*C(25, 0,VmA,0,VpA) -3.*C(27, 1,VmA,0,VmA) +1.5*C(28, 0,VpA,0,VmA) -1.5*C(30, 0,VmA,0,VpA) +1.5*C(32,0,VmA,0,VpA);
    A(9) = 3.*C(23, 0,VmA,1,VpA) -1.5*C(23, 0,VmA,0,VmA) -3.*C(26, 1,VpA,0,VmA) +1.5*C(26, 0,VmA,0,VmA) -1.5*C(29, 0,VmA,0,VmA) +1.5*C(31,0,VmA,0,VmA);
    A(10) = 3.*C(24, 0,VmA,1,VpA) -1.5*C(25, 0,VmA,0,VmA) -3.*C(27, 1,VpA,0,VmA) +1.5*C(28,0,VmA,0,VmA) -1.5*C(30, 0,VmA,0,VmA) +1.5*C(32, 0,VmA,0,VmA);
#undef C

    for(int i=0;i<10;i++) out[i].value(t) = out[i].value(t) * (1.0/sqrt(6.0));
  }
  return out;
}

typedef doubleJackknifeDistribution<std::complex<double> > doubleJackknifeDistributionZD;
typedef dataSeries<int, doubleJackknifeDistributionZD> doubleJackknifeTimeSeriesZD; 

NumericVector<doubleJackknifeTimeSeriesZD> doubleJackknifeResample(const NumericVector<rawTimeSeriesZD> &v){
  int nelem = v.size();
  int Lt = v[0].size();
  
  NumericVector<doubleJackknifeTimeSeriesZD> out(nelem, doubleJackknifeTimeSeriesZD(Lt));
  for(int i=0;i<nelem;i++)
    for(int t=0;t<Lt;t++){
      out[i].value(t).resample(v[i].value(t));
      out[i].coord(t) = v[i].coord(t);
    }  
  return out;
}



int main(void){
  int Lt = 64;
  int traj_start = 636;
  int traj_inc = 4;
  int traj_lessthan = 640;
  int tsep_k_pi = 10;
  int tsep_pipi = 4;
  std::string data_dir = ".";
  
  type1234Data type1 = readType(1, traj_start, traj_inc, traj_lessthan, tsep_k_pi, tsep_pipi, Lt, data_dir);
  type1234Data type2 = readType(2, traj_start, traj_inc, traj_lessthan, tsep_k_pi, tsep_pipi, Lt, data_dir);
  type1234Data type3 = readType(3, traj_start, traj_inc, traj_lessthan, tsep_k_pi, tsep_pipi, Lt, data_dir);
  type1234Data type4 = readType(4, traj_start, traj_inc, traj_lessthan, tsep_k_pi, tsep_pipi, Lt, data_dir);
  
  bubbleDataAllMomenta raw_bubble_data;
  readBubble(raw_bubble_data, data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan);

  srcAvgType1234Data type4_srcavg(type4);
  NumericVector<rawTimeSeriesZD> A0_type4 = computeAmplitudeType4(type4_srcavg);

  NumericVector<doubleJackknifeTimeSeriesZD> A0_type4_dj = doubleJackknifeResample(A0_type4);

  return 0;
}
