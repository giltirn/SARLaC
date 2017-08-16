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

#include <fit_pipi_gparity/data_containers.h>
#include <fit_pipi_gparity/mom_data_containers.h>
#include <fit_pipi_gparity/read_data.h>

struct contractionData{
  std::vector<rawDataDistributionD> d;

  contractionData(const int nsample): d(8,rawDataDistributionD(nsample,0.)){
    mapping(0,'V',0,'A'); //touch mapping
  }

  void parse(std::istream &in, const int sample){
    double im; //discard imaginary part
    for(int i=0;i<8;i++)
      if(!(in >> d[i].sample(sample) >> im)) error_exit(std::cout << "type1234Data::contractions::contractionData::parse failed to parse contraction data for config " << sample << "\n");
  }
    
  inline int map_off(const int i, const char A, const int j, const char B) const{
    int iA = A == 'V' ? 0 : 1;
    int iB = B == 'V' ? 0 : 1;
    return i + 2*(iA + 2*(j + 2*iB));      
  }
      
  int mapping(const int i, const char A, const int j, const char B) const{
    static int map[16];
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
    return map[map_off(i,A,j,B)];
  }
    
  inline rawDataDistributionD & operator()(const int i, const char A, const int j, const char B){
    return d[mapping(i,A,j,B)];
  }
  inline const rawDataDistributionD & operator()(const int i, const char A, const int j, const char B) const{
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
    for(int i=0;i<8;i++)
      for(int s=0;s<d[i].size();s++)
	if(d[i].sample(s) != 0.) return false;
    return true;
  }
    
};
    
struct contractions{
  std::vector<contractionData> con;
  std::vector<rawDataDistributionD> extra; //F0 g5, -F1 g5

  contractions(){}
  contractions(const int nsample, const int ncontract, const int nextra = 0): con(ncontract, contractionData(nsample)), extra(nextra,rawDataDistributionD(nsample,0.)){}

  //Parse the ncontract * 8 * 2 double entries for a given tK, t, sample
  void parse(std::istream &in, const int sample){
    for(int c=0;c<con.size();c++) con[c].parse(in,sample);

    double im; //discard imaginary part
    for(int e=0;e<extra.size();e++)
      if(!(in >> extra[e].sample(sample) >> im) ) error_exit(std::cout << "type1234Data::contractions::parse failed to parse extra idx " << e << " data for config " << sample << "\n");    
  }

  const contractionData & C(const int i) const{ assert(i>=0 && i<con.size()); return con[i]; }
  const rawDataDistributionD & mix() const{ return extra[0]; } //only use the F0 g5 spin-flavor structure in practise 
  
  inline bool isZero() const{ 
    for(int i=0;i<con.size();i++) if(!con[i].isZero()) return false;

    for(int e=0;e<extra.size();e++)
      for(int s=0;s<extra[e].size();s++)
	if(extra[e].sample(s) != 0.) return false;
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

  const contractions & operator()(const int tK, const int t) const{
    return data[tK][t];
  }

  std::vector<int> getNonZeroKaonTimeslices() const{
    std::vector<int> out;
    for(int tK=0;tK<Lt;tK++){
      bool allzero = true;
      for(int t=0;t<Lt;t++) if(!data[tK][t].isZero()){ allzero=false; break; }
      if(!allzero) out.push_back(tK);
    }
    return out;
  }  
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
#ifdef DAIQIAN_COMPATIBILITY_MODE
    if(type == 2 || type == 3) typedata = typedata * 0.5; //correct for missing coefficient
#endif
    return typedata;
  }else{
    error_exit(std::cout << "readType invalid type " << type << std::endl);
  }
}


NumericTensor<rawDataDistributionD,1> readA2projectedBubble(const int traj_start, const int traj_inc, const int traj_lessthan, const int tsep_pipi, const int Lt, const std::string &data_dir){
  bubbleDataAllMomenta raw_bubble_data;
  readBubble(raw_bubble_data, data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan);

  NumericTensor<rawDataDistributionD,1> out({Lt});
  for(int t=0;t<Lt;t++) out({t}) = ( raw_bubble_data({1,1,1})(t)  + raw_bubble_data({-1,-1,-1})(t)
				     + raw_bubble_data({-1,1,1})(t) + raw_bubble_data({1,-1,-1})(t)
				     + raw_bubble_data({1,-1,1})(t) + raw_bubble_data({-1,1,-1})(t)
				     + raw_bubble_data({1,1,-1})(t) + raw_bubble_data({-1,-1,1})(t) )/8.;
  return out;
}

enum LR{ VpA, VmA }; 
rawDataDistributionD computeLRcontraction(const int cidx, const int i, const LR g1, const int j, const LR g2, const contractions &from){
  //(V+aA)(V+bA) = V bA + aA V
  double a = g1 == VpA ? 1 : -1;
  double b = g2 == VpA ? 1 : -1;
  return b * from.C(cidx)(i,'V',j,'A') + a * from.C(cidx)(i,'A',j,'V');
}

//For a given type, what is the first diagram index
inline int idxOffset(const int type){
  static const int offs[4] = {1,7,13,23}; 
  return offs[type-1];
}

struct computeAmplitudeAlltKtensorControls{
  typedef type1234Data inputType;
  typedef NumericTensor<rawDataDistributionD,3> outputType; //(Qidx,tK,t)

  int Lt;
  outputType out;
  const inputType &in;
  
  computeAmplitudeAlltKtensorControls(const inputType &_in): Lt(_in.getLt()), in(_in), out({10,_in.getLt(),_in.getLt()}){}

  inline int size(){ return in.getLt()*in.getLt(); }
  
  inline rawDataDistributionD & A(const int I, const int tt){ return out({I,tt/Lt, tt%Lt}); }
  
  inline rawDataDistributionD C(const int IDX, const int i, const LR g1, const int j, const LR g2, const int tt){
    int tK=tt/Lt, t=tt%Lt;
    return computeLRcontraction(IDX, i,g1,j,g2, in(tK, t));
  }

  inline void normalize(const double nrm, const int tt){
    int tK=tt/Lt, t=tt%Lt;
    for(int i=0;i<10;i++) out({i,tK,t}) = out({i,tK,t}) * nrm;
  }
};  





template<typename Controls>
typename Controls::outputType computeAmplitudeType1(const typename Controls::inputType &in){
  Controls controls(in);
  int off = idxOffset(1);
    
  for(int t=0;t<controls.size();t++){
#define A(I) controls.A(I-1,t)
#define C(IDX,I,G1,J,G2) controls.C(IDX-off,I,G1,J,G2,t)

    A(1) = C(1, 1,VpA,0,VmA) - C(4, 0,VmA,1,VpA) -2.*C(4, 0,VmA,0,VmA);
    A(2) = C(2, 1,VpA,0,VmA) - C(5, 0,VmA,1,VpA) -2.*C(6, 0,VmA,0,VmA);
    A(3) = -3.*C(4, 0,VmA,1,VpA) -3.*C(4, 0,VmA,0,VmA);
    A(4) = C(3, 0,VmA,0,VmA) -3.*C(6, 0,VmA,0,VmA) +C(2, 1,VpA,0,VmA) -3.*C(5, 0,VmA,1,VpA);      
    A(5) = -3.*C(4, 0,VmA,1,VmA) -3.*C(4, 0,VmA,0,VpA);
    A(6) = C(3, 0,VpA,0,VmA) -3.*C(6, 0,VmA,0,VpA) + C(2, 1,VmA,0,VmA) -3.*C(5, 0,VmA,1,VmA);
    A(7) = C(1, 1,VmA,0,VmA) -0.5*C(1, 0,VpA,0,VmA) -1.5*C(4, 0,VmA,0,VpA);
    A(8) = -0.5*C(3, 0,VpA,0,VmA) -1.5*C(6, 0,VmA,0,VpA) + C(2, 1,VmA,0,VmA);
    A(9) = C(1, 1,VpA,0,VmA) -0.5*C(1, 0,VmA,0,VmA) -1.5*C(4, 0,VmA,0,VmA);
    A(10) = -0.5*C(3, 0,VmA,0,VmA) -1.5*C(6, 0,VmA,0,VmA) + C(2, 1,VpA,0,VmA);
#undef C

    controls.normalize(1./sqrt(6.),t);
  }
  return controls.out;
}

template<typename Controls>
typename Controls::outputType computeAmplitudeType2(const typename Controls::inputType &in){
  Controls controls(in);
  int off = idxOffset(2);
    
  for(int t=0;t<controls.size();t++){
#define A(I) controls.A(I-1,t)
#define C(IDX,I,G1,J,G2) controls.C(IDX-off,I,G1,J,G2,t)

    A(1) = 3.*C(7, 0,VmA,1,VpA) -3.*C(10, 1,VpA,0,VmA);
    A(2) = 3.*C(8, 0,VmA,1,VpA) -3.*C(11, 1,VpA,0,VmA);
    A(3) = 3.*C(7, 0,VmA,1,VpA) +3.*C(7, 0,VmA,0,VmA) -3.*C(10, 1,VpA,0,VmA) -3.*C(10, 0,VmA,0,VmA);
    A(4) = 3.*C(8, 0,VmA,1,VpA) +3.*C(9, 0,VmA,0,VmA) -3.*C(11, 1,VpA,0,VmA) -3.*C(12, 0,VmA,0,VmA);
    A(5) = 3.*C(7, 0,VmA,1,VmA) +3.*C(7, 0,VmA,0,VpA) -3.*C(10, 1,VmA,0,VmA) -3.*C(10, 0,VpA,0,VmA);
    A(6) = 3.*C(8, 0,VmA,1,VmA) +3.*C(9, 0,VmA,0,VpA) -3.*C(11, 1,VmA,0,VmA) -3.*C(12, 0,VpA,0,VmA);
    A(7) = 3.*C(7, 0,VmA,1,VmA) -1.5*C(7, 0,VmA,0,VpA) -3.*C(10, 1,VmA,0,VmA) +1.5*C(10, 0,VpA,0,VmA);
    A(8) = 3.*C(8, 0,VmA,1,VmA) -1.5*C(9, 0,VmA,0,VpA) -3.*C(11, 1,VmA,0,VmA) +1.5*C(12, 0,VpA,0,VmA);
    A(9) = 3.*C(7, 0,VmA,1,VpA) -1.5*C(7, 0,VmA,0,VmA) -3.*C(10, 1,VpA,0,VmA) +1.5*C(10, 0,VmA,0,VmA); 
    A(10) = 3.*C(8, 0,VmA,1,VpA) -1.5*C(9, 0,VmA,0,VmA) -3.*C(11, 1,VpA,0,VmA) +1.5*C(12, 0,VmA,0,VmA);
#undef C

    controls.normalize(1./sqrt(6.),t);
  }
  return controls.out;
}

template<typename Controls>
typename Controls::outputType computeAmplitudeType3(const typename Controls::inputType &in){
  Controls controls(in);
  int off = idxOffset(3);
    
  for(int t=0;t<controls.size();t++){
#define A(I) controls.A(I-1,t)
#define C(IDX,I,G1,J,G2) controls.C(IDX-off,I,G1,J,G2,t)

    A(1) = 3.*C(13, 0,VmA,1,VpA) -3.*C(16, 1,VpA,0,VmA);      
    A(2) = 3.*C(14, 0,VmA,1,VpA) -3.*C(17, 1,VpA,0,VmA);      
    A(3) = 3.*C(13, 0,VmA,1,VpA) + 3.*C(13, 0,VmA,0,VmA) -3.*C(16, 1,VpA,0,VmA) -3.*C(16, 0,VmA,0,VmA) + 3.*C(19, 0,VmA,0,VmA) -3.*C(21, 0,VmA,0,VmA);    
    A(4) = 3.*C(14, 0,VmA,1,VpA) +3.*C(15, 0,VmA,0,VmA) -3.*C(17, 1,VpA,0,VmA) -3.*C(18, 0,VmA,0,VmA) +3.*C(20, 0,VmA,0,VmA) -3.*C(22, 0,VmA,0,VmA); 
    A(5) = 3.*C(13, 0,VmA,1,VmA) + 3.*C(13, 0,VmA,0,VpA) -3.*C(16, 1,VmA,0,VmA) -3.*C(16, 0,VpA,0,VmA) + 3.*C(19, 0,VmA,0,VpA) -3.*C(21, 0,VmA,0,VpA);      
    A(6) = 3.*C(14, 0,VmA,1,VmA) +3.*C(15, 0,VmA,0,VpA) -3.*C(17, 1,VmA,0,VmA) -3.*C(18, 0,VpA,0,VmA) +3.*C(20, 0,VmA,0,VpA) -3.*C(22, 0,VmA,0,VpA);
    A(7) = 3.*C(13, 0,VmA,1,VmA) -1.5*C(13, 0,VmA,0,VpA) -3.*C(16, 1,VmA,0,VmA) +1.5*C(16, 0,VpA,0,VmA) -1.5*C(19, 0,VmA,0,VpA) +1.5*C(21, 0,VmA,0,VpA);                                             
    A(8) = 3.*C(14, 0,VmA,1,VmA) -1.5*C(15, 0,VmA,0,VpA) -3.*C(17, 1,VmA,0,VmA) +1.5*C(18, 0,VpA,0,VmA) -1.5*C(20, 0,VmA,0,VpA) +1.5*C(22, 0,VmA,0,VpA);
    A(9) = 3.*C(13, 0,VmA,1,VpA) -1.5*C(13, 0,VmA,0,VmA) -3.*C(16, 1,VpA,0,VmA) +1.5*C(16, 0,VmA,0,VmA) -1.5*C(19, 0,VmA,0,VmA) +1.5*C(21, 0,VmA,0,VmA);
    A(10) = 3.*C(14, 0,VmA,1,VpA) -1.5*C(15, 0,VmA,0,VmA) -3.*C(17, 1,VpA,0,VmA) +1.5*C(18, 0,VmA,0,VmA) -1.5*C(20, 0,VmA,0,VmA) +1.5*C(22, 0,VmA,0,VmA);
#undef C

    controls.normalize(1./sqrt(6.),t);
  }
  return controls.out;
}


template<typename Controls>
typename Controls::outputType computeAmplitudeType4(const typename Controls::inputType &in){
  Controls controls(in);
  int off = idxOffset(4);
    
  for(int t=0;t<controls.size();t++){
#define A(I) controls.A(I-1,t)
#define C(IDX,I,G1,J,G2) controls.C(IDX-off,I,G1,J,G2,t)

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

    controls.normalize(1./sqrt(6.),t);
  }
  return controls.out;
}


//Functor for tensor reduction by average
template<typename T, int Rank>
struct averageDimensionFunctor{
  const int dim;
  std::vector<int> const* use;
  
  averageDimensionFunctor(const int _dim, std::vector<int> const*_use = NULL): dim(_dim), use(_use){}
  
  void operator()(T &o, int const *coord, const NumericTensor<T,Rank> &from) const{
    int full_coord[Rank];
    int i=0; for(int ii=0;ii<Rank;ii++) if(ii!=dim) full_coord[ii] = coord[i++];    
    zeroit(o);
    if(use != NULL){
      assert(use->size()> 0);
      full_coord[dim] = use->at(0);
      o = from(full_coord);      
      for(int i=1;i<use->size();i++){
	full_coord[dim] = use->at(i);
	o = o + from(full_coord);
      }
      o = o/double(use->size());
    }else{
      assert(from.size(dim)>0);
      full_coord[dim] = 0;
      o = from(full_coord);
      for(int i=1;i<from.size(dim);i++){
	full_coord[dim] = i;
	o = o + from(full_coord);
      }
      o = o/double(from.size(dim));
    }
  }
};
    
template<typename Resampled, typename Raw>
struct resampleFunctor{
  inline Resampled operator()(int const* coord,const Raw &from) const{
    Resampled o(from.size());
    o.resample(from);
    return o;
  }
};

template<typename distributionType>
struct iterate;

template<typename T>
struct iterate<doubleJackknifeDistribution<T> >{
  static inline int size(const doubleJackknifeDistribution<T> &from){ return from.size() * (from.size()-1); } //j + (from.size()-1)*i
  static inline const T& at(const int i, const doubleJackknifeDistribution<T> &from){
    const int nn = from.size()-1;
    return from.sample(i/nn).sample(i%nn);
  }
  static inline T & at(const int i, doubleJackknifeDistribution<T> &from){
    const int nn = from.size()-1;
    return from.sample(i/nn).sample(i%nn);
  }    
};
template<typename T>
struct iterate<rawDataDistribution<T> >{
  static inline int size(const rawDataDistribution<T> &from){ return from.size(); } 
  static inline const T& at(const int i, const rawDataDistribution<T> &from){
    return from.sample(i);
  }
  static inline T & at(const int i, rawDataDistribution<T> &from){
    return from.sample(i);
  }    
};

  

template<typename distributionType>
struct multiplyBubbleFunctor{
  const NumericTensor<distributionType,1> &bubble;
  int tsep_k_pi;
  int tsep_pipi;
  int Lt;
  int tK_dim;
  multiplyBubbleFunctor(const NumericTensor<distributionType,1> &_bubble, const int _tsep_k_pi, const int _tsep_pipi, const int _Lt, const int _tK_dim): bubble(_bubble), tsep_k_pi(_tsep_k_pi), tsep_pipi(_tsep_pipi), Lt(_Lt), tK_dim(_tK_dim){}
  inline distributionType operator()(int const* coord, const distributionType &from) const{
    typedef iterate<distributionType> iter;
    int tK = coord[tK_dim]; int tB = (tK + tsep_k_pi + tsep_pipi) % Lt;
    distributionType out(from);
    for(int i=0;i<iter::size(out);i++)
      iter::at(i,out) = iter::at(i,out) * iter::at(i,bubble({tB}));
    return out;
  }
};



int main(void){
  int Lt = 64;
  int traj_start = 636;
  int traj_inc = 4;
  int traj_lessthan = 648;
  int tsep_k_pi = 10;
  int tsep_pipi = 4;
  std::string data_dir = ".";
  
  type1234Data type1 = readType(1, traj_start, traj_inc, traj_lessthan, tsep_k_pi, tsep_pipi, Lt, data_dir);
  type1234Data type2 = readType(2, traj_start, traj_inc, traj_lessthan, tsep_k_pi, tsep_pipi, Lt, data_dir);
  type1234Data type3 = readType(3, traj_start, traj_inc, traj_lessthan, tsep_k_pi, tsep_pipi, Lt, data_dir);
  type1234Data type4 = readType(4, traj_start, traj_inc, traj_lessthan, tsep_k_pi, tsep_pipi, Lt, data_dir);

  std::vector<int> type1_nonzerotK = type1.getNonZeroKaonTimeslices();
  std::vector<int> type2_nonzerotK = type2.getNonZeroKaonTimeslices();
  std::vector<int> type3_nonzerotK = type3.getNonZeroKaonTimeslices();
  std::vector<int> type4_nonzerotK = type4.getNonZeroKaonTimeslices();

  NumericTensor<rawDataDistributionD,1> bubble = readA2projectedBubble(traj_start,traj_inc,traj_lessthan,tsep_pipi,Lt,data_dir);
  NumericTensor<doubleJackknifeDistributionD,1> bubble_dj = bubble.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());

  //Get the type 4 contribution without the bubble and double-jackknife resample
  NumericTensor<rawDataDistributionD,3> A0_type4_alltK = computeAmplitudeType4<computeAmplitudeAlltKtensorControls>(type4); //[Qidx][tK][t]
  NumericTensor<doubleJackknifeDistributionD,3> A0_type4_alltK_dj = A0_type4_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());

  //Get and double-jackknife resample the mix4 (no bubble) diagram
  NumericTensor<rawDataDistributionD,2> mix4_alltK({Lt,Lt}); //[tK][t]
  for(int tK=0;tK<Lt;tK++) for(int t=0;t<Lt;t++) mix4_alltK({tK,t}) = type4(tK,t).mix();
  NumericTensor<doubleJackknifeDistributionD,2> mix4_alltK_dj = mix4_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());
  
  //tK-average the type4 and mix4 (no bubble) double-jackknife data then double-jackknife resample  
  NumericTensor<doubleJackknifeDistributionD,2> A0_type4_srcavg_dj = A0_type4_alltK_dj.reduce(1, averageDimensionFunctor<doubleJackknifeDistributionD,3>(1, &type4_nonzerotK)); //[Qidx][t]
  NumericTensor<doubleJackknifeDistributionD,1> mix4_srcavg_dj = mix4_alltK_dj.reduce(0, averageDimensionFunctor<doubleJackknifeDistributionD,2>(0, &type4_nonzerotK)); //[t]
   
  //Compute alpha from the above
  NumericTensor<doubleJackknifeDistributionD,2> alpha({10,Lt}, [&](int const *c){ return A0_type4_srcavg_dj({c[0],c[1]})/mix4_srcavg_dj({c[1]}); }); //[Qidx][t]

  //Compute vacuum subtractions
  NumericTensor<doubleJackknifeDistributionD,3> A0_type4_alltK_vacsub = A0_type4_alltK_dj.transform(multiplyBubbleFunctor<doubleJackknifeDistributionD>(bubble_dj,tsep_k_pi,tsep_pipi,Lt,1));//[Qidx][tK][t]
  NumericTensor<doubleJackknifeDistributionD,2> mix4_alltK_vacsub = mix4_alltK_dj.transform(multiplyBubbleFunctor<doubleJackknifeDistributionD>(bubble_dj,tsep_k_pi,tsep_pipi,Lt,0)); //[tK][t]

  //Source average the vacuum subtractions
  NumericTensor<doubleJackknifeDistributionD,2> A0_type4_srcavg_vacsub = A0_type4_alltK_vacsub.reduce(1, averageDimensionFunctor<doubleJackknifeDistributionD,3>(1, &type4_nonzerotK)); //[Qidx][t]
  NumericTensor<doubleJackknifeDistributionD,1> mix4_srcavg_vacsub = mix4_alltK_vacsub.reduce(0, averageDimensionFunctor<doubleJackknifeDistributionD,2>(0, &type4_nonzerotK)); //[t]
  
  //Apply the bubble diagram to the type4 and mix4 data
  A0_type4_alltK = A0_type4_alltK.transform(multiplyBubbleFunctor<rawDataDistributionD>(bubble,tsep_k_pi,tsep_pipi,Lt,1));
  mix4_alltK = mix4_alltK.transform(multiplyBubbleFunctor<rawDataDistributionD>(bubble,tsep_k_pi,tsep_pipi,Lt,0));

  //Recompute double-jackknife resamplings of full type4 and mix4
  A0_type4_alltK_dj = A0_type4_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());
  mix4_alltK_dj = mix4_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());

  //Recompute double-jackknife resamplings of full, source-averaged type4 and mix4
  A0_type4_srcavg_dj = A0_type4_alltK_dj.reduce(1, averageDimensionFunctor<doubleJackknifeDistributionD,3>(1, &type4_nonzerotK)); //[Qidx][t]
  mix4_srcavg_dj = mix4_alltK_dj.reduce(0, averageDimensionFunctor<doubleJackknifeDistributionD,2>(0, &type4_nonzerotK)); //[t]

  //Get, double-jackknife resample and source-average the type1, type2, type3 and mix3 contributions
  NumericTensor<rawDataDistributionD,3> A0_type1_alltK = computeAmplitudeType1<computeAmplitudeAlltKtensorControls>(type1); //[Qidx][tK][t]
  NumericTensor<doubleJackknifeDistributionD,3> A0_type1_alltK_dj = A0_type1_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());
  NumericTensor<doubleJackknifeDistributionD,2> A0_type1_srcavg_dj = A0_type1_alltK_dj.reduce(1, averageDimensionFunctor<doubleJackknifeDistributionD,3>(1, &type1_nonzerotK)); //[Qidx][t]
  
  NumericTensor<rawDataDistributionD,3> A0_type2_alltK = computeAmplitudeType2<computeAmplitudeAlltKtensorControls>(type2); //[Qidx][tK][t]
  NumericTensor<doubleJackknifeDistributionD,3> A0_type2_alltK_dj = A0_type2_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());
  NumericTensor<doubleJackknifeDistributionD,2> A0_type2_srcavg_dj = A0_type2_alltK_dj.reduce(1, averageDimensionFunctor<doubleJackknifeDistributionD,3>(1, &type2_nonzerotK)); //[Qidx][t]

  NumericTensor<rawDataDistributionD,3> A0_type3_alltK = computeAmplitudeType3<computeAmplitudeAlltKtensorControls>(type3); //[Qidx][tK][t]
  NumericTensor<doubleJackknifeDistributionD,3> A0_type3_alltK_dj = A0_type3_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());
  NumericTensor<doubleJackknifeDistributionD,2> A0_type3_srcavg_dj = A0_type3_alltK_dj.reduce(1, averageDimensionFunctor<doubleJackknifeDistributionD,3>(1, &type3_nonzerotK)); //[Qidx][t]

  NumericTensor<rawDataDistributionD,2> mix3_alltK({Lt,Lt}); //[tK][t]
  for(int tK=0;tK<Lt;tK++) for(int t=0;t<Lt;t++) mix3_alltK({tK,t}) = type3(tK,t).mix();
  NumericTensor<doubleJackknifeDistributionD,2> mix3_alltK_dj = mix3_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());
  NumericTensor<doubleJackknifeDistributionD,1> mix3_srcavg_dj = mix3_alltK_dj.reduce(0, averageDimensionFunctor<doubleJackknifeDistributionD,2>(0, &type3_nonzerotK)); //[t]
  
  //Subtract the pseudoscalar operators and mix4 vacuum term
  A0_type3_srcavg_dj = A0_type3_srcavg_dj.transform([&](int const* coord, const doubleJackknifeDistributionD &from){ return doubleJackknifeDistributionD(from - alpha(coord)*mix3_srcavg_dj(coord+1)); }); 
  A0_type4_srcavg_dj = A0_type4_srcavg_dj.transform(
						    [&](int const* coord, const doubleJackknifeDistributionD &from){
						      return doubleJackknifeDistributionD(from - alpha(coord)*( mix4_srcavg_dj(coord+1) + mix4_srcavg_vacsub(coord+1) ) );
						    }
						    ); 
  
  //Perform the type 4 vacuum subtraction
  A0_type4_srcavg_dj = A0_type4_srcavg_dj - A0_type4_srcavg_vacsub;

  //Get the full double-jackknife amplitude
  NumericTensor<doubleJackknifeDistributionD,2> A0_full_srcavg_dj = A0_type1_srcavg_dj + A0_type2_srcavg_dj + A0_type3_srcavg_dj + A0_type4_srcavg_dj;
  
  //Get the single-elimination jackknife distributions
  NumericTensor<jackknifeDistributionD,2> A0_type1_srcavg_j = A0_type1_srcavg_dj.transform([](int const* coord, const doubleJackknifeDistributionD &from){ return from.toJackknife(); });
  NumericTensor<jackknifeDistributionD,2> A0_type2_srcavg_j = A0_type2_srcavg_dj.transform([](int const* coord, const doubleJackknifeDistributionD &from){ return from.toJackknife(); });
  NumericTensor<jackknifeDistributionD,2> A0_type3_srcavg_j = A0_type3_srcavg_dj.transform([](int const* coord, const doubleJackknifeDistributionD &from){ return from.toJackknife(); });
  NumericTensor<jackknifeDistributionD,2> A0_type4_srcavg_j = A0_type4_srcavg_dj.transform([](int const* coord, const doubleJackknifeDistributionD &from){ return from.toJackknife(); });

  NumericTensor<jackknifeDistributionD,2> A0_full_srcavg_j = A0_type1_srcavg_j + A0_type2_srcavg_j + A0_type3_srcavg_j + A0_type4_srcavg_j;

  NumericTensor<jackknifeDistributionD,2> sigma = A0_full_srcavg_dj.transform([](int const *c, const doubleJackknifeDistributionD &d){ return jackknifeDistributionD(sqrt(doubleJackknifeDistributionD::covariance(d,d))); });

  

  



    
  
  return 0;
}


