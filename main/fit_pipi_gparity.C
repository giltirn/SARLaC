#include <fstream>
#include <algorithm>
#include <minimizer.h>
#include <cost_function.h>
#include <fitfunc.h>
#include <fitfunc_mapping.h>
#include <random.h>
#include <plot.h>
#include <distribution.h>
#include <data_series.h>
#include <parser.h>
#include <common_defs.h>
#include <sstream>

typedef std::complex<double> complexD;

class bubbleData{
  NumericVector<distribution<complexD> > d; //(t).sample(cfg)
  int Lt;
public:
  void setup(const int _Lt, const int _Nsample){ Lt = _Lt; d.resize(_Lt, distribution<complexD>(_Nsample)); }
  
  inline const int getLt() const{ return Lt; }
  inline const int getNsample() const{ return d[0].size(); }
  
  inline distribution<complexD> & operator()(const int t){ return d[t]; }
  inline const distribution<complexD> & operator()(const int t) const { return d[t]; }

  void parse(std::istream &in, const int sample){
    int t;
    for(int t_expect=0;t_expect<Lt;t_expect++){
      if(!(in >> t)) error_exit(std::cout << "bubbleData::parse failed to read t for config " << sample << "\n");
      if(t != t_expect) error_exit(std::cout << "bubbleData::parse t doesn't match expectations: " << t << ":" << t_expect << " for config " << sample << "\n");

      double &re = reinterpret_cast<double(&)[2]>( d[t].sample(sample) )[0];
      double &im = reinterpret_cast<double(&)[2]>( d[t].sample(sample) )[1];
      if(!(in >> re >> im)) error_exit(std::cout << "bubbleData::parse failed to real values for config " << sample << "\n");
    }
  }
  void parse(const std::string &filename, const int sample){
    std::ifstream is(filename.c_str());
    if(!is.good()){ std::cout << "Could not open file \"" << filename << "\"\n"; std::cout.flush(); exit(-1); }
    parse(is,sample);
    if(is.fail() || is.bad()){ std::cout << "Error reading file \"" << filename << "\"\n"; std::cout.flush(); exit(-1); }
    is.close();
  }    
  
};

class figureData{
  NumericMatrix<distribution<complexD> > d; //(tsrc,tsep).sample(cfg)
  int Lt;
  friend std::ostream & operator<<(std::ostream &os, const figureData &f);
public:
  
  typedef figureData ET_tag;
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,figureData>::value, int>::type = 0>
  figureData(U&& expr);

  figureData() = default;
  figureData(const figureData &r) = default;
  figureData(figureData &&r) = default;
  figureData(const int _Lt, const int _Nsample): Lt(_Lt), d(_Lt,distribution<complexD>(_Nsample)){}
  
  figureData &operator=(const figureData &r) = default;
  figureData &operator=(figureData &&r) = default;

  void zero(){
    for(int i=0;i<Lt;i++)
      for(int j=0;j<Lt;j++)
	for(int s=0;s<d(i,j).size();s++)
	  d(i,j).sample(s) = complexD(0.);
  }
    
  
  void setup(const int _Lt, const int _Nsample){ Lt = _Lt; d.resize(_Lt, distribution<complexD>(_Nsample)); }
  
  void parseCDR(std::istream &in, const int sample){
    const int nelems = Lt*Lt;

    int tsrc,tsep;
    for(int e=0;e<nelems;e++){
      int tsep_expect = e % Lt;
      int tsrc_expect = e / Lt;

      if(!(in >> tsrc >> tsep)) error_exit(std::cout << "FigureData::parseCDR failed to read tsrc, tsep for config " << sample << "\n");
      if(tsep != tsep_expect || tsrc != tsrc_expect) error_exit(std::cout << "FigureData tsrc tsep don't match expectations: "
								<< tsrc << ":" << tsrc_expect << " " << tsep << ":" << tsep_expect
								<< " for config " << sample << "\n");
      double &re = reinterpret_cast<double(&)[2]>( d(tsrc,tsep).sample(sample) )[0];
      double &im = reinterpret_cast<double(&)[2]>( d(tsrc,tsep).sample(sample) )[1];
      if(!(in >> re >> im)) error_exit(std::cout << "FigureData::parseCDR failed to real values for config " << sample << "\n");
    }
  }
  void parseCDR(const std::string &filename, const int sample){
    std::ifstream is(filename.c_str());
    if(!is.good()){ std::cout << "Could not open file \"" << filename << "\"\n"; std::cout.flush(); exit(-1); }
    parseCDR(is,sample);
    if(is.fail() || is.bad()){ std::cout << "Error reading file \"" << filename << "\"\n"; std::cout.flush(); exit(-1); }
    is.close();
  }    
    
  //B(tsrc + tsep + tsep_pipi, -p1_snk) B(tsrc, p1_src)
  void computeV(const bubbleData & Bmp1_snk, const bubbleData & Bp1_src, const int tsep_pipi){
    for(int tsrc=0;tsrc<Lt;tsrc++)
      for(int tsep=0;tsep<Lt;tsep++)
	d(tsrc,tsep) = Bmp1_snk( (tsrc + tsep + tsep_pipi) % Lt ) * Bp1_src( tsrc );
  }
  
  inline const int getLt() const{ return Lt; }
  inline const int getNsample() const{ return d(0,0).size(); }
  
  inline distribution<complexD> & operator()(const int tsrc, const int tsep){ return d(tsrc,tsep); }
  inline const distribution<complexD> & operator()(const int tsrc, const int tsep) const { return d(tsrc,tsep); }

  //Data not measured on every tsrc usually
  bool isZero(const int tsrc) const{    
    for(int tsep=0;tsep<Lt;tsep++)
      for(int sample=0;sample<getNsample();sample++){
	const complexD & v = d(tsrc,tsep).sample(sample);
	if(v.real() != 0.0 || v.imag() != 0.0) return false;
      }
    return true;
  }
};

template<>
struct getElem<figureData>{
  static inline auto elem(figureData &v, const int i)->decltype(v(0,0)){ return v(i/v.getLt(), i%v.getLt()); }
  static inline auto elem(const figureData &v, const int i)->decltype(v(0,0)){ return v(i/v.getLt(), i%v.getLt()); }
  static inline int common_properties(const figureData &v){ return v.getLt(); }
};

template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, figureData::ET_tag>::value && !std::is_same<U,figureData>::value, int>::type>
figureData::figureData(U&& expr): d(expr.common_properties()), Lt(expr.common_properties()){
  for(int i=0;i<Lt*Lt;i++)
    getElem<figureData>::elem(*this, i) = expr[i];
}


std::ostream & operator<<(std::ostream &os, const figureData &f){
  basicPrint<> p(os);
  p << f.d;
  return os;
}



typedef std::array<int,3> threeMomentum;
typedef std::pair<threeMomentum, threeMomentum> sinkSourceMomenta;


template<typename DataPolicies>
class figureDataAllMomentaBase{
public:
  typedef typename DataPolicies::ContainerType ContainerType;
  typedef typename std::map<sinkSourceMomenta, ContainerType>::const_iterator const_iterator;  
private:
  typedef std::map<sinkSourceMomenta, ContainerType> MapType;
  MapType C;
  MapType D;
  MapType R;
  MapType V;
  int Lt;
  int Nsample;
  
  ContainerType & get(const char fig, const sinkSourceMomenta &mom, bool lock){
    MapType* mp;
    switch(fig){
    case 'C':
      mp = &C; break;
    case 'D':
      mp = &D; break;
    case 'R':
      mp = &R; break;
    case 'V':
      mp = &V; break;
    default:
      error_exit(std::cout << "figureDataAllMomenta::get invalid figure " << fig << std::endl);
    }
    typename MapType::iterator it = mp->find(mom);
    if(it == mp->end()){
      if(lock) error_exit(std::cout << "figureDataAllMomenta::get Could not find requested momentum\n");

      it = mp->insert(std::make_pair(mom, figureData())).first;
      DataPolicies::setup(it->second, Lt, Nsample);
    }
    return it->second;
  }
public:
  figureDataAllMomentaBase(const int _Lt, const int _Nsample): Lt(_Lt), Nsample(_Nsample){}
  figureDataAllMomentaBase(){}

  void setup(const int _Lt, const int _Nsample){
    Lt = _Lt; Nsample = _Nsample;
  }
  
  const ContainerType &operator()(const char fig, const sinkSourceMomenta &mom) const{
    figureDataAllMomentaBase<DataPolicies> *t = const_cast<figureDataAllMomentaBase<DataPolicies> *>(this);
    return const_cast<const ContainerType &>( t->get(fig,mom,true) );
  }
  ContainerType &operator()(const char fig, const sinkSourceMomenta &mom){
    return this->get(fig,mom,false);
  }
  inline int getLt() const{ return Lt; }
  inline int getNsample() const{ return Nsample; }
};


struct doubleJackDataPolicy{
  typedef std::vector<doubleJackknifeDistribution<complexD> > ContainerType;
  static void setup(ContainerType &con, const int Lt, const int Nsample){ con.resize(Lt, doubleJackknifeDistribution<complexD>(Nsample)); }
};
struct figureDataPolicy{
  typedef figureData ContainerType;
  static void setup(ContainerType &con, const int Lt, const int Nsample){ con.setup(Lt,Nsample); }
};

typedef figureDataAllMomentaBase<figureDataPolicy> figureDataAllMomenta;
typedef figureDataAllMomentaBase<doubleJackDataPolicy> figureDataDoubleJackAllMomenta;



// class figureDataAllMomenta{
//   std::map<sinkSourceMomenta, figureData> C;
//   std::map<sinkSourceMomenta, figureData> D;
//   std::map<sinkSourceMomenta, figureData> R;
//   std::map<sinkSourceMomenta, figureData> V;
//   int Lt;
//   int Nsample;
  
//   figureData & get(const char fig, const sinkSourceMomenta &mom, bool lock){
//     std::map<sinkSourceMomenta, figureData>* mp;
//     switch(fig){
//     case 'C':
//       mp = &C; break;
//     case 'D':
//       mp = &D; break;
//     case 'R':
//       mp = &R; break;
//     case 'V':
//       mp = &V; break;
//     default:
//       error_exit(std::cout << "figureDataAllMomenta::get invalid figure " << fig << std::endl);
//     }
//     std::map<sinkSourceMomenta, figureData>::iterator it = mp->find(mom);
//     if(it == mp->end()){
//       if(lock) error_exit(std::cout << "figureDataAllMomenta::get Could not find requested momentum\n");

//       it = mp->insert(std::make_pair(mom, figureData())).first;
//       it->second.setup(Lt,Nsample);
//     }
//     return it->second;
//   }
// public:
//   figureDataAllMomenta(const int _Lt, const int _Nsample): Lt(_Lt), Nsample(_Nsample){}
//   figureDataAllMomenta(){}

//   void setup(const int _Lt, const int _Nsample){
//     Lt = _Lt; Nsample = _Nsample;
//   }
  
//   const figureData &operator()(const char fig, const sinkSourceMomenta &mom) const{
//     figureDataAllMomenta *t = const_cast<figureDataAllMomenta *>(this);
//     return const_cast<const figureData &>( t->get(fig,mom,true) );
//   }
//   figureData &operator()(const char fig, const sinkSourceMomenta &mom){
//     return this->get(fig,mom,false);
//   }
//   inline int getLt() const{ return Lt; }
//   inline int getNsample() const{ return Nsample; }
// };

template<typename DataPolicies>
class bubbleDataAllMomentaBase{
public:
  typedef typename DataPolicies::ContainerType ContainerType;
  typedef typename std::map<threeMomentum, ContainerType>::const_iterator const_iterator;  
private:
  std::map<threeMomentum, ContainerType> B;
  int Lt;
  int Nsample;
  
  ContainerType & get(const threeMomentum &mom, bool lock){
    typename std::map<threeMomentum, ContainerType>::iterator it = B.find(mom);
    if(it == B.end()){
      if(lock) error_exit(std::cout << "bubbleDataAllMomenta::get Could not find requested momentum\n");

      it = B.insert(std::make_pair(mom, ContainerType())).first;
      DataPolicies::setup(it->second,Lt,Nsample);
    }
    return it->second;
  }
public:
  bubbleDataAllMomentaBase(const int _Lt, const int _Nsample): Lt(_Lt), Nsample(_Nsample){}
  bubbleDataAllMomentaBase(){}

  void setup(const int _Lt, const int _Nsample){
    Lt = _Lt; Nsample = _Nsample;
  }
  
  const ContainerType &operator()(const threeMomentum &mom) const{
    bubbleDataAllMomentaBase<DataPolicies> *t = const_cast<bubbleDataAllMomentaBase<DataPolicies> *>(this);
    return const_cast<const ContainerType &>( t->get(mom,true) );
  }
  ContainerType &operator()(const threeMomentum &mom){
    return this->get(mom,false);
  }
  inline int getLt() const{ return Lt; }
  inline int getNsample() const{ return Nsample; }

  inline const_iterator begin() const{ return B.begin(); }
  inline const_iterator end() const{ return B.end(); }
};

struct bubbleDataPolicy{
  typedef bubbleData ContainerType;
  static void setup(ContainerType &con, const int Lt, const int Nsample){ con.setup(Lt,Nsample); }
};

typedef bubbleDataAllMomentaBase<bubbleDataPolicy> bubbleDataAllMomenta;
typedef bubbleDataAllMomentaBase<doubleJackDataPolicy> bubbleDataDoubleJackAllMomenta;


inline sinkSourceMomenta momComb(const int snkx, const int snky, const int snkz,
				 const int srcx, const int srcy, const int srcz){
  return sinkSourceMomenta( {snkx,snky,snkz}, {srcx,srcy,srcz} );
}
inline sinkSourceMomenta momComb(const threeMomentum &snk, const threeMomentum &src){
  return sinkSourceMomenta(snk,src);
}


template<typename T, typename Tag>
struct tagged{
  T value;
  inline operator T() const{ return value; }
  inline tagged(const T &v): value(v){}
  inline tagged(){}
};

template<typename DistributionType>
class correlationFunction: public dataSeries< tagged<double,correlationFunction<DistributionType> > , DistributionType>{
  typedef dataSeries< tagged<double,correlationFunction<DistributionType> > , DistributionType> Parent;
public:
  typedef typename Parent::ElementType ElementType;
  typedef correlationFunction<DistributionType> ET_tag;
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,correlationFunction<DistributionType> >::value, int>::type = 0>
  correlationFunction(U&& expr) : Parent(expr.common_properties()){
   for(int i=0;i<this->size();i++)
     (*this)[i] = expr[i];
  }

  correlationFunction(correlationFunction &&r) = default;
  template<typename Initializer>
  inline correlationFunction(const int n, const Initializer &initializer): Parent(n,initializer){}
};

template<typename DistributionType>
struct getElem<correlationFunction<DistributionType> >{
  static inline auto elem(correlationFunction<DistributionType> &v, const int i)->decltype(v[0]){ return v[i]; }
  static inline auto elem(const correlationFunction<DistributionType> &v, const int i)->decltype(v[0]){ return v[i]; }
  static inline int common_properties(const correlationFunction<DistributionType> &v){
    return v.size();
  }
};

correlationFunction<distribution<double> > realSourceAverage(const figureData & data){
  int Lt = data.getLt();
  int nsample = data.getNsample();
  correlationFunction<distribution<double> > into(data.getLt(),
						  [nsample](int i) {  return correlationFunction<distribution<double> >::ElementType(i, distribution<double>(nsample,0.)); }
						  );
  std::vector<int> tsrc_include;
  for(int tsrc=0;tsrc<Lt;tsrc++){
    bool is_nonzero = !data.isZero(tsrc);
    if(is_nonzero)
      tsrc_include.push_back(tsrc);
  }
  double N(tsrc_include.size());
    
  for(int tsep=0;tsep<Lt;tsep++){
    for(int s=0;s<data.getNsample();s++){
      double &v = into.value(tsep).sample(s);
      v = 0.;
      for(int i=0;i<tsrc_include.size();i++)
	v += data(tsrc_include[i],tsep).sample(s).real();
      v /= N;
    }
  }
  return into;
}

typedef std::pair<tagged<double, correlationFunction<distribution<double> > >, distribution<double> >  CFDpair;

inline CFDpair operator*(const int a, const CFDpair &e){
  return CFDpair(e.first, a*e.second);
}
inline CFDpair operator+(const CFDpair &d, const CFDpair &e){
  assert(e.first == d.first);
  return CFDpair(e.first, d.second+e.second);
}
inline CFDpair operator-(const CFDpair &d, const CFDpair &e){
  assert(e.first == d.first);
  return CFDpair(e.first, d.second-e.second);
}

figureData projectA2(const char fig, const figureDataAllMomenta &raw_data){
  threeMomentum R[8] = { {1,1,1}, {-1,-1,-1},
			 {1,1,-1}, {-1,-1,1},
			 {1,-1,1}, {-1,1,-1},
			 {-1,1,1}, {1,-1,-1} };
  
  figureData out(raw_data.getLt(), raw_data.getNsample()); out.zero();
  for(int psnk=0;psnk<8;psnk++)
    for(int psrc=0;psrc<8;psrc++)
      out = out + raw_data(fig, momComb(R[psnk], R[psrc]));
  out = out / double(64);
  return out;
}

inline std::string momStr(const threeMomentum &p){
  std::ostringstream os;
  for(int i=0;i<3;i++)
    os << (p[i] < 0 ? "_" : "") << abs(p[i]);
  return os.str();
}

std::string figureFile(const std::string &data_dir, const char fig, const int traj, const threeMomentum &psnk, const threeMomentum &psrc, const int tsep_pipi){
  std::ostringstream os;
  os << data_dir << "/traj_" << traj << "_Figure" << fig << "_sep" << tsep_pipi << "_mom" << momStr(psrc) << "_mom" << momStr(psnk);
  return os.str();
}

void readFigure(figureDataAllMomenta &raw_data, const char fig, const std::string &data_dir, const int tsep_pipi, const int Lt,
		 const int traj_start, const int traj_inc, const int traj_lessthan){
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
	std::string filename = figureFile(data_dir, fig, traj, R[psnk], R[psrc], tsep_pipi);
	std::cout << "Parsing " << filename << std::endl;
	into.parseCDR(filename, sample);
      }
    }
}

std::string bubbleFile(const std::string &data_dir, const int traj, const threeMomentum &p, const int tsep_pipi){
  std::ostringstream os;
  os << data_dir << "/traj_" << traj << "_FigureVdis_sep" << tsep_pipi << "_mom" << momStr(p);
  return os.str();
}

void readBubble(bubbleDataAllMomenta &raw_data, const std::string &data_dir, const int tsep_pipi, const int Lt,
		 const int traj_start, const int traj_inc, const int traj_lessthan){
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
      std::string filename = bubbleFile(data_dir, traj, R[p], tsep_pipi);
      std::cout << "Parsing " << filename << std::endl;
      into.parse(filename, sample);
    }
  }
}

inline threeMomentum operator-(const threeMomentum &p){
  return threeMomentum({-p[0],-p[1],-p[2]});
}
// template<typename OutputType, typename InputType>
// void computeV(OutputType &raw_data, const InputType &raw_bubble_data, const int tsep_pipi){
//   threeMomentum R[8] = { {1,1,1}, {-1,-1,-1},
// 			 {1,1,-1}, {-1,-1,1},
// 			 {1,-1,1}, {-1,1,-1},
// 			 {-1,1,1}, {1,-1,-1} };
//   for(int psnk=0;psnk<8;psnk++)
//     for(int psrc=0;psrc<8;psrc++){  
//       //B(tsrc + tsep + tsep_pipi, -p1_snk) B(tsrc, p1_src)
//       auto &Bmp1_snk = raw_bubble_data( -R[psnk] );
//       auto &Bp1_src  = raw_bubble_data(  R[psrc] );
//       raw_data('V',momComb(R[psnk],R[psrc])).computeV(Bmp1_snk, Bp1_src, tsep_pipi);
//     }      
// }


void computeV(figureDataAllMomenta &raw_data, const bubbleDataAllMomenta &raw_bubble_data, const int tsep_pipi){
  threeMomentum R[8] = { {1,1,1}, {-1,-1,-1},
			 {1,1,-1}, {-1,-1,1},
			 {1,-1,1}, {-1,1,-1},
			 {-1,1,1}, {1,-1,-1} };
  for(int psnk=0;psnk<8;psnk++)
    for(int psrc=0;psrc<8;psrc++){  
      //B(tsrc + tsep + tsep_pipi, -p1_snk) B(tsrc, p1_src)
      const bubbleData Bmp1_snk = raw_bubble_data( -R[psnk] );
      const bubbleData Bp1_src  = raw_bubble_data(  R[psrc] );
      raw_data('V',momComb(R[psnk],R[psrc])).computeV(Bmp1_snk, Bp1_src, tsep_pipi);
    }      
}

bubbleDataDoubleJackAllMomenta doubleJackknifeResampleBubble(const bubbleDataAllMomenta &bubbles){
  int Lt = bubbles.getLt();
  bubbleDataDoubleJackAllMomenta out(Lt,bubbles.getNsample());
  for(auto it = bubbles.begin(); it != bubbles.end(); it++){
    const threeMomentum & mom = it->first;
    const typename bubbleDataAllMomenta::ContainerType & raw = it->second;
    for(int t=0;t<Lt;t++) out(mom)[t].resample(raw(t));
  }
  return out;
}


struct pooh;
int main(int argc, char* argv[]){
  const std::string data_dir = "/home/ckelly/CPS/build/CPSfit/pipi_data";
  const int tsep_pipi = 4;
  const int Lt = 64;
  const int traj_start = 988;
  const int traj_inc = 4;
  const int traj_lessthan = 996;
  const int nsample = (traj_lessthan - traj_start)/traj_inc;
  
  figureDataAllMomenta raw_data;
  readFigure(raw_data, 'C', data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan);
  readFigure(raw_data, 'D', data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan);
  readFigure(raw_data, 'R', data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan);
  
  bubbleDataAllMomenta raw_bubble_data;
  readBubble(raw_bubble_data, data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan);
  computeV(raw_data, raw_bubble_data, tsep_pipi);

  figureData A2_C = projectA2('C', raw_data);
  figureData A2_D = projectA2('D', raw_data);
  figureData A2_R = projectA2('R', raw_data);
  figureData A2_V = projectA2('V', raw_data);

  auto A2_realavg_C = realSourceAverage(A2_C);
  auto A2_realavg_D = realSourceAverage(A2_D);
  auto A2_realavg_R = realSourceAverage(A2_R);
  auto A2_realavg_V = realSourceAverage(A2_V);

  correlationFunction<distribution<double> > pipi_raw = 2*A2_realavg_D + A2_realavg_C - 6*A2_realavg_R + 3*A2_realavg_V;
  
  correlationFunction<doubleJackknifeDistribution<double> > pipi_dj(Lt,
								    [&pipi_raw,nsample](const int t)
								    {
								      typename correlationFunction<doubleJackknifeDistribution<double> >::ElementType out(t, doubleJackknifeDistribution<double>(nsample));
								      out.second.resample(pipi_raw.value(t));
								      return out;
								    }
								    );
  bubbleDataDoubleJackAllMomenta dj_bubble_data = doubleJackknifeResampleBubble(raw_bubble_data);



  
  return 0;
}
