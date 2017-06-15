#ifndef PIPI_DATA_CONTAINERS_H
#define PIPI_DATA_CONTAINERS_H

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


#endif
