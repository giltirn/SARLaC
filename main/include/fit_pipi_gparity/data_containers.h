#ifndef PIPI_DATA_CONTAINERS_H
#define PIPI_DATA_CONTAINERS_H

typedef std::complex<double> complexD;

struct null_type{};

template<typename DistributionType, typename Policies = null_type>
class bubbleDataBase: public Policies{
  NumericVector<DistributionType> d; //(t).sample(cfg)
  int Lt;
public:
  void setup(const int _Lt, const int _Nsample){ Lt = _Lt; d.resize(_Lt, DistributionType(_Nsample)); }
  
  inline const int getLt() const{ return Lt; }
  inline const int getNsample() const{ return d[0].size(); }
  
  inline DistributionType & operator()(const int t){ return d[t]; }
  inline const DistributionType & operator()(const int t) const { return d[t]; }

  inline DistributionType & at(const int t){ return d[t]; }
  inline const DistributionType & at(const int t) const { return d[t]; }  
};

class bubbleDataPolicies{
  inline bubbleDataBase<distribution<complexD>, bubbleDataPolicies> & upcast(){ return *static_cast< bubbleDataBase<distribution<complexD>, bubbleDataPolicies>* >(this); }
  inline const bubbleDataBase<distribution<complexD>, bubbleDataPolicies> & upcast() const{ return *static_cast< bubbleDataBase<distribution<complexD>, bubbleDataPolicies> const* >(this); }

public:
  void parse(std::istream &in, const int sample){
    bubbleDataBase<distribution<complexD>, bubbleDataPolicies> & me = upcast();
    int t;
    for(int t_expect=0;t_expect<me.getLt();t_expect++){
      if(!(in >> t)) error_exit(std::cout << "bubbleData::parse failed to read t for config " << sample << "\n");
      if(t != t_expect) error_exit(std::cout << "bubbleData::parse t doesn't match expectations: " << t << ":" << t_expect << " for config " << sample << "\n");

      double &re = reinterpret_cast<double(&)[2]>( me(t).sample(sample) )[0];
      double &im = reinterpret_cast<double(&)[2]>( me(t).sample(sample) )[1];
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

typedef bubbleDataBase<distribution<complexD>, bubbleDataPolicies> bubbleData;
typedef bubbleDataBase<doubleJackknifeDistribution<complexD> > bubbleDoubleJackData;

template<typename DistributionType, typename Policies = null_type>
class figureDataBase: public Policies{
  NumericMatrix<DistributionType> d; //(tsrc,tsep).sample(cfg)
  int Lt;
  template<typename T,typename P>
  friend std::ostream & operator<<(std::ostream &os, const figureDataBase<T,P> &f);
public:
  
  typedef figureDataBase<DistributionType,Policies> ET_tag;
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,figureDataBase<DistributionType,Policies> >::value, int>::type = 0>
  figureDataBase<DistributionType,Policies>(U&& expr): d(expr.common_properties()), Lt(expr.common_properties()){
    for(int i=0;i<Lt*Lt;i++)
      getElem<figureDataBase<DistributionType,Policies> >::elem(*this, i) = expr[i];    
  }

  figureDataBase() = default;
  figureDataBase(const figureDataBase<DistributionType,Policies> &r) = default;
  figureDataBase(figureDataBase<DistributionType,Policies> &&r) = default;
  figureDataBase(const int _Lt, const int _Nsample): Lt(_Lt), d(_Lt,DistributionType(_Nsample)){}
  
  figureDataBase<DistributionType,Policies> &operator=(const figureDataBase<DistributionType,Policies> &r) = default;
  figureDataBase<DistributionType,Policies> &operator=(figureDataBase<DistributionType,Policies> &&r) = default;

  void zero(){
    for(int i=0;i<Lt;i++)
      for(int j=0;j<Lt;j++)
	for(int s=0;s<d(i,j).size();s++)
	  d(i,j).sample(s) = complexD(0.);
  }
    
  
  void setup(const int _Lt, const int _Nsample){ Lt = _Lt; d.resize(_Lt, DistributionType(_Nsample)); }
  
  inline const int getLt() const{ return Lt; }
  inline const int getNsample() const{ return d(0,0).size(); }
  
  inline DistributionType & operator()(const int tsrc, const int tsep){ return d(tsrc,tsep); }
  inline const DistributionType & operator()(const int tsrc, const int tsep) const { return d(tsrc,tsep); }

  inline DistributionType & at(const int tsrc, const int tsep){ return d(tsrc,tsep); }
  inline const DistributionType & at(const int tsrc, const int tsep) const { return d(tsrc,tsep); }
};

template<typename DistributionType,typename Policies>
struct getElem<figureDataBase<DistributionType,Policies> >{
  static inline auto elem(figureDataBase<DistributionType,Policies> &v, const int i)->decltype(v(0,0)){ return v(i/v.getLt(), i%v.getLt()); }
  static inline auto elem(const figureDataBase<DistributionType,Policies> &v, const int i)->decltype(v(0,0)){ return v(i/v.getLt(), i%v.getLt()); }
  static inline int common_properties(const figureDataBase<DistributionType,Policies> &v){ return v.getLt(); }
};


template<typename DistributionType,typename Policies>
std::ostream & operator<<(std::ostream &os, const figureDataBase<DistributionType,Policies> &f){
  basicPrint<> p(os);
  p << f.d;
  return os;
}

class figureDataPolicies{
  inline figureDataBase<distribution<complexD>, figureDataPolicies> & upcast(){ return *static_cast< figureDataBase<distribution<complexD>, figureDataPolicies>* >(this); }
  inline const figureDataBase<distribution<complexD>, figureDataPolicies> & upcast() const{ return *static_cast< figureDataBase<distribution<complexD>, figureDataPolicies> const* >(this); }

public:
  void parseCDR(std::istream &in, const int sample){
    figureDataBase<distribution<complexD>, figureDataPolicies> &me = upcast();
    
    const int Lt = me.getLt();
    const int nelems = Lt*Lt;

    int tsrc,tsep;
    for(int e=0;e<nelems;e++){
      int tsep_expect = e % Lt;
      int tsrc_expect = e / Lt;

      if(!(in >> tsrc >> tsep)) error_exit(std::cout << "FigureData::parseCDR failed to read tsrc, tsep for config " << sample << "\n");
      if(tsep != tsep_expect || tsrc != tsrc_expect) error_exit(std::cout << "FigureData tsrc tsep don't match expectations: "
								<< tsrc << ":" << tsrc_expect << " " << tsep << ":" << tsep_expect
								<< " for config " << sample << "\n");
      double &re = reinterpret_cast<double(&)[2]>( me.at(tsrc,tsep).sample(sample) )[0];
      double &im = reinterpret_cast<double(&)[2]>( me.at(tsrc,tsep).sample(sample) )[1];
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

  //Data not measured on every tsrc usually
  bool isZero(const int tsrc) const{
    const figureDataBase<distribution<complexD>, figureDataPolicies> &me = upcast();
    
    for(int tsep=0;tsep<me.getLt();tsep++)
      for(int sample=0;sample<me.getNsample();sample++){
	const complexD & v = me.at(tsrc,tsep).sample(sample);
	if(v.real() != 0.0 || v.imag() != 0.0) return false;
      }
    return true;
  }
};

typedef figureDataBase<distribution<complexD> , figureDataPolicies> figureData;
typedef figureDataBase<doubleJackknifeDistribution<complexD> > figureDoubleJackData;







#endif
