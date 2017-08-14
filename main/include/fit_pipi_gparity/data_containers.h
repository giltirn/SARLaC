#ifndef PIPI_DATA_CONTAINERS_H
#define PIPI_DATA_CONTAINERS_H

#define DAIQIAN_COMPATIBILITY_MODE
struct null_type{};

template<typename DistributionType, typename Policies = null_type>
class bubbleDataBase: public Policies{
  NumericVector<DistributionType> d; //(t).sample(cfg)
  int Lt;

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & d & Lt;
  }
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
  inline bubbleDataBase<rawDataDistributionD, bubbleDataPolicies> & upcast(){ return *static_cast< bubbleDataBase<rawDataDistributionD, bubbleDataPolicies>* >(this); }
  inline const bubbleDataBase<rawDataDistributionD, bubbleDataPolicies> & upcast() const{ return *static_cast< bubbleDataBase<rawDataDistributionD, bubbleDataPolicies> const* >(this); }

public:
  void parse(std::istream &in, const int sample){
    bubbleDataBase<rawDataDistributionD, bubbleDataPolicies> & me = upcast();
    int t;
    for(int t_expect=0;t_expect<me.getLt();t_expect++){
      if(!(in >> t)) error_exit(std::cout << "bubbleData::parse failed to read t for config " << sample << "\n");
      if(t != t_expect) error_exit(std::cout << "bubbleData::parse t doesn't match expectations: " << t << ":" << t_expect << " for config " << sample << "\n");

      double &re = me(t).sample(sample);
      double im; //discard because it is zero
      if(!(in >> re >> im)) error_exit(std::cout << "bubbleData::parse failed to real values for config " << sample << "\n");
#ifdef DAIQIAN_COMPATIBILITY_MODE
      re = -re; //correct for missing minus sign
#endif      
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

typedef bubbleDataBase<rawDataDistributionD, bubbleDataPolicies> bubbleData;
typedef bubbleDataBase<doubleJackknifeDistributionD > bubbleDataDoubleJack;

template<typename _DistributionType, typename Policies = null_type>
class figureDataBase: public Policies{
public:
  typedef _DistributionType DistributionType;
private:
  NumericSquareMatrix<DistributionType> d; //(tsrc,tsep).sample(cfg)
  int Lt;
  template<typename T,typename P>
  friend std::ostream & operator<<(std::ostream &os, const figureDataBase<T,P> &f);

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & d & Lt;
  }
public:
  typedef figureDataBase<DistributionType,Policies> ET_tag;
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,figureDataBase<DistributionType,Policies> >::value, int>::type = 0>
  figureDataBase<DistributionType,Policies>(U&& expr): d(expr.common_properties()), Lt(expr.common_properties()){
#pragma omp parallel for
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
	zeroit(d(i,j));
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
  for(int tsrc=0;tsrc<f.d.size();tsrc++)
    for(int tsep=0;tsep<f.d.size();tsep++)
      os << tsrc << " " << tsep << " " << f(tsrc,tsep) << std::endl;
  return os;
}

class figureDataPolicies{
  inline figureDataBase<rawDataDistributionD, figureDataPolicies> & upcast(){ return *static_cast< figureDataBase<rawDataDistributionD, figureDataPolicies>* >(this); }
  inline const figureDataBase<rawDataDistributionD, figureDataPolicies> & upcast() const{ return *static_cast< figureDataBase<rawDataDistributionD, figureDataPolicies> const* >(this); }

public:
  void parseCDR(std::istream &in, const int sample){
    figureDataBase<rawDataDistributionD, figureDataPolicies> &me = upcast();
    
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
      double &re = me.at(tsrc,tsep).sample(sample);
      double im; //discard because it averages to zero
      if(!(in >> re >> im)) error_exit(std::cout << "FigureData::parseCDR failed to read values for config " << sample << "\n");
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
    const figureDataBase<rawDataDistributionD, figureDataPolicies> &me = upcast();
    
    for(int tsep=0;tsep<me.getLt();tsep++)
      for(int sample=0;sample<me.getNsample();sample++)
	if( me.at(tsrc,tsep).sample(sample) != 0.0 ) return false;
    return true;
  }
};

class figureDataDoubleJackPolicies{
  inline figureDataBase<doubleJackknifeDistributionD, figureDataDoubleJackPolicies> & upcast(){ return *static_cast< figureDataBase<doubleJackknifeDistributionD, figureDataDoubleJackPolicies>* >(this); }
  inline const figureDataBase<doubleJackknifeDistributionD, figureDataDoubleJackPolicies> & upcast() const{ return *static_cast< figureDataBase<doubleJackknifeDistributionD, figureDataDoubleJackPolicies> const* >(this); }

public:
  bool isZero(const int tsrc) const{
    const figureDataBase<doubleJackknifeDistributionD, figureDataDoubleJackPolicies> &me = upcast();
    
    for(int tsep=0;tsep<me.getLt();tsep++)
      for(int sample=0;sample<me.getNsample();sample++)
	for(int sub_sample=0;sub_sample<me.getNsample()-1;sub_sample++)
	  if( me.at(tsrc,tsep).sample(sample).sample(sub_sample) != 0. ) return false;
      
    return true;
  }
};

  

typedef figureDataBase<rawDataDistributionD , figureDataPolicies> figureData;
typedef figureDataBase<doubleJackknifeDistributionD, figureDataDoubleJackPolicies > figureDataDoubleJack;







#endif
