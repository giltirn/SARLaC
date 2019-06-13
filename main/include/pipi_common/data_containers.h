#ifndef PIPI_DATA_CONTAINERS_H
#define PIPI_DATA_CONTAINERS_H

#define DAIQIAN_COMPATIBILITY_MODE

#include<map>

#include<config.h>
#include<utils/macros.h>

#include<distribution.h>
#include<tensors.h>
#include<common.h>

CPSFIT_START_NAMESPACE

enum SourceOrSink { Source, Sink };
inline void write(CPSfit::HDF5writer &writer, const SourceOrSink d, const std::string &tag){ write(writer,(int)d,tag); }
inline void read(CPSfit::HDF5reader &reader, SourceOrSink &d, const std::string &tag){
  int dd; read(reader,dd,tag); 
  d = (SourceOrSink)dd;
}

template<typename _DistributionType, typename Policies = empty_t>
class bubbleDataBase: public Policies{
public:
  typedef _DistributionType DistributionType;

private:
  NumericVector<DistributionType> d; //(t).sample(cfg)
  int Lt;
  int tsep_pipi;
  SourceOrSink src_snk;

public:
  bubbleDataBase() = default;
  bubbleDataBase(const bubbleDataBase &r) = default;
  bubbleDataBase(bubbleDataBase &&r) = default;
  bubbleDataBase(const SourceOrSink _src_snk, const int _Lt, const int _tsep_pipi, const DistributionType &init = DistributionType()): src_snk(_src_snk), Lt(_Lt), tsep_pipi(_tsep_pipi), d(_Lt,init){}

  bubbleDataBase& operator=(const bubbleDataBase &r) = default;
  bubbleDataBase& operator=(bubbleDataBase &&r) = default;

  void setup(const SourceOrSink _src_snk, const int _Lt, const int _tsep_pipi, const DistributionType &init = DistributionType()){ 
    src_snk = _src_snk; 
    tsep_pipi = _tsep_pipi; 
    Lt = _Lt; 
    d.resize(Lt, init); 
  }
  
  inline int getLt() const{ return Lt; }
  inline int getTsepPiPi() const{ return tsep_pipi; }
  inline SourceOrSink getSrcSnk() const{ return src_snk; }

  inline DistributionType & operator()(const int t){ return d[t]; } //time coordinate is for pi1 (earlier pion for sink, later for src)
  inline const DistributionType & operator()(const int t) const { return d[t]; }

  inline DistributionType & at(const int t){ return d[t]; }
  inline const DistributionType & at(const int t) const { return d[t]; }

  void zero(){
    for(int i=0;i<Lt;i++)
      zeroit(d[i]);
  }
  void initializeElements(const DistributionType &init){ 
    for(int i=0;i<Lt;i++)
      d[i] = init;
  }  

  GENERATE_HDF5_SERIALIZE_METHOD((Lt)(tsep_pipi)(d));
};
#ifdef HAVE_HDF5
template<typename D, typename P>
inline void write(HDF5writer &writer, const bubbleDataBase<D,P> &d, const std::string &tag){ d.write(writer,tag); }
template<typename D, typename P>
inline void read(HDF5reader &reader, bubbleDataBase<D,P> &d, const std::string &tag){ d.read(reader,tag); }
#endif


//Useful policies for setting output type to real or re/im
struct setValueRealPart{
  template<typename T>
  static inline void set(T &into, const double re, const double im){ into = re; }
};
struct setValueComplex{
  template<typename T>
  static inline void set(T &into, const double re, const double im){ into.real(re); into.imag(im); }
};

//Policy for parsing raw bubble data
template<typename ValueType, typename ValueTypeSetPolicy>
class bubbleDataPolicies{
  typedef bubbleDataPolicies<ValueType, ValueTypeSetPolicy> MyType;
  typedef rawDataDistribution<ValueType> DistType;
  typedef bubbleDataBase<DistType, MyType> MasterType;

  inline MasterType & upcast(){ return *static_cast<MasterType*>(this); }
  inline const MasterType & upcast() const{ return *static_cast<MasterType const* >(this); }

public:
  void parse(std::istream &in, const int sample){
    auto & me = upcast();
    int t;
    for(int t_expect=0;t_expect<me.getLt();t_expect++){
      if(!(in >> t)) error_exit(std::cout << "bubbleData::parse failed to read t for config " << sample << "\n");
      if(t != t_expect) error_exit(std::cout << "bubbleData::parse t doesn't match expectations: " << t << ":" << t_expect << " for config " << sample << "\n");

      if(me.getSrcSnk() == Sink) t = (t - me.getTsepPiPi() + me.getLt() ) % me.getLt(); //the raw data file time coordinate is that of the later pion

      double re, im;
      if(!(in >> re >> im)) error_exit(std::cout << "bubbleData::parse failed to real values for config " << sample << "\n");
#ifdef DAIQIAN_COMPATIBILITY_MODE
      re = -re; //correct for missing minus sign
      im = -im;
#endif      
      ValueTypeSetPolicy::set(me(t).sample(sample), re, im);
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
    auto & me = upcast();
    const int Lt = me.getLt();
#pragma omp parallel for
    for(int t=0;t<Lt;t++){
      me.at(t) = me.at(t).bin(bin_size);
    }
  }

};

typedef bubbleDataBase<rawDataDistributionD, bubbleDataPolicies<double, setValueRealPart> > bubbleData;
typedef bubbleDataBase<rawDataDistribution<std::complex<double> >, bubbleDataPolicies<std::complex<double>, setValueComplex> > bubbleDataZ;
typedef bubbleDataBase<jackknifeDistributionD > bubbleDataJack;
typedef bubbleDataBase<doubleJackknifeDistributionD > bubbleDataDoubleJack;
typedef bubbleDataBase<blockDoubleJackknifeDistributionD > bubbleDataBlockDoubleJack;
typedef bubbleDataBase<bootstrapDistributionD > bubbleDataBoot;

template<typename DistributionType>
struct _bubbleDataTypeSelector{};
template<>
struct _bubbleDataTypeSelector<rawDataDistributionD>{ typedef bubbleData type; };
template<>
struct _bubbleDataTypeSelector<rawDataDistribution<std::complex<double> > >{ typedef bubbleDataZ type; };
template<>
struct _bubbleDataTypeSelector<jackknifeDistributionD>{ typedef bubbleDataJack type; };
template<>
struct _bubbleDataTypeSelector<doubleJackknifeDistributionD>{ typedef bubbleDataDoubleJack type; };
template<>
struct _bubbleDataTypeSelector<blockDoubleJackknifeDistributionD>{ typedef bubbleDataBlockDoubleJack type; };
template<>
struct _bubbleDataTypeSelector<bootstrapDistributionD>{ typedef bubbleDataBoot type; };

template<typename DistributionType>
using bubbleDataSelect = typename _bubbleDataTypeSelector<DistributionType>::type;



inline bubbleData reIm(const bubbleDataZ &in, const int reim){
  bubbleData out(in.getSrcSnk(), in.getLt(), in.getTsepPiPi());
  for(int t=0;t<in.getLt();t++)
    out(t) = rawDataDistributionD(in(t).size(), [&](const int s){ return reim == 0 ? in(t).sample(s).real() : in(t).sample(s).imag(); });
  return out;
}

template<typename _DistributionType, typename Policies = empty_t>
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
  figureDataBase(const int _Lt, const DistributionType &init = DistributionType()): Lt(_Lt), d(_Lt,init){}
  
  figureDataBase<DistributionType,Policies> &operator=(const figureDataBase<DistributionType,Policies> &r) = default;
  figureDataBase<DistributionType,Policies> &operator=(figureDataBase<DistributionType,Policies> &&r) = default;

  void zero(){
    for(int i=0;i<Lt;i++)
      for(int j=0;j<Lt;j++)
	zeroit(d(i,j));
  }
  void initializeElements(const DistributionType &init){ 
    for(int i=0;i<Lt;i++)
      for(int j=0;j<Lt;j++)
	d(i,j) = init;
  }  
  
  void setup(const int _Lt, const DistributionType &init = DistributionType()){ Lt = _Lt; d.resize(_Lt,init); }
  
  inline int getLt() const{ return Lt; }
  
  inline DistributionType & operator()(const int tsrc, const int tsep){ return d(tsrc,tsep); }
  inline const DistributionType & operator()(const int tsrc, const int tsep) const { return d(tsrc,tsep); }

  inline DistributionType & at(const int tsrc, const int tsep){ return d(tsrc,tsep); }
  inline const DistributionType & at(const int tsrc, const int tsep) const { return d(tsrc,tsep); }

  GENERATE_HDF5_SERIALIZE_METHOD((Lt)(d));
};
#ifdef HAVE_HDF5
template<typename D, typename P>
inline void write(HDF5writer &writer, const figureDataBase<D,P> &d, const std::string &tag){ d.write(writer,tag); }
template<typename D, typename P>
inline void read(HDF5reader &reader, figureDataBase<D,P> &d, const std::string &tag){ d.read(reader,tag); }
#endif

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
  typedef figureDataBase<rawDataDistributionD, figureDataPolicies> MasterType;
  inline MasterType & upcast(){ return *static_cast< MasterType* >(this); }
  inline const MasterType & upcast() const{ return *static_cast< MasterType const* >(this); }

public:
  void parseCDR(std::istream &in, const int sample){
    auto &me = upcast();
    
    const int Lt = me.getLt();
    const int nelems = Lt*Lt;

    int tsrc,tsep;
    for(int e=0;e<nelems;e++){
      int tsep_expect = e % Lt;
      int tsrc_expect = e / Lt;

      if(!(in >> tsrc >> tsep)) error_exit(std::cout << "FigureData::parseCDR (thread " << omp_get_thread_num() << " failed to read tsrc, tsep for config " << sample 
					   << ". Expected tsrc " << tsrc_expect << " tsnk " << tsep_expect <<". Badbit=" << in.bad() << "\n");
      if(tsep != tsep_expect || tsrc != tsrc_expect) error_exit(std::cout << "FigureData tsrc tsep don't match expectations: "
								<< tsrc << ":" << tsrc_expect << " " << tsep << ":" << tsep_expect
								<< " for config " << sample << "\n");
      double &re = me.at(tsrc,tsep).sample(sample);
      double im; //discard because it averages to zero
      if(!(in >> re >> im)) error_exit(std::cout << "FigureData::parseCDR failed to read values for config " << sample 
				       << " tsrc " << tsrc_expect << " tsnk " << tsep_expect << ". Badbit=" << in.bad() << "\n");
    }
  }
  
  //Optional cache passed in of previously read data on this configuration
  typedef std::map< std::string, std::vector<double> > CacheType;

  void parseCDR(const std::string &filename, const int sample, CacheType* sample_cache = NULL){
    auto & me = upcast();

    CacheType::const_iterator it;
    if(sample_cache != NULL && ( (it = sample_cache->find(filename)) != sample_cache->end() )){
      int i = 0;
      for(int tsrc=0;tsrc<me.getLt();tsrc++) 
	for(int t=0;t<me.getLt();t++) 
	  me.at(tsrc, t).sample(sample) = it->second[i++];
      std::cout << "Retrieved " << filename << " from cache" << std::endl;
      return;
    }

    std::ifstream is(filename.c_str());
    if(!is.good()){ std::cout << "Could not open file \"" << filename << "\"\n"; std::cout.flush(); exit(-1); }
    parseCDR(is,sample);
    if(is.fail() || is.bad()){ std::cout << "Error reading file \"" << filename << "\"\n"; std::cout.flush(); exit(-1); }
    is.close();

    if(sample_cache != NULL){
      std::vector<double> v(me.getLt() * me.getLt());
      int i = 0;
      for(int tsrc=0;tsrc<me.getLt();tsrc++) 
	for(int t=0;t<me.getLt();t++) 
	  v[i++] = me.at(tsrc, t).sample(sample);
      (*sample_cache)[filename] = std::move(v);
    }
  }    

  //Data not measured on every tsrc usually
  bool isZero(const int tsrc) const{
    auto &me = upcast();
    
    for(int tsep=0;tsep<me.getLt();tsep++)
      for(int sample=0;sample<me.at(tsrc,tsep).size();sample++)
	if( me.at(tsrc,tsep).sample(sample) != 0.0 ) return false;
    return true;
  }

  void bin(const int bin_size){
    auto &me = upcast();
    const int Lt = me.getLt();
#pragma omp parallel for
    for(int i=0;i<Lt*Lt;i++){
      int tsep = i % Lt;
      int tsrc = i / Lt;
      me.at(tsrc,tsep) = me.at(tsrc,tsep).bin(bin_size);
    }
  }
};

template<typename DistributionType>
class figureDataDistributionPolicies{
  typedef figureDataBase<DistributionType, figureDataDistributionPolicies<DistributionType> > MasterType;

  inline MasterType & upcast(){ return *static_cast<MasterType* >(this); }
  inline const MasterType & upcast() const{ return *static_cast<MasterType const* >(this); }

public:
  bool isZero(const int tsrc) const{
    const auto &me = upcast();
    
    for(int tsep=0;tsep<me.getLt();tsep++){
      const auto &dist = me.at(tsrc,tsep);
      for(int i=0;i<iterate<DistributionType>::size(dist);i++)
	if( iterate<DistributionType>::at(i, dist) != 0. ) return false;
    }      
    return true;
  }
};

  

typedef figureDataBase<rawDataDistributionD , figureDataPolicies> figureData;
typedef figureDataBase<doubleJackknifeDistributionD, figureDataDistributionPolicies<doubleJackknifeDistributionD> > figureDataDoubleJack;
typedef figureDataBase<blockDoubleJackknifeDistributionD, figureDataDistributionPolicies<blockDoubleJackknifeDistributionD> > figureDataBlockDoubleJack;
typedef figureDataBase<jackknifeDistributionD, figureDataDistributionPolicies<jackknifeDistributionD> > figureDataJack;
typedef figureDataBase<bootstrapDistributionD, figureDataDistributionPolicies<bootstrapDistributionD> > figureDataBoot;

template<typename DistributionType>
struct _figureDataTypeSelector{};
template<>
struct _figureDataTypeSelector<rawDataDistributionD>{ typedef figureData type; };
template<>
struct _figureDataTypeSelector<jackknifeDistributionD>{ typedef figureDataJack type; };
template<>
struct _figureDataTypeSelector<doubleJackknifeDistributionD>{ typedef figureDataDoubleJack type; };
template<>
struct _figureDataTypeSelector<blockDoubleJackknifeDistributionD>{ typedef figureDataBlockDoubleJack type; };
template<>
struct _figureDataTypeSelector<bootstrapDistributionD>{ typedef figureDataBoot type; };

template<typename DistributionType>
using figureDataSelect = typename _figureDataTypeSelector<DistributionType>::type;


template<typename _DistributionType, typename Policies = empty_t>
class sigmaSelfContractionBase: public Policies{
public:
  typedef _DistributionType DistributionType;
private:
  NumericVector<DistributionType> d; //(t).sample(cfg)
  int Lt;

public:
  sigmaSelfContractionBase(){}
  sigmaSelfContractionBase(const int _Lt, const DistributionType &init = DistributionType()): Lt(_Lt), d(_Lt,init){}

  void setup(const int _Lt, const DistributionType &init = DistributionType()){ 
    Lt = _Lt; 
    d.resize(_Lt,init); 
  }
  
  inline int getLt() const{ return Lt; }

  inline DistributionType & operator()(const int t){ return d[t]; } //time coordinate is for pi1 (earlier pion for sink, later for src)
  inline const DistributionType & operator()(const int t) const { return d[t]; }

  inline DistributionType & at(const int t){ return d[t]; }
  inline const DistributionType & at(const int t) const { return d[t]; }

  void zero(){
    for(int t=0;t<Lt;t++) zeroit(d(t));
  }

  GENERATE_HDF5_SERIALIZE_METHOD((Lt)(d));
};
#ifdef HAVE_HDF5
template<typename D,typename P>
inline void write(HDF5writer &writer, const sigmaSelfContractionBase<D,P> &d, const std::string &tag){ d.write(writer,tag); }
template<typename D,typename P>
inline void read(HDF5reader &reader, sigmaSelfContractionBase<D,P> &d, const std::string &tag){ d.read(reader,tag); }
#endif

template<typename ValueType, typename ValueTypeSetPolicy>
class sigmaSelfContractionPolicies{
  typedef sigmaSelfContractionPolicies<ValueType, ValueTypeSetPolicy> MyType;
  typedef rawDataDistribution<ValueType> DistType;
  typedef sigmaSelfContractionBase<DistType, MyType> MasterType;

  inline MasterType & upcast(){ return *static_cast<MasterType*>(this); }
  inline const MasterType & upcast() const{ return *static_cast<MasterType const* >(this); }

public:

  void parse(std::istream &in, const int sample){
    auto &me = upcast();

    int t;
    for(int t_expect=0;t_expect<me.getLt();t_expect++){
      if(!(in >> t)) error_exit(std::cout << "sigmaSelfContraction::parse failed to read t for config " << sample << "\n");
      if(t != t_expect) error_exit(std::cout << "sigmaSelfContraction::parse t doesn't match expectations: " << t << ":" << t_expect << " for config " << sample << "\n");
      
      double re, im;
      if(!(in >> re >> im)) error_exit(std::cout << "sigmaSelfContraction::parse failed to real values for config " << sample << "\n");
      ValueTypeSetPolicy::set(me(t).sample(sample), re, im); 
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
    auto &me = upcast();

#pragma omp parallel for
    for(int t=0;t<me.getLt();t++){
      me(t) = me(t).bin(bin_size);
    }
  }
};


typedef sigmaSelfContractionBase<rawDataDistributionD , sigmaSelfContractionPolicies<double, setValueRealPart> > sigmaSelfContraction;
typedef sigmaSelfContractionBase<rawDataDistribution<std::complex<double> > , sigmaSelfContractionPolicies<std::complex<double>, setValueComplex> > sigmaSelfContractionZ;
typedef sigmaSelfContractionBase<jackknifeDistributionD> sigmaSelfContractionJack;
typedef sigmaSelfContractionBase<doubleJackknifeDistributionD> sigmaSelfContractionDoubleJack;
typedef sigmaSelfContractionBase<blockDoubleJackknifeDistributionD> sigmaSelfContractionBlockDoubleJack;
typedef sigmaSelfContractionBase<bootstrapDistributionD> sigmaSelfContractionBoot;

template<typename DistributionType>
struct _sigmaSelfContractionTypeSelector{};
template<>
struct _sigmaSelfContractionTypeSelector<rawDataDistributionD>{ typedef sigmaSelfContraction type; };
template<>
struct _sigmaSelfContractionTypeSelector<rawDataDistribution<std::complex<double> > >{ typedef sigmaSelfContractionZ type; };
template<>
struct _sigmaSelfContractionTypeSelector<jackknifeDistributionD>{ typedef sigmaSelfContractionJack type; };
template<>
struct _sigmaSelfContractionTypeSelector<doubleJackknifeDistributionD>{ typedef sigmaSelfContractionDoubleJack type; };
template<>
struct _sigmaSelfContractionTypeSelector<blockDoubleJackknifeDistributionD>{ typedef sigmaSelfContractionBlockDoubleJack type; };
template<>
struct _sigmaSelfContractionTypeSelector<bootstrapDistributionD>{ typedef sigmaSelfContractionBoot type; };

template<typename DistributionType>
using sigmaSelfContractionSelect = typename _sigmaSelfContractionTypeSelector<DistributionType>::type;




inline sigmaSelfContraction reIm(const sigmaSelfContractionZ &in, const int reim){
  sigmaSelfContraction out(in.getLt());
  for(int t=0;t<in.getLt();t++)
    out(t) = rawDataDistributionD(in(t).size(), [&](const int s){ return reim == 0 ? in(t).sample(s).real() : in(t).sample(s).imag(); });
  return out;
}




CPSFIT_END_NAMESPACE

#endif
