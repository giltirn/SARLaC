#pragma once

#include "data_container_common.h"

CPSFIT_START_NAMESPACE

//A container for figure data (pipi, sigma, pipi->sigma)
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
typedef figureDataBase<bootJackknifeDistributionD, figureDataDistributionPolicies<bootJackknifeDistributionD> > figureDataBootJack;

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
template<>
struct _figureDataTypeSelector<bootJackknifeDistributionD>{ typedef figureDataBootJack type; };


template<typename DistributionType>
using figureDataSelect = typename _figureDataTypeSelector<DistributionType>::type;

CPSFIT_END_NAMESPACE
