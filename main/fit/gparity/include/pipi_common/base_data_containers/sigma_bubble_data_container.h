#pragma once

#include "data_container_common.h"

SARLAC_START_NAMESPACE

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
typedef sigmaSelfContractionBase<bootJackknifeDistributionD> sigmaSelfContractionBootJack;

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
template<>
struct _sigmaSelfContractionTypeSelector<bootJackknifeDistributionD>{ typedef sigmaSelfContractionBootJack type; };


template<typename DistributionType>
using sigmaSelfContractionSelect = typename _sigmaSelfContractionTypeSelector<DistributionType>::type;




inline sigmaSelfContraction reIm(const sigmaSelfContractionZ &in, const int reim){
  sigmaSelfContraction out(in.getLt());
  for(int t=0;t<in.getLt();t++)
    out(t) = rawDataDistributionD(in(t).size(), [&](const int s){ return reim == 0 ? in(t).sample(s).real() : in(t).sample(s).imag(); });
  return out;
}

SARLAC_END_NAMESPACE
