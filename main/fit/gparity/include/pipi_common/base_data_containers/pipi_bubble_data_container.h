#pragma once

#include "data_container_common.h"

CPSFIT_START_NAMESPACE

#define DAIQIAN_COMPATIBILITY_MODE

//A container for pipi bubble data
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
typedef bubbleDataBase<bootJackknifeDistributionD > bubbleDataBootJack;

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
template<>
struct _bubbleDataTypeSelector<bootJackknifeDistributionD>{ typedef bubbleDataBootJack type; };



template<typename DistributionType>
using bubbleDataSelect = typename _bubbleDataTypeSelector<DistributionType>::type;



inline bubbleData reIm(const bubbleDataZ &in, const int reim){
  bubbleData out(in.getSrcSnk(), in.getLt(), in.getTsepPiPi());
  for(int t=0;t<in.getLt();t++)
    out(t) = rawDataDistributionD(in(t).size(), [&](const int s){ return reim == 0 ? in(t).sample(s).real() : in(t).sample(s).imag(); });
  return out;
}

CPSFIT_END_NAMESPACE
