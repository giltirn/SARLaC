#include<config.h>
#include<distribution.h>
#include<random.h>

#ifdef HAVE_HDF5

template<typename DistributionType>
struct setRandom{
  static inline void doit(DistributionType &d){ gaussianRandom(d,1,1); }
};
template<typename T>
struct setRandom<jackknifeCdistribution<T> >{
  static inline void doit(jackknifeCdistribution<T> &d){ gaussianRandom(d,1,1); d.best() = d.mean(); }
};
template<typename T>
struct setRandom<doubleJackknifeDistribution<T> >{
  static inline void doit(doubleJackknifeDistribution<T> &d){
    for(int i=0;i<d.size();i++) gaussianRandom(d.sample(i), 1, 1);
  }
};

template<typename distributionType>
void testBasic(){
  std::cout << "Testing IO for " << printType<distributionType>() << std::endl;
  distributionType vw(50);
  setRandom<distributionType>::doit(vw);

  {
    HDF5writer writer("test.hdf5");
    write(writer,vw,"v");
  }
      
  distributionType vr;
  {
    HDF5reader reader("test.hdf5");
    read(reader,vr,"v");

    assert(vr == vw);
  }
  std::cout << "Testing standard-format IO for " << printType<distributionType>() << std::endl;
  writeParamsStandard(vw, "test_std.hdf5");
  readParamsStandard(vr, "test_std.hdf5");
  assert(vr == vw);
  
  std::cout << "Testing standard-format vector IO for " << printType<distributionType>() << std::endl;
  std::vector<distributionType> vvw(2);
  vvw[0] = vvw[1] = vw;
  writeParamsStandard(vvw, "test_vstd.hdf5");
  
  std::vector<distributionType> vvr(2);
  readParamsStandard(vvr, "test_vstd.hdf5");

  assert(vvr.size() == vvw.size());
  for(int i=0;i<vvw.size();i++)
    assert(vvr[i] == vvw[i]);
    
  std::cout << "Test passed\n\n\n";
}


struct S{
  double a;
  double b;

  inline int size() const{ return 2; }
  inline double & operator()(const int i){ return i==0 ? a : b; }
  inline const double & operator()(const int i) const{ return i==0 ? a : b; }

  bool operator==(const S &r) const{ return a == r.a && b == r.b; }
  bool operator!=(const S &r) const{ return !(*this == r); }
};

template<template<typename> class DistributionType>
struct setRandom<DistributionType<S> >{
  static inline void doit(DistributionType<S> &d){
    DistributionType<double> a(d.size()); setRandom<DistributionType<double> >::doit(a);
    DistributionType<double> b(d.size()); setRandom<DistributionType<double> >::doit(b);
    for(int s=0;s<d.size();s++){
      d.sample(s).a = a.sample(s);
      d.sample(s).b = b.sample(s);
    }
  }
};

void write(HDF5writer &writer, const S &value, const std::string &tag){
  writer.enter(tag);
  write(writer,value.a, "a");
  write(writer,value.b, "b");
  writer.leave();
}
void read(HDF5reader &reader, S &value, const std::string &tag){
  reader.enter(tag);
  read(reader,value.a, "a");
  read(reader,value.b, "b");
  reader.leave();
}
#endif




int main(void){
  RNG.initialize(1245);
#ifdef HAVE_HDF5
  std::cout << "Running test\n\n\n";
  testBasic<distribution<double> >();
  testBasic<rawDataDistribution<double> >();
  testBasic<jackknifeDistribution<double> >();
  testBasic<jackknifeCdistribution<double> >();
  testBasic<doubleJackknifeDistribution<double> >();

  testBasic<jackknifeDistribution<S> >();
#else
  std::cout << "HDF5 not being used\n";
#endif  
  return 0;
}
