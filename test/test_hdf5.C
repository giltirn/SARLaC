#include<distribution.h>
#include<random.h>

using namespace CPSfit;

#ifdef HAVE_HDF5
template<typename DistributionType>
struct createOne{
  static inline DistributionType doit(){ return DistributionType(50); }
};
template<typename T>
struct createOne<superJackknifeDistribution<T> >{
  static inline superJackknifeDistribution<T> doit(){
    superJackknifeLayout *layout = new superJackknifeLayout;
    layout->addEnsemble("A",20);
    layout->addEnsemble("B",20);
    return superJackknifeDistribution<T>(*layout,0.);
  }
};


template<typename DistributionType>
struct setRandom{
  static inline void doit(DistributionType &d){ gaussianRandom(d,1.,1.); }
};
template<typename T, template<typename> class V>
struct setRandom<jackknifeCdistribution<T,V> >{
  static inline void doit(jackknifeCdistribution<T,V> &d){ gaussianRandom(d,1.,1.); d.best() = d.mean(); }
};
template<typename T, template<typename> class V>
struct setRandom<doubleJackknifeDistribution<T,V> >{
  static inline void doit(doubleJackknifeDistribution<T,V> &d){
    for(int i=0;i<d.size();i++) gaussianRandom(d.sample(i), 1., 1.);
  }
};
template<typename T>
struct setRandom<superJackknifeDistribution<T> >{
  static inline void doit(superJackknifeDistribution<T> &d){
    for(int i=-1;i<d.size();i++) gaussianRandom<T>(d.osample(i), 1., 1.);
  }
};



template<typename distributionType>
void testBasic(){
  std::cout << "Testing IO for " << printType<distributionType>() << std::endl;
  distributionType vw = createOne<distributionType>::doit();
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

  {
    std::cout << "Testing default-format vector IO for " << printType<distributionType>() << std::endl;
    std::vector<distributionType> vvw(2);
    vvw[0] = vvw[1] = vw;
    {
      HDF5writer writer("test_vdef.hdf5");
      write(writer,vvw,"vv");
    }
    std::vector<distributionType> vvr;
    {
      HDF5reader reader("test_vdef.hdf5");
      read(reader,vvr,"vv");
    }
    assert(vvr.size() == vvw.size());
    
    for(int i=0;i<vvw.size();i++)
      assert(vvr[i] == vvw[i]);
  }

  {
    std::cout << "Testing default-format vector<vector> IO for " << printType<distributionType>() << std::endl;
    std::vector<std::vector<distributionType> > vvw(2, std::vector<distributionType>(2));
    vvw[0][0] = vvw[0][1] = vvw[1][0] = vvw[1][1] = vw;
    {
      HDF5writer writer("test_vvdef.hdf5");
      write(writer,vvw,"vv");
    }
    std::vector<std::vector<distributionType> > vvr;
    {
      HDF5reader reader("test_vvdef.hdf5");
      read(reader,vvr,"vv");
    }
    assert(vvw.size() == vvr.size());
    for(int i=0;i<vvw.size();i++){
      assert(vvw[i].size() == vvr[i].size());
      for(int j=0;j<vvw[i].size();j++)      
	assert(vvr[i][j] == vvw[i][j]);
    }
  }
  

  {
    std::cout << "Testing conventional-format IO for " << printType<distributionType>() << std::endl;
    writeParamsStandard(vw, "test_std.hdf5");
    readParamsStandard(vr, "test_std.hdf5");
    assert(vr == vw);
  }

  {
    std::cout << "Testing conventional-format vector IO for " << printType<distributionType>() << std::endl;
    std::vector<distributionType> vvw(2);
    vvw[0] = vvw[1] = vw;
    writeParamsStandard(vvw, "test_vstd.hdf5");
  
    std::vector<distributionType> vvr(2);
    readParamsStandard(vvr, "test_vstd.hdf5");

    assert(vvr.size() == vvw.size());
    for(int i=0;i<vvw.size();i++)
      assert(vvr[i] == vvw[i]);
  }    
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

template<template<typename,template<typename> class> class DistributionType, template<typename> class V>
struct setRandom<DistributionType<S,V> >{
  static inline void doit(DistributionType<S,V> &d){
    DistributionType<double,basic_vector> a(d.size()); setRandom<DistributionType<double,basic_vector> >::doit(a);
    DistributionType<double,basic_vector> b(d.size()); setRandom<DistributionType<double,basic_vector> >::doit(b);
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
  testBasic<superJackknifeDistribution<double> >();
  {//test std::map
    std::cout << "Testing IO for std::map<std::string, int>\n";
    std::map<std::string, int> mp;
    mp["hello"] = 0;
    mp["world"] = 1;

    {
      HDF5writer writer("test.hdf5");
      write(writer,mp,"map");
    }
      
    std::map<std::string, int> mpr;
    {
      HDF5reader reader("test.hdf5");
      read(reader,mpr,"map");
    }
    assert(mp.size() == mpr.size());

    for(auto it = mp.begin(); it != mp.end(); it++){
      auto rit = mpr.find(it->first);
      assert(rit != mpr.end());
      assert(rit->second == it->second);
    }
    std::cout << "Test passed\n";
  }

#else
  std::cout << "HDF5 not being used\n";
#endif  
  return 0;
}
