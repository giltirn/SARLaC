#include<distribution.h>
#include<random.h>
#include<containers/constrained_memory_vector.h>
#include<tensors/numeric_tensor.h>

using namespace SARLaC;

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
template<typename T>
struct createOne<superMultiDistribution<T> >{
  static inline superMultiDistribution<T> doit(){
    superMultiLayout *layout = new superMultiLayout;
    layout->addEnsemble("A",MultiType::Jackknife,20);
    layout->addEnsemble("B",MultiType::Bootstrap,20);
    return superMultiDistribution<T>(*layout,0.);
  }
};
template<typename T>
struct createOne<bootstrapDistribution<T> >{
  static inline bootstrapDistribution<T> doit(){ return bootstrapDistribution<T>(bootstrapInitType(50, 23)); } //use different numbers than default
};
template<typename T>
struct createOne<bootJackknifeDistribution<T> >{
  static inline bootJackknifeDistribution<T> doit(){ return bootJackknifeDistribution<T>(bootJackknifeInitType(50,49, 30, 23)); }
};


template<typename DistributionType>
struct setRandom{
  static inline void doit(DistributionType &d){ 
    typedef iterate<DistributionType> iter;
    for(int s=0;s<iter::size(d);s++) gaussianRandom(iter::at(s,d), 1.0, 1.0);
  }
};
template<typename T>
struct setRandom< distribution<T> >{
  static inline void doit(distribution<T> &d){ 
    for(int s=0;s<d.size();s++) gaussianRandom(d.sample(s), 1.0, 1.0);
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

  S(){}
  S(const double a, const double b): a(a), b(b){}

  inline int size() const{ return 2; }
  inline double & operator()(const int i){ return i==0 ? a : b; }
  inline const double & operator()(const int i) const{ return i==0 ? a : b; }

  bool operator==(const S &r) const{ return a == r.a && b == r.b; }
  bool operator!=(const S &r) const{ return !(*this == r); }
};

template<template<typename,template<typename> class> class DistributionType, template<typename> class V>
struct setRandom<DistributionType<S,V> >{
  static inline void doit(DistributionType<S,V> &d){
    typedef iterate<DistributionType<double,basic_vector> > iter_from;
    typedef iterate<DistributionType<S,V> > iter_to;

    DistributionType<double,basic_vector> a(d.getInitializer()); setRandom<DistributionType<double,basic_vector> >::doit(a);
    DistributionType<double,basic_vector> b(d.getInitializer()); setRandom<DistributionType<double,basic_vector> >::doit(b);
    for(int s=0;s<iter_to::size(d);s++){
      iter_to::at(s, d).a = iter_from::at(s, a);
      iter_to::at(s, d).b = iter_from::at(s, b);
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

  {
    std::cout << "Testing tensor" << std::endl;
    NumericTensor<rawDataDistribution<double>, 2> wtensor({10,10});
    NumericTensor<rawDataDistribution<double>, 2> rtensor;
    for(int i=0;i<10;i++){
      for(int j=0;j<10;j++){
	wtensor({i,j}).resize(100);
	gaussianRandom(wtensor({i,j}), 5.0, 1.0);
      }
    }

    {
      HDF5writer writer("test.hdf5");
      write(writer,wtensor,"tensor");
    }
  
    {
      HDF5reader reader("test.hdf5");
      read(reader,rtensor,"tensor");
    }
  
    assert(rtensor.size(0) == 10 && rtensor.size(1)==10);
    for(int i=0;i<10;i++)
      for(int j=0;j<10;j++)
	assert(rtensor({i,j}) == wtensor({i,j}));
  }

  {
    std::cout << "Testing pointer read/write" << std::endl;
    double v = 1234;
    std::string s = "hello";

    double const* vp = &v;
    std::string const *sp = &s;
    double const* np = NULL;

    std::vector<double*> a = { &v, NULL, &v };

    {
      HDF5writer writer("test.hdf5");
      writePointer(writer,vp,"vp");
      writePointer(writer,sp,"sp");
      writePointer(writer,np,"np");
      writePointer(writer,a,"a");
    }
  
    double* vpr, *npr;
    std::string* spr;
    std::vector<double*> ar;
    {
      HDF5reader reader("test.hdf5");
      readPointer(reader,vpr,"vp");
      readPointer(reader,spr,"sp");
      readPointer(reader,npr,"np");
      readPointer(reader,ar, "a");
    }  
    assert(*vpr == v);
    assert(*spr == s);
    assert(npr == NULL);
    assert(ar.size() == 3);
    assert(*ar[0] == v);
    assert(*ar[2] == v);
    assert(ar[1] == NULL);
  }



  testBasic<distribution<double> >();
  testBasic<rawDataDistribution<double> >();
  testBasic<jackknifeDistribution<double> >();
  testBasic<jackknifeCdistribution<double> >();
  testBasic<doubleJackknifeDistribution<double> >();
  testBasic<bootstrapDistribution<double> >();
  testBasic<bootJackknifeDistribution<double> >();

  testBasic<jackknifeDistribution<S> >();
  testBasic<superJackknifeDistribution<double> >();
  testBasic<superMultiDistribution<double> >();
  testBasic<bootstrapDistribution<S> >();

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
  {//test std::complex
    std::cout << "Testing IO for std::complex<double>" << std::endl;
    std::complex<double> v(3.14,-9.22);
    {
      HDF5writer writer("test.hdf5");
      write(writer,v,"value");
    }
    std::complex<double> u;
    {
      HDF5reader reader("test.hdf5");
      read(reader,u,"value");
    }
    assert(v == u);
  }

  {
    std::cout << "Testing complex, vector and distribution of complex" << std::endl;
    std::complex<double> v(2.2,3.3);
    std::vector<std::complex<double> > vv(3, v);
    jackknifeDistribution<std::complex<double> > vvv(3, v);

    std::complex<double> r;
    std::vector<std::complex<double> > rr;
    jackknifeDistribution<std::complex<double> > rrr;

    {
      HDF5writer wr("test.hdf5");
      write(wr, v, "z1");
      write(wr, vv, "zv1");    
      write(wr, vvv, "zd1");
    }
    {
      HDF5reader rd("test.hdf5");    
      read(rd, r, "z1");
      read(rd, rr, "zv1");
      read(rd, rrr, "zd1");
    }
  
    assert(r == v);
    assert(rrr == vvv);
    for(int i=0;i<3;i++)
      assert(rr[i] == vv[i]);

    //Test to ensure complex vector writes un uncompact format can also be read by default read command (ensuring backwards compatibility)
    {
      HDF5writer wr("test.hdf5");
      writeUncompact(wr, vv, "zv1");    
    }
    {
      HDF5reader rd("test.hdf5");    
      read(rd, rr, "zv1");
    }
  
    for(int i=0;i<3;i++)
      assert(rr[i] == vv[i]);
  }

  {
    //Test vector of string
    std::cout << "Test vector of string" << std::endl;
    std::vector<std::string> v = { "hello" , "world" };
    std::vector<std::string> r;
    {
      HDF5writer wr("test.hdf5");
      write(wr, v, "v");
    }
    {
      HDF5reader rd("test.hdf5");
      read(rd, r, "v");
    }
    assert( r == v );
  }

  std::cout << "Passed all tests" << std::endl;
#else
  std::cout << "HDF5 not being used\n";
#endif  
  return 0;
}
