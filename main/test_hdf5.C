#include<config.h>
#include<distribution.h>
#include<random.h>

#ifdef HAVE_HDF5

template<typename distributionType>
void testBasic(){
  std::cout << "Testing IO for " << printType<distributionType>() << std::endl;
  distributionType vw(50);
  gaussianRandom(vw, 1, 1);

  {
    HDF5writer writer("test.hdf5");
    write(writer,vw,"v");
  }
      
  distributionType vr;
  {
    HDF5reader reader("test.hdf5");
    read(reader,vr,"v");

    assert(vr.size() == vw.size());
    for(int i=0;i<vw.size();i++){
      if(vr.sample(i)!=vw.sample(i)){
	std::cout << "Error " << i << " " << vr.sample(i) << " " << vw.sample(i) << std::endl;
	assert(0);
      }
    }
  }
  std::cout << "Test passed\n\n\n";
}

template<typename T>
void testJackknifeCdistribution(){
  std::cout << "Testing IO for " << printType<jackknifeCdistribution<T> >() << std::endl;
  jackknifeCdistribution<T> vw(50);
  gaussianRandom(vw, 1, 1);
  vw.best() = vw.mean();

  {
    HDF5writer writer("test.hdf5");
    write(writer,vw,"v");
  }
      
  jackknifeCdistribution<T> vr;
  {
    HDF5reader reader("test.hdf5");
    read(reader,vr,"v");

    assert(vr.best() == vw.best());
    assert(vr.size() == vw.size());
    for(int i=0;i<vw.size();i++){
      if(vr.sample(i)!=vw.sample(i)){
	std::cout << "Error " << i << " " << vr.sample(i) << " " << vw.sample(i) << std::endl;
	assert(0);
      }
    }
  }
  std::cout << "Test passed\n\n\n";
}

template<typename T>
void testDoubleJackknifeDistribution(){
  std::cout << "Testing IO for " << printType<doubleJackknifeDistribution<T> >() << std::endl;
  doubleJackknifeDistribution<T> vw(50);
  for(int i=0;i<50;i++)
    gaussianRandom(vw.sample(i), 1, 1);

  {
    HDF5writer writer("test.hdf5");
    write(writer,vw,"v");
  }
      
  doubleJackknifeDistribution<T> vr;
  {
    HDF5reader reader("test.hdf5");
    read(reader,vr,"v");

    assert(vr.size() == vw.size());
    for(int i=0;i<vw.size();i++){
      assert(vr.sample(i).size() == vw.sample(i).size());
      for(int j=0;j<vw.sample(i).size();j++){
	if(vr.sample(i).sample(j)!=vw.sample(i).sample(j)){
	  std::cout << "Error " << i << " " << j << " " << vr.sample(i).sample(j) << " " << vw.sample(i).sample(j) << std::endl;
	  assert(0);
	}
      }
    }
  }
  std::cout << "Test passed\n\n\n";
}




#endif




int main(void){
  RNG.initialize(1245);
#ifdef HAVE_HDF5
  std::cout << "Running test\n\n\n";
  testBasic<distribution<double> >();
  testBasic<rawDataDistribution<double> >();
  testBasic<jackknifeDistribution<double> >();
  testJackknifeCdistribution<double>();
  testDoubleJackknifeDistribution<double>();
#else
  std::cout << "HDF5 not being used\n";
#endif  
  return 0;
}
