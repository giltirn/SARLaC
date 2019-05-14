#include<distribution.h>
#include<distribution/block_double_jackknife.h>
#include<random.h>

using namespace CPSfit;

template<typename T>
bool equals(const T &a, const T &b, const double tol){
  if(a.size() != b.size()) return false;
  for(int i=0;i<a.size();i++) if(!equals(a.sample(i),b.sample(i),tol)) return false;
  return true;
}
template<>
bool equals<double>(const double &a, const double &b, const double tol){
  return fabs(a - b) <= tol;
}


int main(void){
  RNG.initialize(1234);
  
  { //Test that a bdj with bin size 1 and a dj are equivalent
    int bin_size = 1;
    int nsample = 100;

    rawDataDistribution<double> raw(nsample);
    for(int i=0;i<nsample;i++) raw.sample(i) = gaussianRandom<double>(1.0, 0.5);

    blockDoubleJackknifeDistribution<double> bdj(raw, bin_size);
    doubleJackknifeDistribution<double> dj(raw);

    assert(bdj == bdj); //test equality!
    
    assert(bdj.size() == dj.size());
    assert(bdj.sample(0).size() == dj.sample(0).size());
    assert(bdj.nSamplesUnbinned() == nsample);
    assert(bdj.binSize() == bin_size);
    
    jackknifeDistribution<double> bdj_mn = bdj.mean();
    jackknifeDistribution<double> dj_mn = dj.mean();
    std::cout << bdj_mn << " " << dj_mn << std::endl;
    assert(equals(bdj_mn, dj_mn, 1e-10));


    jackknifeDistribution<double> bdj_cov = blockDoubleJackknifeDistribution<double>::covariance(bdj, bdj);
    jackknifeDistribution<double> dj_cov = doubleJackknifeDistribution<double>::covariance(dj, dj);
    std::cout << bdj_cov << " " << dj_cov << std::endl;
    assert(equals(bdj_cov, dj_cov, 1e-10));
  }
  {
    //Test bin cropping
    int nsample = 100;
    int bin_size = 3;
    
    rawDataDistribution<double> raw(nsample);
    for(int i=0;i<nsample;i++) raw.sample(i) = gaussianRandom<double>(1.0, 0.5);

    int nexpect = 99;
    int oexpect = 33;
    int iexpect = 96;  // 99 - 3

    blockDoubleJackknifeDistribution<double> bdj(raw, bin_size);

    std::cout << bdj.size() << " " << bdj.sample(0).size() << std::endl;

    assert(bdj.size() == oexpect);
    assert(bdj.nSamplesUnbinned() == nexpect);
    assert(bdj.sample(0).size() == iexpect);
  }
  { //Test ET

    int bin_size = 4;
    int nsample = 100;

    rawDataDistribution<double> a(nsample);
    for(int i=0;i<nsample;i++) a.sample(i) = gaussianRandom<double>(1.0, 0.5);

    rawDataDistribution<double> b(nsample);
    for(int i=0;i<nsample;i++) b.sample(i) = gaussianRandom<double>(3.0, 0.5);
    
    blockDoubleJackknifeDistribution<double> bdj_a(a, bin_size);
    blockDoubleJackknifeDistribution<double> bdj_b(b, bin_size);

    assert(bdj_a != bdj_b); //test equality
    
    blockDoubleJackknifeDistribution<double> expect(nsample, bin_size);
    assert(expect.size() == nsample/bin_size);
    assert(expect.sample(0).size() == nsample - bin_size);

    for(int i=0;i<expect.size();i++)
      for(int j=0;j<expect.sample(i).size();j++)
	expect.sample(i).sample(j) = bdj_a.sample(i).sample(j) * bdj_a.sample(i).sample(j) + bdj_b.sample(i).sample(j);
    
    blockDoubleJackknifeDistribution<double> calc = bdj_a * bdj_a + bdj_b;
  
    assert(equals(calc, expect, 1e-10));
  }
  { //Test HDF5 IO
    
    int bin_size = 1;
    int nsample = 100;

    rawDataDistribution<double> raw(nsample);
    for(int i=0;i<nsample;i++) raw.sample(i) = gaussianRandom<double>(1.0, 0.5);

    blockDoubleJackknifeDistribution<double> bdj(raw, bin_size);
    
    {
      HDF5writer wr("test.hdf5");
      write(wr, bdj, "test_bdj");
    }

    blockDoubleJackknifeDistribution<double> bdj_rd;
    {
      HDF5reader rd("test.hdf5");
      read(rd, bdj_rd, "test_bdj");
    }
    assert(bdj_rd == bdj);

    std::vector<blockDoubleJackknifeDistribution<double> > v_bdj(2);
    for(int i=0;i<2;i++){
      rawDataDistribution<double> raw(nsample);
      gaussianRandom(raw, i+0.3, 0.2);
      
      v_bdj[i].resample(raw, bin_size);
    }
    std::vector<blockDoubleJackknifeDistribution<double> > v_bdj_rd;
    {
      HDF5writer wr("test.hdf5");
      write(wr, v_bdj, "test_bdj"); //don't flatten
    }  
    {
      HDF5reader rd("test.hdf5");
      read(rd, v_bdj_rd, "test_bdj");
    }
    assert(v_bdj_rd.size() == 2);
    for(int i=0;i<2;i++) assert(v_bdj_rd[i] == v_bdj[i]);
  }    
    

  std::cout << "Done" << std::endl;
  return 0;
}
