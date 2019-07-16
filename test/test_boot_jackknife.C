#include<distribution.h>
#include<distribution/boot_jackknife.h>
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
  
  {
    //Test that with a resample table that does not change the order, that the samples produced by covariance are all equal to the jackknife covariance on the original sample
    int nboot = 10;
    int nsample = 100;
    rawDataDistribution<double> a(nsample), b(nsample);
    gaussianRandom(a, 1.0,0.5);
    gaussianRandom(b, 3.0,0.5);
    
    jackknifeDistribution<double> a_j(a);
    jackknifeDistribution<double> b_j(b);
    
    double cov_j = jackknifeDistribution<double>::covariance(a_j, b_j);

    std::vector<std::vector<int> > table(nboot, std::vector<int>(nsample));
    for(int i=0;i<nboot;i++)
      for(int j=0;j<nsample;j++) table[i][j] = j;
    
    bootJackknifeDistribution<double> a_bj(a, table, 68);
    bootJackknifeDistribution<double> b_bj(b, table, 68);

    bootstrapDistribution<double> cov_bj = bootJackknifeDistribution<double>::covariance(a_bj, b_bj);

    assert(cov_bj.size() == nboot);
    assert(equals(cov_bj.best(), cov_j, 1e-12));
    for(int b=0;b<nboot;b++) assert(equals(cov_bj.sample(b), cov_j, 1e-12));
  }

  {
    //Test that when the number of samples in the table is smaller than the original (for example if using a block table and the block size does not evenly divide nsample)
    //that it it equivalent to jackknife with those configurations removed

    int nboot = 10;
    int nsample_orig = 100;
    rawDataDistribution<double> a(nsample_orig), b(nsample_orig);
    gaussianRandom(a, 1.0,0.5);
    gaussianRandom(b, 3.0,0.5);
    
    int nsample_trunc = 90;
    rawDataDistribution<double> a_trunc(nsample_trunc, [&](const int s){ return a.sample(s); });
    rawDataDistribution<double> b_trunc(nsample_trunc, [&](const int s){ return b.sample(s); });

    jackknifeDistribution<double> a_j(a_trunc);
    jackknifeDistribution<double> b_j(b_trunc);
    
    double cov_j = jackknifeDistribution<double>::covariance(a_j, b_j);

    std::vector<std::vector<int> > table(nboot, std::vector<int>(nsample_trunc));
    for(int i=0;i<nboot;i++)
      for(int j=0;j<nsample_trunc;j++) table[i][j] = j;
    
    bootJackknifeDistribution<double> a_bj(a, table, 68);
    bootJackknifeDistribution<double> b_bj(b, table, 68);

    for(int b=0;b<nboot;b++){
      assert(a_bj.sample(b).size() == nsample_trunc);
      assert(b_bj.sample(b).size() == nsample_trunc);
    }

    bootstrapDistribution<double> cov_bj = bootJackknifeDistribution<double>::covariance(a_bj, b_bj);

    assert(cov_bj.size() == nboot);
    assert(equals(cov_bj.best(), cov_j, 1e-12));
    for(int b=0;b<nboot;b++)
      assert(equals(cov_bj.sample(b), cov_j, 1e-12));
  
  }



  { //Test ET

    int nboot = 10;
    int nsample = 100;

    rawDataDistribution<double> a(nsample), b(nsample);
    gaussianRandom(a, 1.0,0.5);
    gaussianRandom(b, 3.0,0.5);

    std::vector<std::vector<int> > table = resampleTable(RNG, nsample, nboot);

    bootJackknifeDistribution<double> a_bj(a, table, 68);
    bootJackknifeDistribution<double> b_bj(b, table, 68);

    assert(a_bj != b_bj); //test equality
    
    bootJackknifeDistribution<double> expect(bootJackknifeInitType(nsample, nboot, 68));
    assert(expect.size() == nboot);
    assert(expect.sample(0).size() == nsample);
    
    for(int i=0;i<expect.size();i++)
      for(int j=0;j<expect.sample(i).size();j++)
	expect.sample(i).sample(j) = a_bj.sample(i).sample(j) * a_bj.sample(i).sample(j) + b_bj.sample(i).sample(j);
    
    bootJackknifeDistribution<double> calc = a_bj*a_bj + b_bj;
  
    assert(equals(calc, expect, 1e-10));
  }
  { //Test HDF5 IO
    int nboot = 10;
    int nsample = 100;

    rawDataDistribution<double> raw(nsample);
    gaussianRandom(raw, 1.0,0.5);

    std::vector<std::vector<int> > table = resampleTable(RNG, nsample, nboot);

    bootJackknifeDistribution<double> bj(raw, table, 68);
    
    {
      HDF5writer wr("test.hdf5");
      write(wr, bj, "test_bj");
    }

    bootJackknifeDistribution<double> bj_rd;
    {
      HDF5reader rd("test.hdf5");
      read(rd, bj_rd, "test_bj");
    }
    assert(bj_rd == bj);

    std::vector<bootJackknifeDistribution<double> > v_bj(2);
    for(int i=0;i<2;i++){
      rawDataDistribution<double> raw(nsample);
      gaussianRandom(raw, i+0.3, 0.2);
      
      v_bj[i].resample(raw, table);
    }
    std::vector<bootJackknifeDistribution<double> > v_bj_rd;
    {
      HDF5writer wr("test.hdf5");
      write(wr, v_bj, "test_bj"); //don't flatten
    }  
    {
      HDF5reader rd("test.hdf5");
      read(rd, v_bj_rd, "test_bj");
    }
    assert(v_bj_rd.size() == 2);
    for(int i=0;i<2;i++) assert(v_bj_rd[i] == v_bj[i]);
  }    

  std::cout << "Done" << std::endl;
  return 0;
}
