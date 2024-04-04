#include<array>
#include<vector>
#include<complex>
#include<iostream>
#include<fstream>

#include<tensors.h>

using namespace SARLaC;

int main(void){

  typedef NumericRectMatrix<double> M;
  typedef NumericVector<double> V;
  typedef NumericSquareMatrix<double> SM;

  {
    M m(2,3,double(0.));
    std::cout << m << std::endl;

    for(int i=0;i<2;i++)
      for(int j=0;j<3;j++)
	assert(m(i,j) == 0.);
  }
  {
    M m(2,3,{0.1,0.2,0.3,1.1,1.2,1.3});
    std::cout << m << std::endl;
  
    for(int i=0;i<2;i++)
      for(int j=0;j<3;j++)
	assert( fabs(m(i,j) - (i+0.1*(j+1))) < 1e-8 );
  }
  {
    M m(2,3,[](const int i,const int j){ return j+3*i+3.14; });
    std::cout << m << std::endl;
  
    for(int i=0;i<2;i++)
      for(int j=0;j<3;j++)
	assert( fabs(m(i,j) - (j+3*i+3.14)) < 1e-8 );

    M mt = m.transpose();
    std::cout << mt << std::endl;
    assert(mt.rows() == 3);
    assert(mt.cols() == 2);
    for(int i=0;i<2;i++)
      for(int j=0;j<3;j++)
	assert(mt(j,i) == m(i,j));

    {
      HDF5writer wr("test.hdf5");
      write(wr,m,"m");
    }    
    {
      HDF5reader rd("test.hdf5");
      M mr;
      read(rd,mr,"m");

      assert(mr.size() == m.size());
      for(int i=0;i<2;i++)
	for(int j=0;j<3;j++)
	  assert(mr(i,j) == m(i,j));
    }    
  }
  {
    M m1(2,3,{0.1,0.2,0.3,1.1,1.2,1.3});
    M m2 = 3.0 * m1;
    std::cout << m1 << std::endl;
    std::cout << m2 << std::endl;
    for(int i=0;i<2;i++)
      for(int j=0;j<3;j++)
	assert(fabs(m2(i,j) - 3.0 * m1(i,j))<1e-8);

    M m3 = m2 + m1;
    std::cout << m3 << std::endl;
    for(int i=0;i<2;i++)
      for(int j=0;j<3;j++)
	assert(fabs(m3(i,j) - 4.0 * m1(i,j))<1e-8);
    
    M m4 = m3 - m2;
    std::cout << m4 << std::endl;
    for(int i=0;i<2;i++)
      for(int j=0;j<3;j++)
	assert(fabs(m4(i,j) -  m1(i,j))<1e-8);

    M m5(3,2,{7.1,6.3,
	  1.2,2.3,
	  -9.1, -1.3});
    
    M m6 = m1*m5;
    assert(m6.rows() == 2);
    assert(m6.cols() == 2);
    for(int i=0;i<2;i++)
      for(int k=0;k<2;k++){
	double v = 0.;	
	for(int j=0;j<3;j++)
	  v += m1(i,j)*m5(j,k);
	assert(fabs(m6(i,k) - v) < 1e-8);
      }

    std::cout << m6 << std::endl;
  }
  {
    M m1(2,3,{0.1,0.2,0.3,1.1,1.2,1.3});
    V v1({-4.1,3.7,0.1});
    V v2 = m1*v1;
    std::cout << v2 << std::endl;

    assert(v2.size() == 2);
    for(int i=0;i<2;i++){
      double v = 0.;
      for(int j=0;j<3;j++)
	v += m1(i,j)*v1(j);
      assert(fabs(v2(i) - v)<1e-8);
    }
  }
  {
    M m1(2,3,{0.1,0.2,0.3,1.1,1.2,1.3});
    SM m2({-1.,-2.,-3.,
	  -4.,-5.,-6.,
	  1.,2.,3.});
    
    M m3 = m1*m2;
    std::cout << m3 << std::endl;

    assert(m3.rows()==2);
    assert(m3.cols()==3);
    
    for(int i=0;i<2;i++){
      for(int k=0;k<3;k++){
	double v = 0.;
	for(int j=0;j<3;j++)
	  v += m1(i,j)*m2(j,k);	
	assert(fabs(m3(i,k) - v)<1e-8);
      }
    }
  }


  return 0;
}
