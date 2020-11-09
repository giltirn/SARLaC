#include <distribution/jackknife.h>
#include <distribution/boot_jackknife.h>
#include <distribution/utils.h>
using namespace CPSfit;

int main(void){
  {
    //test isComplex
    jackknifeDistribution<std::complex<double> > a(10);
    for(int i=0;i<10;i++) a.sample(i) = std::complex<double>(i+1,1e-6);
    assert(!isComplex(a,1e-5));
    assert(isComplex(a,1e-7));

    for(int i=0;i<10;i++) a.sample(i) = std::complex<double>(i+1,i+1);
    assert(isComplex(a));
  }

  { //test real/imag
    jackknifeDistribution<std::complex<double> > a(10);
    for(int i=0;i<10;i++) a.sample(i) = std::complex<double>(i+1,2*(i+1));
    
    jackknifeDistribution<double> b = real(a);
    for(int i=0;i<10;i++) assert(b.sample(i) == double(i+1));
       
    jackknifeDistribution<double> c = imag(a);
    for(int i=0;i<10;i++) assert(c.sample(i) == double(2*(i+1)));
  }
  { //test real/imag bootjack
    bootJackknifeInitType init(10,10);

    bootJackknifeDistribution<std::complex<double> > a(init);
    for(int i=0;i<10;i++)
      for(int j=0;j<10;j++)
	a.sample(i).sample(j) = std::complex<double>(i+j+1,2*(i+j+1));
    
    bootJackknifeDistribution<double> b = real(a);
    for(int i=0;i<10;i++) 
      for(int j=0;j<10;j++)
	assert(b.sample(i).sample(j) == double(i+j+1));

    bootJackknifeDistribution<double> c = imag(a);
    for(int i=0;i<10;i++) 
      for(int j=0;j<10;j++)
	assert(c.sample(i).sample(j) == double(2*(i+j+1)));
  }


  std::cout << "All tests passed" << std::endl;
  return 0;
}
