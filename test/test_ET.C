#include <fstream>

#include<distribution.h>
#include<random.h>
#include<tensors.h>

using namespace CPSfit;


int main(void){
  RNG.initialize(1234);
  int nsample = 4;

  {
    rawDataDistribution<double> a(nsample, 9.0);
    rawDataDistribution<double> b = sqrt(a);
    for(int i=0;i<nsample;i++){
      printf("sqrt(%f) = %f\n",a.sample(i),b.sample(i)); fflush(stdout);
      assert(b.sample(i) == 3.);
    }
  }
  {
    rawDataDistribution<double> a(nsample, 1.0);
    rawDataDistribution<double> b(nsample, 2.0);
    rawDataDistribution<double> c(nsample, 3.0);
    
    rawDataDistribution<double> d = a + b * c;
    for(int i=0;i<nsample;i++){
      printf("%f + %f * %f = %f\n",a.sample(i),b.sample(i),c.sample(i),d.sample(i)); fflush(stdout);
      assert(d.sample(i) == 7.);
    }
  }
  {
    rawDataDistribution<double> a(nsample, 1.0);
    rawDataDistribution<double> b(nsample, 2.0);
    rawDataDistribution<double> c(nsample, 3.0);
    
    rawDataDistribution<double> d = a/b -1.;
    for(int i=0;i<nsample;i++){
      printf("%f/%f - 1 = %f\n",a.sample(i),b.sample(i),d.sample(i)); fflush(stdout);
      assert(d.sample(i) == -0.5);
    }
  }

  {
    NumericTensor<jackknifeDistribution<double>,1> a({3}, [&](const int *c){ return jackknifeDistribution<double>(nsample, double(*c+1)); });
    NumericTensor<jackknifeDistribution<double>,1> b({3}, [&](const int *c){ return jackknifeDistribution<double>(nsample, double(*c+2)); });

    {
      NumericTensor<jackknifeDistribution<double>,1> c({3}, [&](const int *c){ return a(c)/b(c)-b(c); });
      
      for(int i=0;i<3;i++){
	std::cout << a(&i) << "/" << b(&i) << " - " << b(&i) << " = " << c(&i) << std::endl;
	assert(c(&i).sample(0) == a(&i).sample(0)/b(&i).sample(0)-b(&i).sample(0));
      }
    }

    
    {
      NumericTensor<jackknifeDistribution<double>,1> c({3}, [&](const int *c){ return a(c)/b(c)-1.; });
      
      for(int i=0;i<3;i++){
	std::cout << a(&i) << "/" << b(&i) << "-1. = " << c(&i) << std::endl;
	assert(c(&i).sample(0) == a(&i).sample(0)/b(&i).sample(0)-1.);
      }
    }


  }

  {
    NumericSquareMatrix<double> m1(2, 0.);
    m1(0,1) = 1.;
    m1(1,0) = 1.;
    NumericSquareMatrix<double> m2(2, 0.);
    m2(0,0) = 2.;
    m2(1,1) = 2.;

    NumericSquareMatrix<double> m3 = m1 * m2;

    std::cout << m1 << "\n*\n" << m2 << "\n=\n" << m3 << std::endl;
  }

  std::cout << "Normal exit\n"; std::cout.flush();
  return 0;
}

