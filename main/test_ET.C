#include <fstream>
#include <minimizer.h>
#include <cost_function.h>
#include <fitfunc.h>
#include <fitfunc_mapping.h>
#include <random.h>
#include <plot.h>
#include <distribution.h>
#include <data_series.h>


int main(void){
  RNG.initialize(1234);
  int nsample = 4;

  {
    distribution<double> a(nsample, 9.0);
    distribution<double> b = sqrt(a);
    for(int i=0;i<nsample;i++){
      printf("sqrt(%f) = %f\n",a.sample(i),b.sample(i)); fflush(stdout);
      assert(b.sample(i) == 3.);
    }
  }
  {
    distribution<double> a(nsample, 1.0);
    distribution<double> b(nsample, 2.0);
    distribution<double> c(nsample, 3.0);
    
    distribution<double> d = a + b * c;
    for(int i=0;i<nsample;i++){
      printf("%f + %f * %f = %f\n",a.sample(i),b.sample(i),c.sample(i),d.sample(i)); fflush(stdout);
      assert(d.sample(i) == 7.);
    }
  }

  std::cout << "Normal exit\n"; std::cout.flush();
  return 0;
}

