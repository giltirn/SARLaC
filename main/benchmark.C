#include <fstream>
#include <minimizer.h>
#include <cost_function.h>
#include <fitfunc.h>
#include <fitfunc_mapping.h>
#include <random.h>
#include <plot.h>
#include <distribution.h>
#include <data_series.h>
#include <sstream>
#include <boost/timer/timer.hpp>



struct performance{
  boost::timer::cpu_timer timer;

  double Gflops(const double FLOPs) const{
    auto dt = timer.elapsed();
    double time_ns = dt.wall;
    
    return FLOPs / time_ns;
  }
};




template<typename T>
struct mul_flops{
};

template<>
struct mul_flops<double>{
  enum {value = 1};
};
template<>
struct mul_flops<float>{
  enum {value = 1};
};
template<typename T>
struct mul_flops<std::complex<T> >{
  enum {value = 6};
};


template<typename T>
struct add_flops{
};

template<>
struct add_flops<double>{
  enum {value = 1};
};
template<>
struct add_flops<float>{
  enum {value = 1};
};
template<typename T>
struct add_flops<std::complex<T> >{
  enum {value = 2};
};




template< template<typename> class DistributionType, typename T >
void benchmarkDistribution(const size_t nsample, const size_t ntest){
  DistributionType<T> a(nsample);
  DistributionType<T> b(nsample);
  DistributionType<T> c(nsample);
  DistributionType<T> d(nsample);
  gaussianRandom(a,0.5,1.0);
  gaussianRandom(b,3,5.0);
  gaussianRandom(c,-1,9.0);
  gaussianRandom(d,3.14,0.12);

  std::vector<T> av(nsample);
  std::vector<T> bv(nsample);
  gaussianRandom(av, -133.1, 2.1);
  gaussianRandom(bv, 44.5, 12.1);
  {
    performance perf;
    std::vector<T> out(nsample);
    for(size_t i=0;i<ntest;i++){	  
      for(int s=0;s<nsample;s++)
	out[s] = av[s] * bv[s];
    }
    double Gflops = perf.Gflops(ntest * nsample * double(mul_flops<T>::value));
    std::cout << "No threading " << printType<decltype(av)>() << " * " << printType<decltype(bv)>() << " with nsample=" << nsample << " and ntest=" << ntest << " : " << Gflops << " Gflops\n";
  }
  {
    performance perf;
    std::vector<T> out(nsample);
    for(size_t i=0;i<ntest;i++){
#pragma omp parallel for
      for(int s=0;s<nsample;s++)
	out[s] = av[s] * bv[s];
    }
    double Gflops = perf.Gflops(ntest * nsample * double(mul_flops<T>::value));
    std::cout << "Threaded " << printType<decltype(av)>() << " * " << printType<decltype(bv)>() << " with nsample=" << nsample << " and ntest=" << ntest  <<" : " << Gflops << " Gflops\n";
  }
  {
    performance perf;
    DistributionType<T> out(nsample);
    for(size_t i=0;i<ntest;i++){
      out = a * b;
    }
    double Gflops = perf.Gflops(ntest * nsample * double(mul_flops<T>::value));
    std::cout << printType<decltype(a)>() << " * " << printType<decltype(b)>() << " with nsample=" << nsample << " and ntest=" << ntest << " : " << Gflops << " Gflops\n";
  }
  {
    performance perf;
    DistributionType<T> out(nsample);
       
    for(size_t i=0;i<ntest;i++){
      out = a * b + b;
    }
    double Gflops = perf.Gflops(ntest * nsample * double( mul_flops<T>::value + add_flops<T>::value ) );      
    std::cout << printType<decltype(a)>() << " * " << printType<decltype(b)>() << " + " << printType<decltype(b)>() << " with nsample=" << nsample << " and ntest=" << ntest << " : " << Gflops << " Gflops\n";
  }  
}


int main(const int argc, const char* argv[]){

  RNG.initialize(1234);
  size_t nsample = 216;
  size_t ntest = 100000;

  int thr[5] = {1,2,4,8,16};
  for(int t=0;t<5;t++){  
    omp_set_num_threads(thr[t]);
    std::cout << "Number of threads " << omp_get_max_threads() << std::endl;

    benchmarkDistribution<distribution, double>(nsample, ntest);
    benchmarkDistribution<distribution, std::complex<double> >(nsample, ntest);

      
    

  }
  return 0;
};
