#include <fstream>
#include <sstream>

#include<random.h>
#include<distribution.h>
#include<utils/utils/time.h>

using namespace SARLaC;

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

template<typename T, template<typename> class V>
inline void gaussianRandom(doubleJackknifeDistribution<T,V> &v, const double mean, const double stddev){ 
  for(int i=0;i<v.size();i++) gaussianRandom(v.sample(i),mean,stddev);
}
template<typename T, template<typename> class V>
inline void gaussianRandom(blockDoubleJackknifeDistribution<T,V> &v, const double mean, const double stddev){ 
  for(int i=0;i<v.size();i++) gaussianRandom(v.sample(i),mean,stddev);
}
template<typename T, template<typename> class V>
inline void gaussianRandom(bootJackknifeDistribution<T,V> &v, const double mean, const double stddev){ 
  typedef iterate<bootJackknifeDistribution<T,V> > iter;
  for(int i=0;i<iter::size(v);i++) gaussianRandom(iter::at(i,v),mean,stddev);
}



//Number of elements in loop, not necessarily equal to nsample
template< template<typename, template<typename> class> class DistributionType >
struct num_elements{
  static inline size_t value(const int nsample){
    return nsample;
  }
};
template< >
struct num_elements<doubleJackknifeDistribution>{
  static inline size_t value(const int nsample){
    return nsample * (nsample - 1);
  }
};
template< >
struct num_elements<blockDoubleJackknifeDistribution>{
  static inline size_t value(const int nsample){
    return nsample * (nsample - 1);
  }
};
template< >
struct num_elements<bootJackknifeDistribution>{
  static inline size_t value(const int nsample){
    return nsample * (nsample + 1); //nsample+1 jackknifes
  }
};


template<typename T>
struct _setupDistribution{
  static inline T doit(const int nsample){
    return T(nsample);
  }
};
template<typename T, template<typename> class V>
struct _setupDistribution< blockDoubleJackknifeDistribution<T,V> >{
  static inline blockDoubleJackknifeDistribution<T,V> doit(const int nsample){
    return blockDoubleJackknifeDistribution<T,V>(nsample, 1); //bin size 1
  }
};
template<typename T, template<typename> class V>
struct _setupDistribution< bootJackknifeDistribution<T,V> >{
  static inline bootJackknifeDistribution<T,V> doit(const int nsample){
    return bootJackknifeDistribution<T,V>(bootJackknifeInitType(nsample,nsample,nsample)); //nsample boots, nsample samples
  }
};

template<typename DistributionType, typename BaseType, typename DistributionOp, typename BaseOp>
void benchmarkOp(const DistributionOp &dop, const BaseOp &bop, const int FLOPS, const std::string &op_descr,
		 const std::string &dist_descr, const int ntest,
		 const std::vector<DistributionType> &dists, const std::vector<std::vector<BaseType> > &vects){
  int nelem = vects[0].size();

  {
    std::vector<BaseType> out(nelem);
#pragma omp parallel for
    for(int s=0;s<nelem;s++)
      out[s] = bop(vects,s); //touch cache
    
    timer time;
    {
      time.start();
      for(size_t i=0;i<ntest;i++){
#pragma omp parallel for
	for(int s=0;s<nelem;s++)
	  out[s] = bop(vects,s);
      }
      time.stop();
    }

    double Gflops = double(ntest) * double(nelem) * double(FLOPS) / time.elapsed();
    std::cout << "std::vector " << op_descr << " with nelem=" << nelem << " and ntest=" << ntest << " : " << Gflops << " Gflops\n";
  }
  {
    DistributionType out = dists[0];
    out = dop(dists);

    timer time;
    {
      time.start();
      for(size_t i=0;i<ntest;i++){
	out = dop(dists);
      }
      time.stop();
    }
    double Gflops = double(ntest) * double(nelem) * double(FLOPS) / time.elapsed();
    std::cout << dist_descr << " " << op_descr << " with nelem=" << nelem << " and ntest=" << ntest << " : " << Gflops << " Gflops\n";
  }
   
}






template< template<typename, template<typename> class> class DistributionType, typename T, template<typename> class V= basic_vector >
void benchmarkDistribution(const size_t nsample, const size_t ntest, const std::string &descr){
  std::cout << "+++++++++++++++++" << std::endl;
  std::cout << "Benchmarking " << descr << std::endl;
  std::cout << "+++++++++++++++++" << std::endl;

  typedef DistributionType<T,V> D;
  
  int nelem = iterate<D>::size(_setupDistribution<D>::doit(nsample));

  int noperands = 6;
  std::vector<D> dists(noperands);
  std::vector<std::vector<T> > vects(noperands);
  for(int p=0;p<noperands;p++){
    dists[p] = _setupDistribution<D>::doit(nsample);
    gaussianRandom(dists[p],-10.,10.);
    
    vects[p].resize(nelem);
    gaussianRandom(vects[p],-10.,10.);
  }

  benchmarkOp([&](const std::vector<D> &d){ return d[0] * d[1]; },
	      [&](const std::vector<std::vector<T> > &v, const int s){ return v[0][s] * v[1][s]; },
	      mul_flops<T>::value, "a*b",
	      descr, ntest, dists, vects);

  benchmarkOp([&](const std::vector<D> &d){ return d[0] * d[1] + d[2] * d[3]; },
	      [&](const std::vector<std::vector<T> > &v, const int s){ return v[0][s] * v[1][s] + v[2][s] * v[3][s]; },
	      2*mul_flops<T>::value + add_flops<T>::value, "a*b + c*d",
	      descr, ntest, dists, vects);


}


int main(const int argc, const char* argv[]){

  RNG.initialize(1234);
  size_t nsample = 216;
  size_t ntest = 50000;
  
  int nthr = 1;
  for(int i=1;i<argc;i++){
    std::string sargv = argv[i];
    if(sargv == "-nthread"){
      nthr = strToAny<int>(argv[i+1]);
    }else if(sargv == "-ntest"){
      ntest = strToAny<int>(argv[i+1]);
    }else if(sargv == "-nsample"){
      nsample = strToAny<int>(argv[i+1]);
    }
  }

  omp_set_num_threads(nthr);
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "Number of threads " << omp_get_max_threads() << std::endl;
  std::cout << "Number of tests " << ntest << std::endl;
  std::cout << "Number of samples " << nsample << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  
  benchmarkDistribution<rawDataDistribution, double>(nsample, ntest, "raw data distribution (double)");
  benchmarkDistribution<rawDataDistribution, std::complex<double> >(nsample, ntest, "raw data distribution (complex)");    
  benchmarkDistribution<jackknifeDistribution, double>(nsample, ntest, "jackknife distribution (double)");
  benchmarkDistribution<doubleJackknifeDistribution, double>(nsample, ntest, "double jackknife distribution");
  benchmarkDistribution<blockDoubleJackknifeDistribution, double>(nsample, ntest, "block double jackknife distribution");
  benchmarkDistribution<bootJackknifeDistribution, double>(nsample, ntest, "boot jackknife distribution");

  return 0;
};

