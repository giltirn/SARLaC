#ifndef _FIT_KTOPIPI_GPARITY_UTILS_H___
#define _FIT_KTOPIPI_GPARITY_UTILS_H___

#include<distribution.h>

template<typename T, int Size, int Begin>
class IndexedContainer{
  std::array<T,Size> v;
public:
  T & operator()(const int i){ return v[i-Begin]; }
  const T & operator()(const int i) const{ return v[i-Begin]; }
};


//Write results at various stages for debugging against other codes
template<typename DistributionType,int N>
void writeToTextFile(const NumericTensor<DistributionType,N> &m, const std::string &filename){
  std::ofstream of(filename.c_str());
  int NN = iterate<NumericTensor<DistributionType,N> >::size(m);
  for(int i=0;i<NN;i++){
    std::vector<int> coord = iterate<NumericTensor<DistributionType,N> >::unmap(i,m);
    const DistributionType &dist = m(coord.data());
    int SS = iterate<DistributionType>::size(dist);
    for(int s=0;s<SS;s++){
      std::vector<int> sample = iterate<DistributionType>::unmap(s,dist);
      for(int cc=0;cc<coord.size();cc++) of << coord[cc] << " ";
      for(int ss=0;ss<sample.size();ss++) of << sample[ss] << " ";
      of << iterate<DistributionType>::at(s,dist) << std::endl;
    }
  }
  of.close();
}

template<int N>
inline NumericTensor<rawDataDistributionD,N> bin(const NumericTensor<rawDataDistributionD,N> &t, const int bin_size){
  return t.transform([&](int const* c, const rawDataDistributionD &e){ return e.bin(bin_size); });
}

template<typename T, int N>
struct iterate<NumericTensor<T,N> >{
  static inline int size(const NumericTensor<T,N> &from){ 
    int sz = 1;
    for(int i=0;i<N;i++) sz *= from.size(i);
    return sz;
  }
  static inline std::vector<int> unmap(int i, const NumericTensor<T,N> &from){ 
    std::vector<int> coord(N);
    for(int d=N-1;d>=0;d--){
      coord[d] = i % from.size(d);
      i /= from.size(d);
    }
    return coord;
  }    
  static inline const T& at(const int i, const NumericTensor<T,N> &from){
    std::vector<int> coord = unmap(i,from);
    return from(coord.data());
  }
  static inline T & at(const int i, NumericTensor<T,N> &from){
    std::vector<int> coord = unmap(i,from);
    return from(coord.data());
  }
};  

template<typename DistributionType, typename Accessor>
void average(DistributionType & into, const Accessor &data, const int size){
  //boost::timer::auto_cpu_timer total_time("average(DistributionType & into, const Accessor &data, const int size) %w s\n");
  assert(size > 0);
  into = data(0);
  for(int i=1;i<size;i++) into = into+data(i);
  into = into/double(size);
}


struct basic_resampler{
  template<typename DistributionType>
  inline void resample(DistributionType &out, const rawDataDistributionD &in) const{ out.resample(in); }
};

template<typename DistributionType, typename Accessor, typename Resampler>
void resampleAverage(DistributionType & into, const Resampler &resampler, const Accessor &data, const int size){
  assert(size > 0);
  resampler.resample(into, data(0));
  DistributionType tmp;
  for(int i=1;i<size;i++){
    resampler.resample(tmp,data(i));
    into = into+tmp;
  }
  into = into/double(size);
}


#endif