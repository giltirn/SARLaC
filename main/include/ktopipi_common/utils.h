#ifndef _FIT_KTOPIPI_GPARITY_UTILS_H___
#define _FIT_KTOPIPI_GPARITY_UTILS_H___

#include<config.h>
#include<utils/macros.h>

#include<distribution.h>
#include<tensors.h>
#include<common.h>
#include<serialize.h>

CPSFIT_START_NAMESPACE


template<typename T, int Size, int Begin>
class IndexedContainer{
  std::array<T,Size> v;
public:
  GENERATE_HDF5_SERIALIZE_METHOD( (v) );

  T & operator()(const int i){ return v[i-Begin]; }
  const T & operator()(const int i) const{ return v[i-Begin]; }
};
template<typename T, int Size, int Begin>
inline void write(HDF5writer &writer, const IndexedContainer<T,Size,Begin> &d, const std::string &tag){ d.write(writer,tag); }
template<typename T, int Size, int Begin>
inline void read(HDF5reader &reader, IndexedContainer<T,Size,Begin> &d, const std::string &tag){ d.read(reader,tag); }


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

template<typename DistributionType, typename Accessor>
void average(DistributionType & into, const Accessor &data, const int size){
  //boost::timer::auto_cpu_timer total_time("average(DistributionType & into, const Accessor &data, const int size) %w s\n");
  assert(size > 0);
  into = data(0);
  for(int i=1;i<size;i++) into = into+data(i);
  into = into/double(size);
}

//On the fly bin-resample data and average
template<typename DistributionType, typename Accessor, typename BinResampler>
void binResampleAverage(DistributionType & into, const BinResampler &bin_resampler, const Accessor &data, const int size){
  assert(size > 0);

  bin_resampler.binResample(into, data(0));
  DistributionType tmp;
  for(int i=1;i<size;i++){
    bin_resampler.binResample(tmp,data(i));
    into = into+tmp;
  }
  into = into/double(size);
}

CPSFIT_END_NAMESPACE

#endif
