#ifndef _DATA_SERIES_H
#define _DATA_SERIES_H

#include<algorithm>
#include<ostream>

#include<config.h>
#include<utils/template_wizardry.h>
#include<serialize/hdf5_serialize.h>
#include<distribution/distribution_iterate.h>

CPSFIT_START_NAMESPACE

//A class containing a series of data, for example a time series
template<typename _GeneralizedCoordinate, typename _DataType, template<typename,typename> class _PairType = std::pair>
class dataSeries{
public:
  typedef _GeneralizedCoordinate GeneralizedCoordinate;
  typedef _DataType DataType;
  typedef _PairType<GeneralizedCoordinate, DataType> ElementType;
private:
  std::vector<ElementType> series;

#ifdef HAVE_HDF5
  template<typename C, typename D, template<typename,typename> class P>
  friend void write(HDF5writer &writer, const dataSeries<C,D,P> &value, const std::string &tag);
  template<typename C, typename D, template<typename,typename> class P>
  friend void read(HDF5reader &reader, dataSeries<C,D,P> &value, const std::string &tag);
#endif
  
public:
  dataSeries(){}
  explicit dataSeries(const int n): series(n){}
  dataSeries(const int n, const int samples): series(n, ElementType(GeneralizedCoordinate(), DataType(samples)) ){}; 

  inline void resize(const int n){ series.resize(n); }
  inline void resize(const int n, const ElementType &init){ series.resize(n,init); }

  inline void clear(){ series.clear(); }
  
  //Resize where the elements are populated according to an initializer object
  //Initializer must have   ElementType operator()(const int i) const   where 0<=i<n
  template<typename Initializer>
  inline void resize(const int n, const Initializer &initializer){
    series.resize(n);
    for(int i=0;i<n;i++) series[i] = initializer(i);
  }
  //Constructor version of the above
  template<typename Initializer>
  inline dataSeries(const int n, const Initializer &initializer): series(n){
    for(int i=0;i<n;i++) series[i] = initializer(i);
  }
  
  inline int size() const{ return series.size(); }
  
  inline DataType &value(const int i){ return series[i].second; }
  inline const DataType &value(const int i) const{ return series[i].second; }

  inline GeneralizedCoordinate &coord(const int i){ return series[i].first; }
  inline const GeneralizedCoordinate &coord(const int i) const{ return series[i].first; }

  inline const ElementType &operator[](const int i) const{ return series[i]; }
  inline ElementType &operator[](const int i){ return series[i]; }

  inline void push_back(const ElementType &e){ series.push_back(e); }
  inline void push_back(const GeneralizedCoordinate &c, const DataType &d){ series.push_back(ElementType(c,d)); }

  inline void reverse(){ std::reverse(series.begin(),series.end()); }
};

template<typename _GeneralizedCoordinate, typename _DataType, template<typename,typename> class _PairType> 
std::ostream & operator<<(std::ostream & stream, const dataSeries<_GeneralizedCoordinate,_DataType,_PairType> &series){
  stream << "(";
  for(int i=0;i<series.size();i++)
    stream << "{" << series.coord(i) << " " << series.value(i) << "}" << (i != series.size()-1 ? " " : ")");
  return stream;
}
#ifdef HAVE_HDF5

template<typename C, typename D, template<typename,typename> class P>
void write(HDF5writer &writer, const dataSeries<C,D,P> &value, const std::string &tag){
  writer.enter(tag);
  write(writer,value.series,"series");
  writer.leave();
}
template<typename C, typename D, template<typename,typename> class P>
void read(HDF5reader &reader, dataSeries<C,D,P> &value, const std::string &tag){
  reader.enter(tag);
  read(reader,value.series,"series");
  reader.leave();
}

#endif

CPSFIT_END_NAMESPACE

#endif
