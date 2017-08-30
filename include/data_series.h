#ifndef _DATA_SERIES_H
#define _DATA_SERIES_H

#include<ostream>
#include<template_wizardry.h>

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




//Generic filter interface for a data series, allowing for example restriction of data for fitting
//Prepend const to SeriesType if const access required
template<typename SeriesType>
class filteredDataSeries{
public:
  typedef typename SeriesType::GeneralizedCoordinate GeneralizedCoordinate;
  typedef typename SeriesType::DataType DataType;
  std::vector<int> subset_map;
  SeriesType &series;
public:
  template<typename Filter> 
  filteredDataSeries(SeriesType &_series,const Filter &filter):series(_series){
    for(int i=0;i<series.size();i++){
      if(filter.accept(series.coord(i), series.value(i)))
	subset_map.push_back(i);      
    }
  }
  inline int size() const{ return subset_map.size(); }
  inline typename add_const_if<DataType, SeriesType>::type & value(const int i) const{ return series.value(subset_map[i]); }
  inline typename add_const_if<GeneralizedCoordinate, SeriesType>::type & coord(const int i) const{ return series.coord(subset_map[i]); }
};

template<typename GeneralizedCoordinate>
class filterXrange{
  //accept if   x_min<=x<=x_max  
  GeneralizedCoordinate x_min;
  GeneralizedCoordinate x_max;
  bool _invert;
public:
  filterXrange(const GeneralizedCoordinate &_x_min, const GeneralizedCoordinate &_x_max): x_min(_x_min), x_max(_x_max), _invert(false){}

  void invert(){ _invert = !_invert; } //switch accept inside/outside of range
  
  template<typename T>
  bool accept(const GeneralizedCoordinate &x, const T &y) const{
    bool cond = x >= x_min && x <= x_max;
    return _invert ? !cond : cond;    
  }
};


//An accessor class that allows access to individual samples of a series of distributions with different coordinates in some general array
//This is used because the fitter expects single values for each coordinate
//Prepend const to seriesType if const access required
template<typename SeriesType>
class sampleSeries{
public:
  typedef typename SeriesType::DataType ElementType; 
  typedef typename SeriesType::GeneralizedCoordinate GeneralizedCoordinate;
  typedef typename ElementType::DataType DataType;
private:
  const int sample;
  SeriesType &series;
public:
  sampleSeries(SeriesType &_series, const int _sample):series(_series), sample(_sample){ }
  inline int size() const{ return series.size(); }
  inline typename add_const_if<DataType, SeriesType>::type & value(const int i) const{ return series.value(i).sample(sample); }
  inline typename add_const_if<GeneralizedCoordinate, SeriesType>::type & coord(const int i) const{ return series.coord(i); }
};



template<typename ResampledDataSeriesType, typename RawDataSeriesType>
inline void resample(ResampledDataSeriesType &out, const RawDataSeriesType &in){
  out.resize(in.size());
  for(int i=0;i<in.size();i++){
    out.coord(i) = in.coord(i);
    out.value(i).resample(in.value(i));    
  }  
}

#endif
