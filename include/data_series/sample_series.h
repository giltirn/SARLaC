#ifndef _SAMPLE_SERIES_H_
#define _SAMPLE_SERIES_H_

//An accessor class that allows access to individual samples of a series of distributions with different coordinates in some general array
//This is used because the fitter expects single values for each coordinate
//Prepend const to seriesType if const access required

#include<utils/macros.h>
#include<utils/template_wizardry.h>
#include<distribution/distribution_iterate.h>

SARLAC_START_NAMESPACE

template<typename SeriesType, typename SFINAE = void>
class sampleSeries{};

template<typename SeriesType>
class sampleSeries<SeriesType, typename std::enable_if<hasSampleMethod<typename SeriesType::DataType>::value && !hasSampleMethod<typename SeriesType::GeneralizedCoordinate>::value , void>::type >{
public:
  typedef typename SeriesType::GeneralizedCoordinate GeneralizedCoordinate;
  typedef typename SeriesType::DataType::DataType DataType;
private:
  const int sample;
  SeriesType &series;
public:
  sampleSeries(SeriesType &_series, const int _sample):series(_series), sample(_sample){ }
  inline int size() const{ return series.size(); }
  inline typename add_const_if<DataType, SeriesType>::type & value(const int i) const{ return iterate<typename SeriesType::DataType>::at(sample, series.value(i)); }
  inline typename add_const_if<GeneralizedCoordinate, SeriesType>::type & coord(const int i) const{ return series.coord(i); }
};

template<typename SeriesType>
class sampleSeries<SeriesType, typename std::enable_if<hasSampleMethod<typename SeriesType::DataType>::value && hasSampleMethod<typename SeriesType::GeneralizedCoordinate>::value , void>::type >{
public:
  typedef typename SeriesType::GeneralizedCoordinate::DataType GeneralizedCoordinate;
  typedef typename SeriesType::DataType::DataType DataType;
private:
  const int sample;
  SeriesType &series;
public:
  sampleSeries(SeriesType &_series, const int _sample):series(_series), sample(_sample){ }
  inline int size() const{ return series.size(); }
  inline typename add_const_if<DataType, SeriesType>::type & value(const int i) const{ return iterate<typename SeriesType::DataType>::at(sample, series.value(i)); }
  inline typename add_const_if<GeneralizedCoordinate, SeriesType>::type & coord(const int i) const{ return iterate<typename SeriesType::GeneralizedCoordinate>::at(sample, series.coord(i)); }
};

template<typename SeriesType, typename D> 
std::ostream & operator<<(std::ostream & stream, const sampleSeries<SeriesType,D> &series){
  stream << "(";
  for(int i=0;i<series.size();i++)
    stream << "{" << series.coord(i) << " " << series.value(i) << "}" << (i != series.size()-1 ? " " : ")");
  return stream;
}



SARLAC_END_NAMESPACE
#endif
